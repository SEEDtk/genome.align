/**
 *
 */
package org.theseed.genome.align;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Contig;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.proteins.Function;
import org.theseed.proteins.FunctionMap;
import org.theseed.sequence.ExtendedProteinRegion;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.RegionList;
import org.theseed.sequence.Sequence;

/**
 * This command performs a BLAST to splice the DNA from a gegnome (the source) into the DNA from another genome (the
 * reference).  The presumption is that the source is a mutation of the reference, and the reference is more complete.  The
 * result will hopefully be a DNA sequence representing a more complete version of the mutant.
 *
 * The positional parameters are the name of the source genome file and the name of the reference genome file.  The combined file will be
 * written to the standard output.
 *
 * The basic strategy will be to extract all the extended protein regions from the source and use them to replace the matching extended protein
 * regions from the reference.
 *
 * A "map.tbl" file will be produced in the working directory containing the mapping between features in the genome.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	show more detailed progress messages
 * -m	maximum kmer distance for a region to be aligned
 * -K	kmer size for computing distances
 * -u	maximum upstream distance for protein regions
 *
 * --workDir	working directory name; the default is "Temp" in the current directory
 *
 * @author Bruce Parrello
 *
 */
public class SpliceProcessor extends BaseAlignProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SpliceProcessor.class);
    /** function map */
    private FunctionMap funMap;
    /** map of function IDs to region lists from the reference genome */
    private Map<String, RegionList> refRegionMap;
    /** sorted set of region pairs containing matches between the reference and source genomes */
    private SortedSet<RegionList> matchedRegions;
    /** reference genome */
    private Genome refGenome;
    /** next region in reference genome */
    private ExtendedProteinRegion refRegion;
    /** location of next region in reference genome */
    private Location refLocation;
    /** replacement region in source genome */
    private ExtendedProteinRegion sourceRegion;
    /** output file for mapping */
    private PrintWriter mapStream;

    // COMMAND-LINE OPTIONS

    /** maximum upstream distance for protein neighborhoods */
    @Option(name = "-u", aliases = { "--upstream" }, metaVar = "1000", usage = "maximum upstream distance for protein neighborhoods")
    private int maxUpstream;

    /** source genome file */
    @Argument(index = 0, metaVar = "source.gto", usage = "genome to splice into reference")
    private File sourceFile;

    /** reference genome file */
    @Argument(index = 1, metaVar = "reference.gto", usage = "reference genome into which source will be spliced")
    private File refFile;

    @Override
    protected void setProcessDefaults() {
        this.maxUpstream = 0;
    }

    @Override
    protected void validateProcessParms() throws IOException, ParseFailureException {
        if (this.maxUpstream < 0)
            throw new ParseFailureException("Upstream length must be non-negative.");
        if (! this.sourceFile.canRead())
            throw new FileNotFoundException("Source genome file " + this.sourceFile + " not found or unreadable.");
        if (! this.refFile.canRead())
            throw new FileNotFoundException("Reference genome file " + this.refFile + " not found or unreadable.");
        // Create the map output file.
        File mapFile = new File(this.getWorkDir(), "map.tbl");
        log.info("Map will be written to {}.", mapFile);
        this.mapStream = new PrintWriter(mapFile);
    }

    @Override
    protected void runCommand() throws Exception {
        // Start the map output.
        this.mapStream.println("source_fid\treference_fid\tfunction");
        // Read in the reference genome.
        this.refGenome = new Genome(this.refFile);
        // Create the function map.
        this.funMap = new FunctionMap();
        // Get the map of extended protein regions.
        log.info("Scanning proteins in {}.", this.refGenome);
        this.refRegionMap = RegionList.createMap(this.funMap, this.refGenome, this.maxUpstream);
        // Create the storage set for mappings.
        this.matchedRegions = new TreeSet<RegionList>();
        // Read in the source genome.
        Genome sourceGenome = new Genome(this.sourceFile);
        // Get all the extended protein regions for the source genome and match them to the reference genome.
        this.processSource(sourceGenome);
        // Now we want to write the reference genome contigs, plugging in the source regions for the matching
        // reference regions.  To do this, we need an iterator through the matched regions.
        Iterator<RegionList> iter = this.matchedRegions.iterator();
        if (! iter.hasNext())
            throw new IllegalArgumentException("Source genome does not match any part of reference genome.");
        // The FASTA will be written to the standard output.
        try (FastaOutputStream outStream = new FastaOutputStream(System.out)) {
            List<Contig> refContigs = this.refGenome.getContigs().stream().sorted().collect(Collectors.toList());
            // Position on the first matched region pair.
            this.getNextRegion(iter);
            // Loop through the contigs.
            for (Contig contig : refContigs) {
                log.info("Processing contig {}.", contig.getId());
                // We will assemble the new version of this contig in here.
                StringBuffer newSequence = new StringBuffer(contig.length());
                String contigId = contig.getId();
                // Denote we are at the beginning of the contig and then loop through it.
                int pos = 1;
                while (pos <= contig.length()) {
                    int compare = contigId.compareTo(this.refLocation.getContigId());
                    if (compare < 0) {
                        // The next region is in a subsequent contig, so move to the next one.
                        this.flushContig(newSequence, contig, pos, outStream);
                        pos = contig.length() + 1;
                    } else if (compare > 0)
                        throw new IllegalArgumentException("Reference region has invalid location " + this.refLocation.toString() + " for feature " +
                                this.refRegion.getFeature().getId() + ".");
                    else {
                        // The next region is in this contig.  Write out everything in front of it.
                        this.writeContig(newSequence, contig, pos, this.refLocation.getLeft());
                        // Now put in the source region.
                        this.writeSourceRegion(newSequence);
                        // Denote we're past this region.
                        pos = this.refLocation.getRight() + 1;
                        // Get the next region.
                        this.getNextRegion(iter);
                    }
                }
            }
        } finally {
            // Insure we close the map stream.
            this.mapStream.close();
        }
        log.info("Processing complete.");
    }

    /**
     * Write the current source region into the sequence buffer.
     *
     * The only tricky part here is that if we are replacing a sequence from the minus strand, we need to reverse complement.
     *
     * @param newSequence	output sequence buffer
     */
    private void writeSourceRegion(StringBuffer newSequence) {
        String dna = this.sourceRegion.getSequence();
        if (this.sourceRegion.getFullLocation().getDir() == '-')
            dna = Contig.reverse(dna);
        newSequence.append(dna);
    }

    /**
     * Write the specified region of the contig into the sequence buffer.
     *
     * @param newSequence	output sequence buffer
     * @param contig		source contig
     * @param pos			leftmost position
     * @param end			point past the rightmost position
     */
    private void writeContig(StringBuffer newSequence, Contig contig, int pos, int end) {
        String seq = contig.getSequence();
        newSequence.append(StringUtils.substring(seq, pos - 1, end - 1));
    }

    /**
     * Flush the current contig to the output.  The remainder of the contig is added to the
     * sequence buffer and then the sequence buffer is written to the output.
     *
     * @param contig		contig to flush
     * @param pos			current position in the contig
     * @param outStream 	output stream for sequence
     *
     * @throws IOException
     */
    private void flushContig(StringBuffer newSequence, Contig contig, int pos, FastaOutputStream outStream) throws IOException {
        String seq = contig.getSequence();
        if (pos <= seq.length())
            newSequence.append(StringUtils.substring(seq, pos - 1, seq.length()));
        Sequence contigSeq = new Sequence(contig.getId(), contig.getDescription(), newSequence.toString());
        outStream.write(contigSeq);
        log.info("Contig {} written.", contig.getId());
    }

    /**
     * Set up the next region substitution.
     *
     * @param iter	iterator through the region pairs
     */
    private void getNextRegion(Iterator<RegionList> iter) {
        if (! iter.hasNext()) {
            // Here we are at the end, so we need to create a dummy trailer.
            this.refLocation = Location.create("\u00FF", 1, 1);
            this.refRegion = null;
            this.sourceRegion = null;
        } else {
            RegionList pair = iter.next();
            this.refRegion = pair.get(0);
            this.refLocation = this.refRegion.getFullLocation();
            this.sourceRegion = pair.get(1);
        }
    }

    /**
     * This method runs through the source genome, matching each protein region to the closest region in
     * the reference genome.
     *
     * @param sourceGenome	source genome to process
     */
    private void processSource(Genome sourceGenome) {
        // This will count the regions placed.
        int placeCount = 0;
        int lostCount = 0;
        // Scan the source genome for regions.
        log.info("Scanning proteins in {}.", sourceGenome);
        RegionList sourceList = new RegionList(sourceGenome, this.maxUpstream);
        for (ExtendedProteinRegion region : sourceList) {
            // Find the regions with the same function as this one.
            String function = region.getFeature().getPegFunction();
            Function fun = this.funMap.getByName(function);
            if (fun == null) {
                lostCount++;
                this.mapStream.format("%s\t\t%s%n", region.getFeature().getId(), function);
            } else {
                // Look for the closest region with the same function.
                RegionList regions = this.refRegionMap.get(fun.getId());
                ExtendedProteinRegion closest = regions.getClosest(region, this.getMaxDist());
                if (closest == null) {
                    lostCount++;
                    this.mapStream.format("%s\t\t%s%n", region.getFeature().getId(), function);
                } else {
                    // Post this region as a replacement for the reference region.
                    RegionList pair = new RegionList();
                    pair.add(closest);
                    pair.add(region);
                    this.matchedRegions.add(pair);
                    placeCount++;
                    // Create this pairing's map entry.
                    this.mapStream.format("%s\t%s\t%s%n", region.getFeature().getId(), closest.getFeature().getId(), function);
                }
            }
        }
        log.info("{} regions placed, {} lost.", placeCount, lostCount);
    }

}
