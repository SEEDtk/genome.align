/**
 *
 */
package org.theseed.genome.align;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.proteins.Function;
import org.theseed.proteins.FunctionMap;
import org.theseed.reports.MultiAlignReporter;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.clustal.ClustalPipeline;

/**
 * This command creates alignments for a list of genomes.  The DNA sequences will be organized by functional assignment and then
 * aligned.  Only functions that occur at least three times will participate in the alignment, and only protein-encoding genes
 * will be considered.  At least one sequence must be in the first genome, as well.
 *
 * The positional parameters are the name of the file containing the first genome, then the names of the files containing the other
 * genomes.  The alignments will be written to the specified output file.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	show more detailed progress messages
 * -m	maximum kmer distance for a feature to be aligned
 * -K	kmer size for computing distances
 * -o	output file or directory
 *
 * --format		output report format
 * --workDir	working directory name; the default is "Temp" in the current directory
 * --alt		ID of a genome other than the base that is to be used as an alternate base; a snip is only output if it does not
 * 				match the base or any of the alternates
 * --upstream	check for illusory indels
 *
 * @author Bruce Parrello
 *
 */
public class GtoAlignProcessor extends BaseAlignProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GtoAlignProcessor.class);
    /** role map of functional assignments */
    private FunctionMap functionMap;
    /** hash of function IDs to sequence lists */
    private Map<String, SequenceList> sequenceMap;
    /** base genome */
    private Genome baseGenome;
    /** map of genome IDs to other genomes */
    private Map<String, Genome> genomeMap;
    /** number of upstream region recaptures */
    private int recaptures;

    // COMMAND-LINE OPTIONS

    /** output format */
    @Option(name = "--format", usage = "output format")
    private MultiAlignReporter.Type outputFormat;

    /** output file */
    @Option(name = "--output", aliases = { "-o" }, metaVar = "outFile.html", usage = "output file (or directory, if appropriate)",
            required = true)
    private File outFile;

    /** if specified, the upstream will be checked for illusory indels */
    @Option(name = "--upstream", usage = "if specified, an upstream check will be made for an illusory indel")
    private boolean upstreamCheck;

    /** alternate base genome IDs */
    @Option(name = "--alt", metaVar = "83333.1", usage = "alternate base genome ID")
    private String[] altBases;

    /** first genome */
    @Argument(index = 0, metaVar = "baseGTO", usage = "primary genome file", required = true)
    private File gtoBaseFile;

    /** comparative genomes */
    @Argument(index = 1, metaVar = "gto2 gto3 ...", usage = "genomes to align to primary", multiValued = true, required = true)
    private List<File> gtoFiles;

    protected void setProcessDefaults() {
        this.outputFormat = MultiAlignReporter.Type.TEXT;
        this.altBases = new String[0];
        this.upstreamCheck = false;
    }

    protected void validateProcessParms() throws IOException {
        // Verify the genome files.
        if (! this.gtoBaseFile.canRead())
            throw new FileNotFoundException("Genome file " + this.gtoBaseFile + " not found or unreadable.");
        for (File gtoFile : this.gtoFiles) {
            if (! gtoFile.canRead())
                throw new FileNotFoundException("Genome file " + gtoFile + " not found or unreadable.");
        }
    }

    @Override
    protected void runCommand() throws Exception {
        // Create the temporary FASTA file.
        File tempFile = File.createTempFile("align", ".fna", this.getWorkDir());
        tempFile.deleteOnExit();
        // Create the genome map.
        this.genomeMap = new HashMap<String, Genome>(this.gtoFiles.size());
        // Create the function maps.
        log.info("Initializing function maps.");
        this.functionMap = new FunctionMap();
        this.sequenceMap = new HashMap<String, SequenceList>(3000);
        this.recaptures = 0;
        try (MultiAlignReporter reporter = this.outputFormat.create(this.outFile)) {
            // Process the base genome to compute the functions of interest.
            this.processBase(this.gtoBaseFile);
            // Initialize the output report.
            reporter.openReport(this.baseGenome, this.altBases);
            // Loop through the other genomes.
            for (File gtoFile : this.gtoFiles) {
                Genome genome = new Genome(gtoFile);
                reporter.registerGenome(genome);
                log.info("Scanning genome {}.", genome);
                if (this.upstreamCheck)
                    this.genomeMap.put(genome.getId(), genome);
                // Loop through the genome's pegs.
                int kept = 0;
                for (Feature feat : genome.getPegs()) {
                    // Get the function and look for a sequence list.
                    String function = feat.getFunction();
                    Function fun = this.functionMap.getByName(function);
                    if (fun != null) {
                        SequenceList seqs = this.sequenceMap.get(fun.getId());
                        if (seqs != null) {
                            if (this.addFeature(seqs, feat))
                                kept++;
                        }
                    }
                }
                log.info("{} sequences kept from {}.", kept, genome);
            }
            int alignCount = 0;
            // Loop through the alignments
            for (Map.Entry<String, SequenceList> alignRequest : this.sequenceMap.entrySet()) {
                // Verify that this alignment is big enough and has variations.
                SequenceList seqs = alignRequest.getValue();
                if (seqs.size() >= 3 && seqs.getMaxDist() > 0.0) {
                    String function = this.functionMap.getName(alignRequest.getKey());
                    log.info("Processing alignment for {}.", function);
                    // Save the sequences.
                    seqs.save(tempFile);
                    // Process the alignment.
                    ClustalPipeline aligner = new ClustalPipeline(tempFile);
                    List<Sequence> alignment = aligner.run();
                    if (this.upstreamCheck)
                        this.checkUpstream(alignment);
                    reporter.writeAlignment(seqs.getBaseFid(), function, alignment);
                    alignCount++;
                }
            }
            // Finish the report.
            reporter.closeReport();
            log.info("{} alignments output.", alignCount);
            if (this.upstreamCheck)
                log.info("{} upstream regions recaptured.", this.recaptures);
        }
    }

    /**
     * Add upstream sequences to fix indels caused by differing start calls.
     *
     * @param alignment		alignment to verify
     */
    private void checkUpstream(List<Sequence> alignment) {
        // Loop through the alignment looking for an indel at the start.
        List<Sequence> fronted = new ArrayList<Sequence>(alignment.size());
        List<Sequence> full = new ArrayList<Sequence>(alignment.size());
        for (Sequence seq : alignment) {
            if (StringUtils.startsWith(seq.getSequence(), "-"))
                fronted.add(seq);
            else
                full.add(seq);
        }
        // If we have both fronted and full, we should check the upstreams.
        if (fronted.size() > 0 && full.size() > 0) {
            // Loop through the fronted sequences.
            for (Sequence frontedSeq : fronted) {
                // What we do now is find the upstream sequence that corresponds to the indels at the front.  If it corresponds with
                // the prefix of a full sequence, we plug it in.
                int indelLength = StringUtils.indexOfAnyBut(frontedSeq.getSequence(), '-');
                // Find the genome for this sequence.
                Genome genome = this.genomeMap.get(Feature.genomeOf(frontedSeq.getLabel()));
                // Compute the upstream sequence.
                Location loc = Location.fromString(frontedSeq.getComment());
                Location loc2 = loc.upstream(indelLength);
                String upstream = genome.getDna(loc2);
                if (upstream.length() == indelLength) {
                    // Here we have a valid upstream region.  Search for it in the full sequences.
                    boolean found = false;
                    for (int i = 0; i < full.size() && ! found; i++) {
                        Sequence fullSeq = full.get(i);
                        if (StringUtils.startsWith(fullSeq.getSequence(), upstream))
                            found = true;
                    }
                    if (found) {
                        // We found a match, so update the fronted sequence.
                        String tail = StringUtils.substring(frontedSeq.getSequence(), indelLength);
                        frontedSeq.setSequence(upstream + tail);
                        log.info("Upstream region added to {}.", frontedSeq.getLabel());
                        // Update the location, too.
                        int contigLen = genome.getContig(loc.getContigId()).length();
                        loc2 = loc.expandUpstream(indelLength, contigLen);
                        frontedSeq.setComment(loc2.toString());
                        // Count the recapture.
                        this.recaptures++;
                    }
                }
            }
        }
    }

    /**
     * Process the base genome.  This includes isolating all the non-hypothetical functions and storing the
     * DNA sequences in the sequence lists.
     *
     * @param gtoFile	file containing the base genome
     *
     * @throws IOException
     */
    private void processBase(File gtoFile) throws IOException {
        // This will count the total features kept.
        int kept = 0;
        // Read in the genome.
        Genome genome = new Genome(gtoFile);
        log.info("Scanning features from {}.", genome);
        this.baseGenome = genome;
        if (this.upstreamCheck)
            this.genomeMap.put(genome.getId(), genome);
        // Loop through the pegs.
        for (Feature feat : genome.getPegs()) {
            String function = feat.getFunction();
            if (! function.contentEquals("hypothetical protein")) {
                String funId = this.functionMap.findOrInsert(function).getId();
                SequenceList seqList = this.sequenceMap.get(funId);
                if (seqList != null) {
                    // Here we have an existing function, so add this feature to the sequence list.
                    if (this.addFeature(seqList, feat)) kept++;
                } else {
                    // Here we need to create a new sequence list.
                    Location loc = feat.getLocation();
                    seqList = new SequenceList(feat.getId(), loc.toString(), genome.getDna(loc));
                    this.sequenceMap.put(funId, seqList);
                    kept++;
                }
            }
        }
        log.info("{} functions found in {}, {} features kept.", this.sequenceMap.size(), genome, kept);
    }

    /**
     * Add a new feature to a sequence list.
     *
     * @param seqList	target sequence list
     * @param feat		feature whose sequence should be added
     *
     * @return TRUE if the new feature was eligible, else FALSE
     */
    private boolean addFeature(SequenceList seqList, Feature feat) {
        Location loc = feat.getLocation();
        String dna = feat.getDna();
        boolean retVal = seqList.add(feat.getId(), loc.toString(), dna, this.getMaxDist());
        return retVal;
    }

}
