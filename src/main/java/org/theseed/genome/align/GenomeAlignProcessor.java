/**
 *
 */
package org.theseed.genome.align;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.proteins.Function;
import org.theseed.proteins.FunctionMap;
import org.theseed.reports.SnipReporter;
import org.theseed.sequence.ExtendedProteinRegion;
import org.theseed.sequence.RegionList;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.clustal.ClustalPipeline;

/**
 * This command will read the genomes in a directory and output the snips.  One or more genomes will be identified as the wild
 * types.  A snip is only considered relevant if it differs from all the wild types.  The first wild type will be treated as
 * the base genome, and all the snips will be displayed relative to that genome's DNA.
 *
 * The positional parameters are the name of the input directory, the file name of the base genome, and the IDs of the other wild
 * genomes.  The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more detailed progress messages
 * -m	maximum kmer distance for a region to be aligned
 * -K	kmer size for computing distances
 * -u	maximum upstream distance for protein regions
 *
 * --format		output report format
 * --workDir	working directory name; the default is "Temp" in the current directory
 *
 *
 * @author Bruce Parrello
 *
 */
public class GenomeAlignProcessor extends BaseAlignProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GenomeAlignProcessor.class);
    /** function definitions */
    private FunctionMap funMap;
    /** list of regions to align by base genome feature */
    private Map<Feature, RegionList> alignMap;
    /** reporting object */
    private SnipReporter reporter;
    /** ID of the base genome */
    private String baseId;

    // COMMAND-LINE OPTIONS

    /** output format */
    @Option(name = "--format", usage = "output format")
    private SnipReporter.Type outputFormat;

    /** maximum upstream distance for protein neighborhoods */
    @Option(name = "-u", aliases = { "--upstream" }, metaVar = "1000", usage = "maximum upstream distance for protein neighborhoods")
    private int maxUpstream;

    /** input genome directory */
    @Argument(index = 0, metaVar = "inDir", usage = "input genome directory", required = true)
    private File inDir;

    /** base genome ID */
    @Argument(index = 1, metaVar = "baseGtoFile", usage = "base genome file", required = true)
    private File baseGto;

    /** other wild-type genome IDs */
    @Argument(index = 2, metaVar = "222222.2 333333.3 ...", usage = "other wild-type genome IDs for filtering", multiValued = true)
    private List<String> altIds;

    @Override
    protected void setProcessDefaults() {
        this.outputFormat = SnipReporter.Type.TEXT;
    }

    @Override
    protected void validateProcessParms() throws IOException {
        // Verify the upstream distance.
        if (this.maxUpstream < 0)
            throw new IllegalArgumentException("Upstream distance must be 0 or more.");
        // Verify the input directory.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        // Verify the base genome.
        if (! this.baseGto.canRead())
            throw new FileNotFoundException("Base genome file " + this.baseGto + " not found or unreadable.");
        // Create the report object.
        this.reporter = this.outputFormat.create(System.out);
    }

    @Override
    protected void runCommand() throws Exception {
        // Create the function map.
        this.funMap = new FunctionMap();
        // Create the region list map.  It is keyed by feature with sorting by location, so the output is in chromosome order.
        this.alignMap = new TreeMap<Feature, RegionList>(new Feature.LocationComparator());
        // Build all the alignment lists.
        this.buildAlignments();
        // Denote we are done registering genomes.
        this.reporter.initializeOutput();
        // Create our temporary file.
        File tempFile = File.createTempFile("align", ".fa", this.getWorkDir());
        tempFile.deleteOnExit();
        // Loop through the alignments.
        log.info("Processing alignments.");
        for (Map.Entry<Feature, RegionList> alignEntry : this.alignMap.entrySet()) {
            RegionList regions = alignEntry.getValue();
            if (regions.size() >= 3) {
                // Here we have enough data to do an alignment.
                Feature feat = alignEntry.getKey();
                log.info("Performing alignment on {}.", feat);
                regions.save(tempFile);
                ClustalPipeline aligner = new ClustalPipeline(tempFile);
                List<Sequence> alignment = aligner.run();
                // Output the alignment.
                this.reporter.processAlignment(feat.getFunction(), regions, alignment);
            }
        }
        // Finish the report.
        this.reporter.close();
    }

    /**
     * This method reads in all the genomes and builds the alignment lists.
     * @throws IOException
     */
    private void buildAlignments() throws IOException {
        // Read in the base genome and sort the regions by function.  When we process the other genomes, we will use the
        // function-to-region map to find the best feature for alignment.
        Genome base = new Genome(this.baseGto);
        log.info("Scanning base genome {}.", base);
        this.baseId = base.getId();
        Map<String, RegionList> baseMap = RegionList.createMap(this.funMap, base, this.maxUpstream);
        // Now prime the alignment lists from the base map.
        log.info("Sorting features by location.");
        for (RegionList regions : baseMap.values()) {
            for (ExtendedProteinRegion region : regions) {
                RegionList singleton = new RegionList();
                singleton.add(region);
                this.alignMap.put(region.getFeature(), singleton);
            }
        }
        log.info("{} features with {} functions processed for base genome.", this.alignMap.size(), baseMap.size());
        // Register the base genome for the report.
        this.reporter.register(base);
        // Now read the other genomes.  If we find the base, we skip it.
        log.info("Scanning input directory {}.", this.inDir);
        GenomeDirectory genomes = new GenomeDirectory(this.inDir);
        for (Genome genome : genomes) {
            if (this.baseId.contentEquals(genome.getId()))
                log.info("Base genome found in input directory-- skipped.");
            else {
                log.info("Processing input genome {}.", genome);
                // Register the genome if it is NOT one of the wild strains.  Everything is aligned, but only the base and the
                // non-wilds are output.
                if (! this.altIds.contains(genome.getId()))
                    this.reporter.register(genome);
                // Process the regions for this genome.
                int foundCount = 0;
                int badFunCount = 0;
                int tooFarCount = 0;
                ExtendedProteinRegion.GenomeIterator iter = new ExtendedProteinRegion.GenomeIterator(genome, this.maxUpstream);
                while (iter.hasNext()) {
                    // Get this region and try to find the function in the function map.
                    ExtendedProteinRegion region = iter.next();
                    Feature feat = region.getFeature();
                    Function fun = this.funMap.getByName(feat.getFunction());
                    if (fun == null)
                        badFunCount++;
                    else {
                        // Here we found it.  If we didn't find it, then it won't match anything anyway.
                        RegionList regions = baseMap.get(fun.getId());
                        ExtendedProteinRegion closest = regions.getClosest(region, this.getMaxDist());
                        if (closest == null)
                            tooFarCount++;
                        else {
                            // Here we have an eligible close feature.  Add this region to its alignment list.
                            feat = closest.getFeature();
                            this.alignMap.get(feat).add(region);
                            foundCount++;
                        }
                    }
                }
                log.info("{} regions queued for alignment.  {} had unusual functions, {} were too far to align.",
                        foundCount, badFunCount, tooFarCount);
            }
        }
    }

}
