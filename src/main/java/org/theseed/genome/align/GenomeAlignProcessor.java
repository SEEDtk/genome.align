/**
 *
 */
package org.theseed.genome.align;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.FeatureFilter;
import org.theseed.proteins.Function;
import org.theseed.proteins.FunctionMap;
import org.theseed.reports.HtmlSnipReporter;
import org.theseed.reports.HtmlSnipReporter.Sort;
import org.theseed.reports.SnipReporter;
import org.theseed.sequence.ExtendedProteinRegion;
import org.theseed.sequence.MarkedRegionList;
import org.theseed.sequence.RegionList;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.clustal.ClustalPipeline;
import org.theseed.utils.ParseFailureException;

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
 * -w	maximum cell character width for HTML tables
 * -o	output file (if not STDOUT)
 *
 * --filter		types of filters; features that fail the filter will not be added to the pending alignments
 * --format		output report format
 * --workDir	working directory name; the default is "Temp" in the current directory
 * --sort		sort order for HTML output
 * --gFile		file containing a list of genome IDs, used to determine the order of the output columns
 * --groups		name of a tab-delimited file containing group membership information (modulons, subsystems)
 * 				for the base genome features
 *
 * @author Bruce Parrello
 *
 */
public class GenomeAlignProcessor extends BaseAlignProcessor implements FeatureFilter.IParms, SnipReporter.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GenomeAlignProcessor.class);
    /** function definitions */
    private FunctionMap funMap;
    /** list of regions to align by base genome feature */
    private Map<Feature, MarkedRegionList> alignMap;
    /** ID of the base genome */
    private String baseId;
    /** list of filters */
    private List<FeatureFilter> filters;
    /** genome ID reordering list */
    private List<String> genomeIds;
    /** map of feature IDs to group descriptors */
    private Map<String, List<String>> groupMap;
    /** output stream */
    private OutputStream outStream;

    // COMMAND-LINE OPTIONS

    /** output format */
    @Option(name = "--format", usage = "output format")
    private SnipReporter.Type outputFormat;

    /** sort order for HTML tables */
    @Option(name = "--sort", usage = "sort order for HTML tables")
    private HtmlSnipReporter.Sort sortOrder;

    /** filtering scheme for features */
    @Option(name = "--filter", usage = "type of feature filtering")
    private List<FeatureFilter.Type> filterTypes;

    /** maximum cell width for HTML tables */
    @Option(name = "-w", aliases = { "--cellWidth" }, metaVar = "10", usage = "maximum character width for a snip display cell")
    private int cellWidth;

    /** maximum upstream distance for protein neighborhoods */
    @Option(name = "-u", aliases = { "--upstream" }, metaVar = "1000", usage = "maximum upstream distance for protein neighborhoods")
    private int maxUpstream;

    /** file of genome IDs specifying the column ordering on the output */
    @Option(name = "--gFile", metaVar = "genomeIDs.txt", usage = "file of genome IDs specifying the output column order")
    private File gFile;

    /** file of grouping information */
    @Option(name = "--groups", metaVar = "groups.tbl", usage = "file of group information by feature ID")
    private File groupFile;

    /** grouping output file */
    @Option(name = "--groupOut", metaVar = "groups.snips.tbl", usage = "output file for group snip information")
    private File groupOutFile;

    /** output file (if not STDOUT) */
    @Option(name = "-o", aliases = { "--output" }, metaVar = "outFile.html", usage = "output file (if not STDOUT)")
    private File outFile;

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
        this.filterTypes = new ArrayList<FeatureFilter.Type>();
        this.cellWidth = 20;
        this.sortOrder = HtmlSnipReporter.Sort.CHANGES;
        this.gFile = null;
        this.genomeIds = null;
        this.maxUpstream = 100;
        this.groupFile = null;
        this.outFile = null;
        this.groupOutFile = null;
    }

    @Override
    protected void validateProcessParms() throws IOException, ParseFailureException {
        // Verify the upstream distance.
        if (this.maxUpstream < 0)
            throw new ParseFailureException("Upstream distance must be 0 or more.");
        // Verify the input directory.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        // Verify the base genome.
        if (! this.baseGto.canRead())
            throw new FileNotFoundException("Base genome file " + this.baseGto + " not found or unreadable.");
        // Create the filter list.
        this.filters = new ArrayList<FeatureFilter>(this.filterTypes.size());
        for (FeatureFilter.Type filterType : this.filterTypes)
            this.filters.add(filterType.create(this));
        // Create the genome ID ordering list (if any).
        if (this.gFile != null) {
            if (! this.gFile.canRead())
                throw new FileNotFoundException("Genome ID ordering file " + this.gFile + " not found or unreadable.");
            try (TabbedLineReader gStream = new TabbedLineReader(this.gFile, 1)) {
                this.genomeIds = new ArrayList<String>();
                for (TabbedLineReader.Line line : gStream)
                    this.genomeIds.add(line.get(0));
                log.info("{} genome IDs read from ordering file {}.", this.genomeIds.size(), this.gFile);
            }
        }
        // Create the group map (if any).
        if (this.groupFile == null)
            this.groupMap = Collections.emptyMap();
        else {
            if (! this.groupFile.canRead())
                throw new FileNotFoundException("Group file " + this.groupFile + " not found or unreadable.");
            this.groupMap = new HashMap<String, List<String>>(1000);
            try (TabbedLineReader gStream = new TabbedLineReader(this.groupFile)) {
                for (TabbedLineReader.Line line : gStream) {
                    String fid = line.get(0);
                    // Here we build the description string for the feature.
                    List<String> groupList = new ArrayList<String>(10);
                    int arNum = line.getInt(2);
                    groupList.add(String.format("AR%d", arNum));
                    String mods = line.get(1);
                    if (! mods.isEmpty())
                        groupList.addAll(Arrays.asList(StringUtils.split(mods, ",")));
                    this.groupMap.put(fid, groupList);
                }
            }
            log.info("{} group records read from {}.", this.groupMap.size(), this.groupFile);
        }
        // Connect to the output file.
        if (this.outFile == null) {
            log.info("Results will be written to the standard output.");
            this.outStream = System.out;
        } else {
            log.info("Results will be written to {}.", this.outFile);
            this.outStream = new FileOutputStream(this.outFile);
        }


    }

    @Override
    protected void runCommand() throws Exception {
        try (SnipReporter reporter = this.outputFormat.create(this.outStream, this)) {
            // Create the function map.
            this.funMap = new FunctionMap();
            // Create the region list map.  It is keyed by feature with sorting by location, so the output is in chromosome order.
            this.alignMap = new TreeMap<Feature, MarkedRegionList>(new Feature.LocationComparator());
            // Set up the group output file, if any.
            if (this.groupOutFile != null)
                reporter.setupFeatureOutput(this.groupOutFile);
            // Build all the alignment lists.
            this.buildAlignments(reporter);
            // Fix the ordering (if necessary).
            if (this.genomeIds != null)
                reporter.reorder(this.genomeIds);
            // Denote we are done registering genomes.
            reporter.initializeOutput();
            // Create our temporary file.
            File tempFile = File.createTempFile("align", ".fa", this.getWorkDir());
            tempFile.deleteOnExit();
            // Loop through the alignments.
            log.info("Processing alignments.");
            for (Map.Entry<Feature, MarkedRegionList> alignEntry : this.alignMap.entrySet()) {
                MarkedRegionList regions = alignEntry.getValue();
                // We only align these regions if at least one has a functional difference.
                if (regions.getCounter() > 0) {
                    // Here we have enough data to do an alignment.
                    Feature feat = alignEntry.getKey();
                    log.info("Performing alignment on {}.", feat);
                    regions.save(tempFile);
                    ClustalPipeline aligner = new ClustalPipeline(tempFile);
                    List<Sequence> alignment = aligner.run();
                    // Output the alignment.
                    reporter.processAlignment(feat, regions, alignment);
                }
            }
            // Finish the report.
            reporter.finishReport();
        } finally {
            if (this.outFile != null)
                this.outStream.close();
        }
    }

    /**
     * This method reads in all the genomes and builds the alignment lists.
     *
     * @throws IOException
     */
    private void buildAlignments(SnipReporter reporter) throws IOException {
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
                MarkedRegionList singleton = new MarkedRegionList();
                singleton.add(region);
                this.alignMap.put(region.getFeature(), singleton);
            }
        }
        log.info("{} features with {} functions processed for base genome.", this.alignMap.size(), baseMap.size());
        // Register the base genome for the report.
        reporter.register(base);
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
                boolean wild = this.altIds.contains(genome.getId());
                if (! wild)
                    reporter.register(genome);
                // Process the regions for this genome.
                int foundCount = 0;
                int badFunCount = 0;
                int tooFarCount = 0;
                int filterCount = 0;
                int sameCount = 0;
                ExtendedProteinRegion.GenomeIterator iter = new ExtendedProteinRegion.GenomeIterator(genome, this.maxUpstream);
                while (iter.hasNext()) {
                    // Get this region and try to find the function in the function map.
                    ExtendedProteinRegion region = iter.next();
                    if (! this.checkFilters(region))
                        filterCount++;
                    else {
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
                                // Here we have an eligible close feature.  Add this region to its alignment list.  We always
                                // count a wild strain, but we skip anything that has the same protein.
                                Feature featBase = closest.getFeature();
                                String protBase = featBase.getProteinTranslation();
                                String upstreamDnaBase = closest.getUpstreamDna();
                                MarkedRegionList alignment = this.alignMap.get(featBase);
                                alignment.add(region);
                                if (! wild) {
                                    if (protBase.contentEquals(feat.getProteinTranslation()) && upstreamDnaBase.contentEquals(region.getUpstreamDna()))
                                        sameCount++;
                                    else {
                                        alignment.increment();
                                        foundCount++;
                                    }
                                }
                            }
                        }
                    }
                }
                log.info("{} regions queued for alignment.  {} had unusual functions, {} were too far to align.",
                        foundCount, badFunCount, tooFarCount);
                log.info("{} regions rejected by filtering.  {} were functionally identical to the base.", filterCount, sameCount);
            }
        }
    }

    /**
     * @return TRUE if this region is acceptable, FALSE if it is rejected by the feature filters
     *
     * @param region	region to check
     */
    private boolean checkFilters(ExtendedProteinRegion region) {
        Feature feat = region.getFeature();
        boolean retVal = this.filters.stream().allMatch(x -> x.filter(feat));
        return retVal;
    }

    @Override
    public int getCellWidth() {
        return this.cellWidth;
    }

    @Override
    public Sort getSort() {
        return this.sortOrder;
    }

    @Override
    public List<String> getGroups(String fid) {
        return this.groupMap.get(fid);
    }



}
