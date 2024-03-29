/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.sequence.ExtendedProteinRegion;
import org.theseed.sequence.RegionList;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.clustal.ISnipItem;
import org.theseed.sequence.clustal.SnipColumn;
import org.theseed.sequence.clustal.SnipIterator;

/**
 * This is the base class for snip reports.  The snip report presumes a multiple-sequence alignment involving a base genome and a fixed set of aligned genomes.
 * Numerous alignments will be processed, each based on a single base genome feature and its upstream region.
 *
 * This method is also responsible for writing the feature data output file, which indicates which genomes have significant snips
 * in each feature.  For each feature, the output file has one indicator for upstream and one for instream.  An "M" (mutation)
 * indicates only nucleotide changes.  A "D" (deletion) indicates the presence of gap characters.
 *
 * @author Bruce Parrello
 *
 */
public abstract class SnipReporter extends BaseReporter {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SnipReporter.class);
    /** IDs of the aligned genomes */
    private List<GenomeLabel> genomeLabels;
    /** feature data output file */
    private PrintWriter fDataOut;
    /** controlling processor */
    private IParms processor;
    /** map of genome IDs to names */
    private Map<String, String> gNameMap;
    /** feature data output flag string for unmodified features */
    private String fDataUnmodified;


    /**
     * Enum for different report formats
     */
    public static enum Type {
        TEXT, HTML, MAJORPROTEIN, MAJORUPSTREAM;

        /**
         * @return a reporting object of this type
         *
         * @param output		output stream to receive the report
         * @param processor		constructing client
         */
        public SnipReporter create(OutputStream output, IParms processor) {
            SnipReporter retVal = null;
            switch (this) {
            case TEXT :
                retVal = new TextSnipReporter(output, processor);
                break;
            case HTML :
                retVal = new HtmlSnipReporter(output, processor);
                break;
            case MAJORPROTEIN :
                retVal = new MajorSnipReporter.Protein(output, processor);
                break;
            case MAJORUPSTREAM :
                retVal = new MajorSnipReporter.Upstream(output, processor);
                break;
            }
            return retVal;
        }
    }

    /**
     * Interface used by constructing clients to pass along parameters.
     */
    public interface IParms {

        /**
         * @return the maximum width for text in an alignment cell
         */
        int getCellWidth();

        /**
         * @return the sort order for the HTML tables
         */
        HtmlSnipReporter.Sort getSort();

        /**
         * @return the group list for the specified feature, or NULL if it is not in a group
         *
         * @param fid	feature of interest
         */
        List<String> getGroups(String fid);

        /**
         * @return the set of IDs for special genomes in the MAJOR report
         */
        Set<String> getSpecial();

    }

    /**
     * Create a snip report on a specified output stream.
     *
     * @param output		output stream to receive the report
     * @param processor		controlling command processor
     */
    public SnipReporter(OutputStream output, IParms processor) {
        super(output);
        this.genomeLabels = new ArrayList<GenomeLabel>();
        this.fDataOut = null;
        this.processor = processor;
        // Create the genome map and ID set.
        this.gNameMap = new HashMap<String, String>();
    }

    /**
     * Register the feature data output file.
     *
     * @param	outFile		feature data output file
     *
     * @throws FileNotFoundException
     */
    public void setupFeatureOutput(File outFile) throws FileNotFoundException {
        this.fDataOut = new PrintWriter(outFile);
    }
    /**
     * Register an aligned genome.
     *
     * @param genome	genome to register
     */
    public void register(Genome genome) {
        this.genomeLabels.add(new GenomeLabel(genome));
        this.gNameMap.put(genome.getId(), genome.getName());
        this.registerGenome(genome);
    }

    /**
     * Extract useful data from an aligned genome.
     *
     * @param genome	genome to register
     */
    protected abstract void registerGenome(Genome genome);

    /**
     * Reorder the genome IDs according to the input list.  The first genome ID is not changed, since that's the base.
     *
     * @param orderingList	list of genome IDs representing the desired order
     */
    public void reorder(List<GenomeLabel> orderingList) {
        // Create a set of all the existing genomes.
        Set<GenomeLabel> oldGenomes = this.genomeLabels.stream().collect(Collectors.toSet());
        // Copy the base genome to the output list.
        GenomeLabel baseGenome = this.genomeLabels.get(0);
        this.genomeLabels.clear();
        this.genomeLabels.add(baseGenome);
        oldGenomes.remove(baseGenome);
        // Now loop through the ordering list, adding all the genomes that are legitimate.
        for (GenomeLabel genome : orderingList) {
            if (oldGenomes.contains(genome)) {
                this.genomeLabels.add(genome);
                oldGenomes.remove(genome);
            }
        }
        // Now add the residual.
        this.genomeLabels.addAll(oldGenomes);
        // Write the names to the feature-data output.
        if (this.fDataOut != null) {
            for (GenomeLabel genomeData : this.genomeLabels) {
                String genomeId = genomeData.getId();
                this.fDataOut.println(genomeId + "\t" + this.getGName(genomeId));
            }
            this.fDataOut.println("//");
        }
    }

    /**
     * Begin the alignment portion of the report.
     */
    public void initializeOutput() {
        this.openReport(this.genomeLabels);
        // Create a record of empty fields for the feature-data output file.
        this.fDataUnmodified = StringUtils.repeat("  \t", this.genomeLabels.size() - 1) + "  ";
    }

    /**
     * Start the alignment portion of the report.
     *
     * @param genomeIdList	list of aligned genomes, in order
     */
    protected abstract void openReport(List<GenomeLabel> genomeIdList);

    /**
     * Process an alignment.
     *
     * @param feat			base genome feature
     * @param regions		list of regions that were aligned
     * @param alignment		list of aligned sequences
     */
    public void processAlignment(Feature feat, RegionList regions, List<Sequence> alignment) {
        // Start this section of the report.
        String function = feat.getPegFunction();
        this.openAlignment(function, regions, feat);
        // We need to compute the wild strain genomes.  The first region is for a wild genome, and all sequences
        // in the alignment that are not in the main genome ID list are wild as well.  Note the the wild set
        // may change between alignments if a wild strain is missing a region in this alignment run.
        Set<String> wildSet = new TreeSet<String>();
        // Get the genome IDs in order.
        List<String> genomeIds = this.genomeLabels.stream().map(x -> x.getId()).collect(Collectors.toList());
        // Start with the base genome.
        wildSet.add(this.genomeLabels.get(0).getId());
        // Add the genomes not being displayed.
        regions.stream().map(x -> Feature.genomeOf(x.getLabel())).filter(g -> ! genomeIds.contains(g)).forEach(g -> wildSet.add(g));
        // Now we need to iterate through the snips.
        SnipIterator.Run snipRun = new SnipIterator.Run(regions, alignment, wildSet, genomeIds);
        // This will track the genomes that have significant snips.
        String[] modifiedGenomes = new String[this.genomeLabels.size()];
        Arrays.fill(modifiedGenomes, "  ");
        int count = 0;
        for (SnipColumn snipCol : snipRun) {
            this.processSnips(snipCol);
            String baseChars = snipCol.getSnip(0);
            for (int i = 1; i < snipCol.getRows(); i++)
                modifiedGenomes[i] = this.charCode(snipCol.getItem(i), regions.get(snipCol.getFid(i)), baseChars, modifiedGenomes[i]);
            count++;
        }
        String fid = feat.getId();
        log.debug("{} snips found in alignment for {}.", count, fid);
        // Form the flags for the featureData file.
        String flags = IntStream.range(0, this.genomeLabels.size()).mapToObj(i -> modifiedGenomes[i])
                .collect(Collectors.joining("\t"));
        // Write out the featureData record.
        writeFeatureData(fid, flags);
        // Finish off the alignment.
        this.closeAlignment();
    }

    /**
     * Write the featureData record for this feature.
     *
     * @param fid		ID of the feature to output
     * @param flags		flag string of changes associated with the feature
     */
    private void writeFeatureData(String fid, String flags) {
        if (this.fDataOut != null) {
            List<String> groupList = new ArrayList<String>();
            List<String> mainList = this.processor.getGroups(fid);
            if (mainList != null)
                groupList.addAll(mainList);
            List<String> others = this.getOtherGroups(fid);
            if (others != null)
                groupList.addAll(others);
            this.fDataOut.format("%s\t%s\t%s%n", fid, StringUtils.join(groupList, ","), flags);
        }
    }

    /**
     * Write the featureData record for an unmodified feature.
     *
     * @param fid		ID of the feature to output
     */
    public void writeFeatureData(String fid) {
        this.writeFeatureData(fid, this.fDataUnmodified);
    }

    /**
     * @return the character display code for this type of snip
     *
     * @param snip			snip item to check
     * @param region		region containing the snip (or NULL if virtual)
     * @param baseChars		base genome sequence for this snip
     * @param original		original display code for this snip
     */
    private String charCode(ISnipItem snip, ExtendedProteinRegion region, String baseChars, String original) {
        String retVal = original;
        // A character code of "D" overrides everything else.  Also, if this change is not
        // significant, we skip it.
        if (snip.isSignificant()) {
            // Here we want to look for "D" (gap), "M" (change).
            String snipString = snip.getChars();
            int snipLoc = snip.getOffset();
            int snipLen = snip.getLen();
            if (snip.isReal(baseChars, region)) {
                // Determine if we are upstream (0) or instream (1). The result is an index into
                // the character relevant to this change.
                int type = (snipLoc + snipLen > region.getUpstreamDistance() ? 1 : 0);
                // If we are already a "D", we stay that way. Otherwise we test.
                char[] myChars = original.toCharArray();
                if (myChars[type] != 'D') {
                    if (snipString.contains("-"))
                        myChars[type] = 'D';
                    else
                        myChars[type] = 'M';
                    retVal = String.valueOf(myChars);
                }
            }
        }
        return retVal;
    }

    /**
     * Begin the section of the report for a new alignment.
     *
     * @param title		title of alignment
     * @param regions	regions being aligned
     * @param feat		base genome feature
     */
    protected abstract void openAlignment(String title, RegionList regions, Feature feat);

    /**
     * Process a single set of snips.
     *
     * @param snipCol	descriptor for the snips found
     */
    protected abstract void processSnips(SnipColumn snipCol);

    /**
     * Finish the report section for the current alignment.
     */
    protected abstract void closeAlignment();

    /**
     * Finish the entire report.
     */
    protected abstract void closeReport();

    public void finishReport() {
        this.closeReport();
        this.getWriter().flush();
        if (this.fDataOut != null)
            this.fDataOut.close();
    }

    /**
     * @return the processor
     */
    public IParms getProcessor() {
        return processor;
    }

    /**
     * @return the name of the specified genome
     *
     * @param id	ID of genome whose name is desired
     */
    public String getGName(String id) {
        return this.gNameMap.get(id);
    }

    /**
     * @return a list of additional groups, or NULL if there are none
     *
     * This must be overridden by the subclass to return anything.
     */
    protected List<String> getOtherGroups(String fid) {
        return null;
    }

}
