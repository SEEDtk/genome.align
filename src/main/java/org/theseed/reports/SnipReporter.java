/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.BitSet;
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
import org.theseed.sequence.RegionList;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.clustal.SnipColumn;
import org.theseed.sequence.clustal.SnipIterator;

/**
 * This is the base class for snip reports.  The snip report presumes a multiple-sequence alignment involving a base genome and a fixed set of aligned genomes.
 * Numerous alignments will be processed, each based on a single base genome feature and its upstream region.
 *
 * This method is also responsible for writing the feature data output file, which indicates which genomes have significant snips
 * in each feature.
 *
 * @author Bruce Parrello
 *
 */
public abstract class SnipReporter extends BaseReporter {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SnipReporter.class);
    /** IDs of the aligned genomes */
    private List<String> genomeIds;
    /** feature data output file */
    private PrintWriter fDataOut;
    /** controlling processor */
    private IParms processor;
    /** map of genome IDs to names */
    private Map<String, String> gNameMap;


    /**
     * Enum for different report formats
     */
    public static enum Type {
        TEXT, HTML;

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

    }

    /**
     * Create a snip report on a specified output stream.
     *
     * @param output		output stream to receive the report
     * @param processor		controlling command processor
     */
    public SnipReporter(OutputStream output, IParms processor) {
        super(output);
        this.genomeIds = new ArrayList<String>();
        this.fDataOut = null;
        this.processor = processor;
        // Create the genome map.
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
        this.genomeIds.add(genome.getId());
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
    public void reorder(List<String> orderingList) {
        // Create a set of all the existing genomes.
        Set<String> oldGenomes = this.genomeIds.stream().collect(Collectors.toSet());
        // Copy the base genome to the output list.
        String baseGenome = this.genomeIds.get(0);
        this.genomeIds.clear();
        this.genomeIds.add(baseGenome);
        oldGenomes.remove(baseGenome);
        // Now loop through the ordering list, adding all the genomes that are legitimate.
        for (String genome : orderingList) {
            if (oldGenomes.contains(genome)) {
                this.genomeIds.add(genome);
                oldGenomes.remove(genome);
            }
        }
        // Now add the residual.
        this.genomeIds.addAll(oldGenomes);
        // Write the names to the feature-data output.
        if (this.fDataOut != null) {
            for (String genomeId : this.genomeIds)
                this.fDataOut.println(genomeId + "\t" + this.getGName(genomeId));
            this.fDataOut.println("//");
        }
    }

    /**
     * Begin the alignment portion of the report.
     */
    public void initializeOutput() {
        this.openReport(this.genomeIds);
    }

    /**
     * Start the alignment portion of the report.
     *
     * @param genomeIdList	list of aligned genomes, in order
     */
    protected abstract void openReport(List<String> genomeIdList);

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
        this.openAlignment(function, regions);
        // We need to compute the wild strain genomes.  The first region is a wild genome, and all sequences
        // in the alignment that are not in the main genome ID list are wild as well.  Note the the wild set
        // may change between alignments if a wild strain is missing a region in this alignment run.
        Set<String> wildSet = new TreeSet<String>();
        // Start with the base genome.
        wildSet.add(this.genomeIds.get(0));
        // Add the genomes not being displayed.
        regions.stream().map(x -> Feature.genomeOf(x.getLabel())).filter(g -> ! this.genomeIds.contains(g)).forEach(g -> wildSet.add(g));
        // Now we need to iterate through the snips.
        SnipIterator.Run snipRun = new SnipIterator.Run(regions, alignment, wildSet, this.genomeIds);
        // This will track the genomes that have significant snips.
        BitSet modifiedGenomes = new BitSet(this.genomeIds.size());
        int count = 0;
        for (SnipColumn snipCol : snipRun) {
            this.processSnips(snipCol);
            String baseChars = snipCol.getSnip(0);
            for (int i = 1; i < snipCol.getRows(); i++)
                modifiedGenomes.set(i, snipCol.isSignificant(i) && ! snipCol.getSnip(i).contentEquals(baseChars));
            count++;
        }
        String fid = feat.getId();
        log.debug("{} snips found in alignment for {}.", count, fid);
        if (this.fDataOut != null) {
            // Here we need to update the featureData file.
            String flags = IntStream.range(0, this.genomeIds.size()).mapToObj(i -> (modifiedGenomes.get(i) ? "X" : ""))
                    .collect(Collectors.joining("\t"));
            List<String> groupList = new ArrayList<String>();
            List<String> mainList = this.processor.getGroups(fid);
            if (mainList != null)
                groupList.addAll(mainList);
            List<String> others = this.getOtherGroups(fid);
            if (others != null)
                groupList.addAll(others);
            this.fDataOut.format("%s\t%s\t%s%n", fid, StringUtils.join(groupList, ","), flags);
        }
        // Finish off the alignment.
        this.closeAlignment();
    }

    /**
     * Begin the section of the report for a new alignment.
     *
     * @param title		title of alignment
     * @param regions	regions being aligned
     */
    protected abstract void openAlignment(String title, RegionList regions);

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
