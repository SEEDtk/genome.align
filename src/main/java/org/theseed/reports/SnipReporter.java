/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

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
 * @author Bruce Parrello
 *
 */
public abstract class SnipReporter extends BaseReporter {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SnipReporter.class);
    /** IDs of the aligned genomes */
    private List<String> genomeIds;

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

    }

    /**
     * Create a snip report on a specified output stream.
     *
     * @param output	output stream to receive the report
     */
    public SnipReporter(OutputStream output) {
        super(output);
        this.genomeIds = new ArrayList<String>();
    }

    /**
     * Register an aligned genome.
     *
     * @param genome	genome to register
     */
    public void register(Genome genome) {
        this.genomeIds.add(genome.getId());
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
     * @param regions		list of regions that were aligned
     * @param alignment		list of aligned sequences
     */
    public void processAlignment(String title, RegionList regions, List<Sequence> alignment) {
        // Start this section of the report.
        this.openAlignment(title, regions);
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
        int count = 0;
        for (SnipColumn snipCol : snipRun) {
            this.processSnips(snipCol);
            count++;
        }
        log.debug("{} snips found in alignment for {}.", count, title);
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

    /**
     *
     */
    public void finishReport() {
        this.closeReport();
    }

}
