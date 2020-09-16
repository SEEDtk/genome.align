/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

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
        TEXT;

        /**
         * @return a reporting object of this type
         *
         * @param output	output stream to receive the report
         */
        public SnipReporter create(OutputStream output) {
            SnipReporter retVal = null;
            switch (this) {
            case TEXT :
                retVal = new TextSnipReporter(output);
                break;
            }
            return retVal;
        }
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
     * Finish the report section for the current alignment.
     */
    protected abstract void closeAlignment();

    /**
     * Process a single set of snips.
     *
     * @param snipCol	descriptor for the snips found
     */
    protected abstract void processSnips(SnipColumn snipCol);

    /**
     * Begin the section of the report for a new alignment.
     *
     * @param title		title of alignment
     * @param regions	regions being aligned
     */
    protected abstract void openAlignment(String title, RegionList regions);


}
