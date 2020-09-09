/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import java.util.List;

import org.theseed.genome.Genome;
import org.theseed.sequence.Sequence;

/**
 * This is the base class for a multiple-alignment report.
 *
 * @author Bruce Parrello
 *
 */
public abstract class MultiAlignReporter extends BaseReporter {

    /**
     * Construct a report object for a specified output stream.
     *
     * @param outStream		output stream to receive the report
     */
    public MultiAlignReporter(OutputStream outStream) {
        super(outStream);
    }

    /**
     * Initialize the report.
     *
     * @param genome	base genome
     * @param altBases 	array of IDs for additional base genomes
     */
    public abstract void openReport(Genome genome, String[] altBases);

    /**
     * Output a particular alignment.
     *
     * @param title			name of the alignment
     * @param alignment	aligned sequences
     */
    public abstract void writeAlignment(String title, List<Sequence> alignment);

    /**
     * Finish the report.
     */
    public abstract void closeReport();

    /**
     * Enum for report types
     */
    public static enum Type {
            TEXT, HTML, SNIPS;

        public MultiAlignReporter create(OutputStream outStream) {
            MultiAlignReporter retVal = null;
            switch (this) {
            case TEXT :
                retVal = new TextMultiAlignReporter(outStream);
                break;
            case HTML :
                retVal = new HtmlMultiAlignReporter(outStream);
                break;
            case SNIPS :
                retVal = new SnipMultiAlignReporter(outStream);
            }
            return retVal;
        }
    }

}
