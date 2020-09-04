/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import java.util.List;

import org.theseed.genome.Genome;
import org.theseed.sequence.Sequence;

/**
 * This reporting object produces a simple text version of the alignments.
 *
 * @author Bruce Parrello
 *
 */
public class TextMultiAlignReporter extends MultiAlignReporter {

    public TextMultiAlignReporter(OutputStream outStream) {
        super(outStream);
    }

    @Override
    public void openReport(Genome genome) {
    }

    @Override
    public void writeAlignment(String title, List<Sequence> alignment) {
        this.println(title);
        this.println();
        for (Sequence seq : alignment)
            this.print("%s\t%s\t%s", seq.getLabel(), seq.getComment(), seq.getSequence());
        this.println();
    }

    @Override
    public void closeReport() {
    }

}
