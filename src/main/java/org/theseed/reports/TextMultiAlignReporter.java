/**
 *
 */
package org.theseed.reports;

import java.io.File;
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

    public TextMultiAlignReporter(File outFile) {
        super(outFile);
    }

    @Override
    public void openReport(Genome genome, String[] altBases) {
    }

    @Override
    public void writeAlignment(String fid, String title, List<Sequence> alignment) {
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
