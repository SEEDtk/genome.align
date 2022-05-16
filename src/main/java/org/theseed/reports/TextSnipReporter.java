/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import java.util.List;
import java.util.stream.Collectors;

import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.sequence.RegionList;
import org.theseed.sequence.clustal.SnipColumn;

/**
 * This is a simple text version of the snip reporter.  There are no embellishments, it simply presents each snip column in a
 * single row.
 *
 * @author Bruce Parrello
 */
public class TextSnipReporter extends SnipReporter {

    // FIELDS
    /** title of the current alignment */
    private String title;
    /** buffer for building output lines */
    private StringBuilder buffer;

    public TextSnipReporter(OutputStream output, IParms processor) {
        super(output, processor);
        this.buffer = new StringBuilder(100);
    }

    @Override
    protected void registerGenome(Genome genome) {
    }

    @Override
    protected void openReport(List<GenomeLabel> genomeLabels) {
        // Here we can write the column headers.
        this.println(genomeLabels.stream().map(x -> x.getHeader()).collect(Collectors.joining("\t", "function\t", "")));
    }

    @Override
    protected void closeAlignment() {
    }

    @Override
    protected void processSnips(SnipColumn snipCol) {
        // Here we have a single row.
        buffer.setLength(0);
        buffer.append(this.title);
        for (int i = 0; i < snipCol.getRows(); i++)
            buffer.append('\t').append(snipCol.getSnip(i));
        this.println(buffer.toString());
    }

    @Override
    protected void openAlignment(String title, RegionList regions, Feature feat) {
        this.title = title;
    }

    @Override
    protected void closeReport() {
    }

}
