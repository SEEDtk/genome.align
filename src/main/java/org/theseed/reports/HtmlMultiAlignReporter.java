/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.util.List;
import static j2html.TagCreator.*;

import org.theseed.genome.Genome;
import org.theseed.sequence.Sequence;

import j2html.tags.ContainerTag;

/**
 * This is an HTML version of the alignment output.  It produces a set of HTML tables with coloring, and relies on the
 * SEEDtk website style sheets.
 *
 * @author Bruce Parrello
 *
 */
public class HtmlMultiAlignReporter extends MultiAlignReporter {

    // FIELDS
    /** color scheme */
    private AlignColoring scheme;

    public HtmlMultiAlignReporter(File outFile) {
        super(outFile);
        this.scheme = new AlignColoring.Consensus();
    }

    @Override
    public void openReport(Genome genome, String[] altBases) {
        // Write the page header.
        this.println(document().render());
        this.println("<html>");
        this.println(head(title("Alignments"), link().withRel("stylesheet").withHref("https://core.theseed.org/SEEDtk/css/Basic.css")).render());
        // Start the page body.
        this.println("<body>");
    }

    @Override
    public void writeAlignment(String fid, String title, List<Sequence> alignment) {
        ContainerTag alignTable = CoreHtmlUtilities.alignmentTable(alignment, this.scheme);
        this.println(alignTable.render());
    }

    @Override
    public void closeReport() {
        // Close off the page.
        this.println("</body></html>");
    }

}
