/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import static j2html.TagCreator.*;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Genome;
import org.theseed.sequence.ExtendedProteinRegion;
import org.theseed.sequence.RegionList;
import org.theseed.sequence.clustal.RealSnipItem;
import org.theseed.sequence.clustal.SnipColumn;

import j2html.tags.DomContent;

/**
 * This class generates a snip report in HTML.  The use of HTML allows significantly more density in information
 * presentation.  While the basic format is equivalent to that of TextSnipReporter, additional information is
 * provided using tooltips and colors.
 *
 * @author Bruce Parrello
 *
 */
public class HtmlSnipReporter extends SnipReporter {

    // FIELDS
    /** current alignment table */
    private HtmlTable<Key.Null> alignment;
    /** list of genome IDs */
    private List<String> genomeIds;
    /** list of HTML segments */
    private List<DomContent> sections;
    /** alignment table column specification array */
    private ColSpec[] cols;
    /** current section title */
    private String title;
    /** page writer */
    private PageWriter writer;
    /** cell width */
    private int cellWidth;
    /** current region list */
    private RegionList regions;
    /** background color for snip differences */
    public static final Color DIFF_COLOR = new Color(0.50, 0.73, 1.0);
    /** style for snip differences */
    public static final String DIFF_STYLE = "background-color: " + DIFF_COLOR.html();
    /** background color for invisible snip differences */
    public static final Color INVISI_COLOR = new Color(0.50, 1.0, 0.73);
    /** style for invisible snip differences */
    public static final String INVISI_STYLE = "background-color: " + INVISI_COLOR.html();
    /** background color for an upstream change */
    public static final Color UPSTREAM_COLOR = new Color(1.0, 0.73, 0.50);
    /** style for upstream snip differences */
    public static final String UPSTREAM_STYLE = "background-color: " + UPSTREAM_COLOR.html();
    /** background color for a gap change */
    public static final Color GAP_COLOR = new Color(1.0, 1.0, 0.0);
    /** style for a gap change */
    public static final String GAP_STYLE = "background-color: " + GAP_COLOR.html();
    /** HTML encoding for a hard break */
    private static final String BREAK_RENDER = br().render();



    public HtmlSnipReporter(OutputStream output, IParms processor) {
        super(output);
        // Create the section list and prime it with the legend.
        this.sections = new ArrayList<DomContent>(100);
        DomContent legend = p(text("Color scheme: "), span("Upstream difference. ").withStyle(UPSTREAM_STYLE),
                span("Gap-related difference. ").withStyle(GAP_STYLE), span("Invisible difference. ").withStyle(INVISI_STYLE),
                span("Protein-modifying difference.").withStyle(DIFF_STYLE));
        this.sections.add(legend);
        // Create the HTML writer.
        this.writer = new FreePageWriter();
        // Get the cell width.
        this.cellWidth = processor.getCellWidth();
    }

    @Override
    protected void registerGenome(Genome genome) {
    }

    @Override
    protected void openReport(List<String> genomeIdList) {
        // Save the genome ID list.  Note that the first genome is the base.
        this.genomeIds = genomeIdList;
        List<ColSpec> colSpecs = this.genomeIds.stream().map(x -> new ColSpec.Aligned(x)).collect(Collectors.toList());
        ColSpec[] cols = new ColSpec[colSpecs.size()];
        this.cols = colSpecs.toArray(cols);
    }

    @Override
    protected void openAlignment(String title, RegionList regions) {
        // Each alignment has its own table.  We start by saving the title.
        this.title = title;
        // Start the table.
        this.alignment = new HtmlTable<Key.Null>(cols);
        // Save the region list.
        this.regions = regions;
    }

    @Override
    protected void processSnips(SnipColumn snipCol) {
        // Count the differences.  If there are none, we skip the row.
        int diffCount = 0;
        // Create a buffer for building table cells.  We need space for each character plus possible
        // colorings.
        StringBuffer buffer = new StringBuffer(snipCol.getWidth() * 5);
        // Get the base-sequence snip.
        String baseSnip = snipCol.getSnip(0).toUpperCase();
        ExtendedProteinRegion baseRegion = this.regions.get(0);
        String[] baseAA = ((RealSnipItem) snipCol.getItem(0)).getProteinMap(baseRegion);
        // Our prototype cells go in here.
        List<DomContent> row = new ArrayList<DomContent>(snipCol.getRows());
        row.add(this.breakUp(baseSnip));
        // Loop through the aligned snips.
        for (int i = 1; i < snipCol.getRows(); i++) {
            // Skip this column if it is not significant.
            if (! snipCol.isSignificant(i))
                row.add(rawHtml("&nbsp;"));
            else {
                // Get the snip text and the protein map.
                String snip = snipCol.getSnip(i).toUpperCase();
                ExtendedProteinRegion region = this.regions.get(snipCol.getFid(i));
                String[] thisAA = ((RealSnipItem) snipCol.getItem(i)).getProteinMap(region);
                // Prepare the output buffer.
                buffer.setLength(0);
                int width = 0;
                // Insert the snip letters, coloring each one that differs from the base.
                for (int p = 0; p < snipCol.getWidth(); p++) {
                    char c = snip.charAt(p);
                    if (width >= this.cellWidth) {
                        buffer.append(BREAK_RENDER);
                        width = 0;
                    }
                    char cBase = baseSnip.charAt(p);
                    if (c == cBase)
                        buffer.append(c);
                    else if (c == '-' || cBase == '-') {
                        buffer.append(markedLetter(c, GAP_STYLE));
                        diffCount++;
                    } else {
                        // Here we have a character change.  We need to determine the amino acid difference.
                        String oldAA = baseAA[p];
                        if (oldAA.contentEquals("upstream")) {
                            buffer.append(markedLetter(c, UPSTREAM_STYLE));
                            diffCount++;
                        } else {
                            String newAA = thisAA[p];
                            // Check for an invisible change.
                            if (! oldAA.equals("upstream") && oldAA.contentEquals(newAA))
                                buffer.append(markedLetter(c, INVISI_STYLE));
                            else {
                                // Here we have a modifying change.  This requires a tooltip AND coloring.
                                DomContent letter = span(rawHtml(markedLetter(c, DIFF_STYLE)
                                        + span(oldAA + " => " + newAA).withClass("tip").render())).withClass("tt");
                                buffer.append(letter.render());
                                diffCount++;
                            }
                        }
                    }
                    width++;
                }
                // Create the cell with the snip in it.
                row.add(rawHtml(buffer.toString()));
            }
        }
        // Only add the row if we found one visible difference.
        if (diffCount > 0) {
            HtmlTable<Key.Null>.Row trow = this.alignment.new Row(Key.NONE);
            row.stream().forEach(x -> trow.add(x));
        }
    }

    /**
     * @return the rendered HTML for a highlighted letter
     *
     * @param c			letter to highlight
     * @param style		highlight style
     */
    protected String markedLetter(char c, String style) {
        return mark(Character.toString(c)).withStyle(style).render();
    }

    /**
     * @return the incoming string broken into pieces based on the cell width
     *
     * @param snipText	string to break up.
     */
    private DomContent breakUp(String snipText) {
        int snipLen = snipText.length();
        StringBuffer breakBuffer = new StringBuffer(snipLen + BREAK_RENDER.length() * (snipLen / this.cellWidth + 1));
        breakBuffer.append(StringUtils.substring(snipText, 0, this.cellWidth));
        int consumed = breakBuffer.length();
        while (consumed < snipLen) {
            int end = consumed + this.cellWidth;
            breakBuffer.append(BREAK_RENDER).append(StringUtils.substring(snipText, consumed, end));
            consumed = end;
        }
        return rawHtml(breakBuffer.toString());
    }

    @Override
    protected void closeAlignment() {
        // Output the table.
        if (this.alignment.getHeight() > 0) {
            this.sections.add(h2(title));
            this.sections.add(this.alignment.output());
        }
    }

    @Override
    protected void closeReport() {
        // Form the tables into a page.
        DomContent block = this.writer.highlightBlock(this.sections);
        this.writer.writePage("Snip Alignments", h1("Snip Alignments"), block);
    }

}
