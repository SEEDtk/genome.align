/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import static j2html.TagCreator.*;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
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
    /** list of genome IDs */
    private List<String> genomeIds;
    /** current table section */
    private TableEntry table;
    /** list of table sections */
    private List<TableEntry> sections;
    /** alignment table column specification array */
    private ColSpec[] cols;
    /** page writer */
    private PageWriter writer;
    /** cell width */
    private int cellWidth;
    /** tracker for number of genome differences */
    private BitSet diffs;
    /** current region list */
    private RegionList regions;
    /** map of genome IDs to names */
    private Map<String, String> gNameMap;
    /** sorter to use for output */
    private Comparator<TableEntry> sorter;
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

    /**
     * Utility class for sorting the output by number of differences.
     */
    private class TableEntry {

        /** difference count */
        private int diffCount;
        /** base feature location */
        private Location loc;
        /** section title */
        private String title;
        /** section table */
        private HtmlTable<Key.Null> alignment;

        /**
         * Create a new table entry.
         *
         * @param title		title of this table
         * @param regions	source region list
         */
        protected TableEntry(String title, RegionList regions) {
            this.diffCount = 0;
            Feature feat = regions.get(0).getFeature();
            this.loc = feat.getLocation();
            this.title = title;
            this.alignment = new HtmlTable<Key.Null>(HtmlSnipReporter.this.cols);
        }

        /**
         * Add a row to the alignment table.
         *
         * @param locString		location string for the snip in this row
         * @param row			row of snip HTML strings to add
         */
        public void addRow(String locString, List<DomContent> row) {
            HtmlTable<Key.Null>.Row trow = this.alignment.new Row(Key.NONE).add(locString);
            row.stream().forEach(x -> trow.add(x));
        }

        /**
         * @return the number of rows in the alignment table.
         */
        public int getHeight() {
            return this.alignment.getHeight();
        }

        /**
         * @return the output for this table entry
         */
        public DomContent output() {
            return div(h2(this.title), this.alignment.output());
        }

        /**
         * Update the difference count.
         *
         * @param diff	number of modified genomes in this alignment
         */
        public void setDiffCount(int diff) {
            this.diffCount = diff;
        }

    }

    /**
     * Sort types for table entries.
     */
    public static enum Sort {
        LOCATION, CHANGES;

        /**
         * @return a sorter of the specified type
         */
        protected Comparator<TableEntry> create() {
            Comparator<TableEntry> retVal = null;
            switch (this) {
            case LOCATION :
                retVal = new LocationSort();
                break;
            case CHANGES :
                retVal = new DiffSort();
                break;
            }
            return retVal;
        }

    }

    /**
     * Table entry sorter for location-based sorts
     */
    public static class LocationSort implements Comparator<TableEntry> {

        @Override
        public int compare(TableEntry o1, TableEntry o2) {
            int retVal = o1.loc.compareTo(o2.loc);
            if (retVal == 0)
                retVal = o1.title.compareTo(o2.title);
            return retVal;
        }

    }

    /**
     * Table entry sorter for difference-based sorts
     */
    public static class DiffSort implements Comparator<TableEntry> {

        @Override
        public int compare(TableEntry o1, TableEntry o2) {
            int retVal = o2.diffCount - o1.diffCount;
            if (retVal == 0) {
                retVal = o1.loc.compareTo(o2.loc);
                if (retVal == 0)
                    retVal = o1.title.compareTo(o2.title);
            }
            return retVal;
        }

    }


    public HtmlSnipReporter(OutputStream output, IParms processor) {
        super(output);
        // Create the HTML writer.
        this.writer = new FreePageWriter();
        // Get the cell width.
        this.cellWidth = processor.getCellWidth();
        // Create the genome map.
        this.gNameMap = new HashMap<String, String>();
        // Create the table list.
        this.sections = new ArrayList<TableEntry>(1000);
        // Create the output sorter.
        this.sorter = processor.getSort().create();
    }

    @Override
    protected void registerGenome(Genome genome) {
        this.gNameMap.put(genome.getId(), genome.getName());
    }

    @Override
    protected void openReport(List<String> genomeIdList) {
        // Save the genome ID list.  Note that the first genome is the base, and the location column precedes it.
        this.genomeIds = genomeIdList;
        List<ColSpec> colSpecs = new ArrayList<ColSpec>(this.genomeIds.size() + 1);
        colSpecs.add(new ColSpec.Normal("Location"));
        this.genomeIds.stream().forEach(x -> colSpecs.add(new ColSpec.Aligned(x).setTip(this.gNameMap.get(x))));
        ColSpec[] cols = new ColSpec[colSpecs.size()];
        this.cols = colSpecs.toArray(cols);
        // Create the difference bitmap.
        this.diffs = new BitSet(this.genomeIds.size());
    }

    @Override
    protected void openAlignment(String title, RegionList regions) {
        // Each alignment has its own table.
        this.table = this.new TableEntry(title, regions);
        // Save the region list.
        this.regions = regions;
        // Clear the difference bitmap.
        this.diffs.clear();
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
                        buffer.append(markedLetter(c, GAP_STYLE).render());
                        diffCount++;
                        this.diffs.set(i);
                    } else {
                        // Here we have a character change.  We need to determine the amino acid difference.
                        String oldAA = baseAA[p];
                        if (oldAA.contentEquals("upstream")) {
                            buffer.append(markedLetter(c, UPSTREAM_STYLE).render());
                            diffCount++;
                            this.diffs.set(i);
                        } else {
                            String newAA = thisAA[p];
                            // Check for an invisible change.
                            if (! oldAA.equals("upstream") && oldAA.contentEquals(newAA))
                                buffer.append(markedLetter(c, INVISI_STYLE).render());
                            else {
                                // Here we have a modifying change.  This requires a tooltip AND coloring.
                                DomContent letter = CoreHtmlUtilities.toolTip(markedLetter(c, DIFF_STYLE),
                                        oldAA + " => " + newAA);
                                buffer.append(letter.render());
                                diffCount++;
                                this.diffs.set(i);
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
            this.table.addRow(snipCol.getLocString(0), row);
        }
    }

    /**
     * @return the rendered HTML for a highlighted letter
     *
     * @param c			letter to highlight
     * @param style		highlight style
     */
    protected DomContent markedLetter(char c, String style) {
        return mark(Character.toString(c)).withStyle(style);
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
        // Output the table.  Note we record the number of changed genomes, which is used for sorting.
        if (this.table.getHeight() > 0) {
            this.table.setDiffCount(this.diffs.cardinality());
            this.sections.add(this.table);
        }
    }

    @Override
    protected void closeReport() {
        // Form the tables into a page.
        DomContent legend = p(text("Color scheme: "), span("Upstream difference. ").withStyle(UPSTREAM_STYLE),
                span("Gap-related difference. ").withStyle(GAP_STYLE), span("Invisible difference. ").withStyle(INVISI_STYLE),
                span("Protein-modifying difference.").withStyle(DIFF_STYLE));
        this.sections.sort(this.sorter);
        DomContent tables = div().with(this.sections.stream().map(x -> x.output()));
        DomContent block = this.writer.highlightBlock(legend, tables);
        this.writer.writePage("Snip Alignments", h1("Snip Alignments"), block);
    }

}
