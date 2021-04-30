/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import static j2html.TagCreator.*;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.magic.MagicMap;
import org.theseed.magic.MagicObject;
import org.theseed.sequence.ExtendedProteinRegion;
import org.theseed.sequence.RegionList;
import org.theseed.sequence.clustal.RealSnipItem;
import org.theseed.sequence.clustal.SnipColumn;
import org.theseed.web.ColSpec;
import org.theseed.web.HtmlTable;
import org.theseed.web.Key;
import org.theseed.web.Row;

import j2html.tags.ContainerTag;
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
    /** sorter to use for output */
    private Comparator<TableEntry> sorter;
    /** subsystem ID map */
    private MagicMap<Subsystem> subMap;
    /** map of feature IDs to subsystem IDs */
    private Map<String, List<String>> fidSubMap;
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
    /** background color for a gap off the edge of the contig */
    private static final Color EDGE_COLOR = new Color(1.0, 0.73, 0.73);
    /** style for a gap off the edge of the contig */
    public static final String EDGE_STYLE = "background-color: " + EDGE_COLOR.html();
    /** HTML encoding for a hard break */
    private static final String BREAK_RENDER = br().render();
    /** location of the group page */
    private static final String GROUP_URL = "http://core.theseed.org/SEEDtk/rna.cgi/groups?group=";

    /**
     * Utility class for tracking subsystem links.
     */
    private static class Subsystem extends MagicObject {

        /** serialization ID */
        private static final long serialVersionUID = 6147612938895322808L;

        /**
         * Construct a subsystem object from a subsystem name.
         *
         * @param name	subsystem name
         */
        public Subsystem(String name) {
            super(null, name);
        }

        @Override
        protected String normalize() {
            return this.getName();
        }

        /**
         * @return a link to the subsystem group page
         */
        public DomContent getLink() {
            String title = LinkObject.Core.cleanSubsystemName(this.getName());
            DomContent retVal = a(this.getName()).withHref(GROUP_URL + this.getId() + ";title=" + title).withTarget("_blank");
            return retVal;
        }

    }

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
        /** title link */
        private DomContent titleLink;
        /** subsystem list */
        private DomContent subsystems;
        /** group list */
        private List<String> groups;
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
            this.titleLink = feat.getParent().getLinker().featureLink(feat.getId(), text(title));
            Collection<String> subsysList = feat.getSubsystems();
            if (subsysList.isEmpty()) {
                this.subsystems = null;
            } else {
                // This is very complicated.  We need to create a bunch of subsystem links and we need
                // the subsystem IDs for both the links and for eventual output to the groups file.
                List<String> subIds = new ArrayList<String>(subsysList.size());
                List<Subsystem> subsToLink = new ArrayList<Subsystem>(subsysList.size());
                for (String subsysName : subsysList) {
                    Subsystem sub = HtmlSnipReporter.this.subMap.getByName(subsysName);
                    if (sub == null) {
                        sub = new Subsystem(subsysName);
                        HtmlSnipReporter.this.subMap.put(sub);
                    }
                    subIds.add(sub.getId());
                    subsToLink.add(sub);
                }
                // Create the hyperlinked display string for the subsystems.
                this.subsystems = HtmlSnipReporter.this.subsystemList(subsToLink);
                // Save the subsystem ID list so the base class can put it in the group file.
                HtmlSnipReporter.this.fidSubMap.put(feat.getId(), subIds);
            }
            this.alignment = new HtmlTable<Key.Null>(HtmlSnipReporter.this.cols);
            this.groups = HtmlSnipReporter.this.getProcessor().getGroups(feat.getId());
        }

        /**
         * Add a row to the alignment table.
         *
         * @param locString		location string for the snip in this row
         * @param row			row of snip HTML strings to add
         */
        public void addRow(DomContent locString, List<DomContent> row) {
            Row<Key.Null> trow = new Row<Key.Null>(this.alignment, Key.NONE).add(locString);
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
            ContainerTag modList = ul();
            if (this.groups != null)
                modList.with(li(createGroupLinks(this.groups)));
            if (this.subsystems != null)
                modList.with(li(this.subsystems));
            return div(h2(this.titleLink), modList, this.alignment.output());
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
        super(output, processor);
        // Create the HTML writer.
        this.writer = new FreePageWriter(this.getWriter());
        // Get the cell width.
        this.cellWidth = processor.getCellWidth();
        // Create the table list.
        this.sections = new ArrayList<TableEntry>(1000);
        // Create the output sorter.
        this.sorter = processor.getSort().create();
        // Create the subsystem maps.
        this.subMap = new MagicMap<Subsystem>(new Subsystem(""));
        this.fidSubMap = new HashMap<String, List<String>>(1000);
    }

    /**
     * @return a list of subsystem links for the specified subsystems
     *
     * @param subsToLink	list of subsystems
     */
    public DomContent subsystemList(List<Subsystem> subsToLink) {
        return HtmlUtilities.joinDelimited(subsToLink.stream().map(x -> x.getLink()).collect(Collectors.toList()), " | ");
    }

    /**
     * Create a list of groups linked to the appropriate group display pages
     *
     * @param groups	list of groups for the current feature
     *
     * @return a comma-delimited list of hyperlinks for the groups containing the current feature
     */
    public DomContent createGroupLinks(List<String> groups) {
        DomContent retVal = HtmlUtilities.joinDelimited(groups.stream()
                .map(x -> a(x).withHref(GROUP_URL + x).withTarget("_blank")), ", ");
        return retVal;
    }

    @Override
    protected void registerGenome(Genome genome) {
    }

    @Override
    protected void openReport(List<String> genomeIdList) {
        // Save the genome ID list.  Note that the first genome is the base, and the location column precedes it.
        this.genomeIds = genomeIdList;
        List<ColSpec> colSpecs = new ArrayList<ColSpec>(this.genomeIds.size() + 1);
        colSpecs.add(new ColSpec.Normal("Location"));
        this.genomeIds.stream().forEach(x -> colSpecs.add(new ColSpec.Aligned(x).setTip(this.getGName(x))));
        ColSpec[] cols = new ColSpec[colSpecs.size()];
        this.cols = colSpecs.toArray(cols);
        // Create the difference bitmap.
        this.diffs = new BitSet(this.genomeIds.size());
    }

    @Override
    protected void openAlignment(String title, RegionList regions, Feature feat) {
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
        Feature baseFeat = baseRegion.getFeature();
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
                // This tracks our real region offset.
                int regionOffset = snipCol.getOffset(i);
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
                        String style = GAP_STYLE;
                        // Check for a virtual location, indicating the difference is edge-related.
                        if (c == '-' && region.isVirtual(regionOffset))
                            style = EDGE_STYLE;
                        buffer.append(markedLetter(c, style).render());
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
                    if (c != '-') regionOffset++;
                    width++;
                }
                // Create the cell with the snip in it.
                row.add(rawHtml(buffer.toString()));
            }
        }
        // Only add the row if we found a visible difference.
        if (diffCount > 0) {
            // Determine if this is an upstream or instream snip.
            DomContent label;
            String labelString = snipCol.getLocString(0);
            Location snipLoc = snipCol.getLoc(0);
            if (snipLoc.isOverlapping(baseFeat.getLocation()))
                label = b(labelString);
            else
                label = i(labelString);
            this.table.addRow(label, row);
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
                span("Protein-modifying difference. ").withStyle(DIFF_STYLE),
                span("Contig edge difference. ").withStyle(EDGE_STYLE));
        this.sections.sort(this.sorter);
        DomContent tables = div().with(this.sections.stream().map(x -> x.output()));
        DomContent block = this.writer.highlightBlock(legend, tables);
        this.writer.writePage("Snip Alignments", h1("Snip Alignments"), block);
    }

    @Override
    protected List<String> getOtherGroups(String fid) {
        return this.fidSubMap.get(fid);
    }
}
