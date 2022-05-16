/**
 *
 */
package org.theseed.reports;

import java.io.IOException;
import java.io.OutputStream;
import java.io.UncheckedIOException;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.sequence.ExtendedProteinRegion;
import org.theseed.sequence.RegionList;
import org.theseed.sequence.clustal.SnipColumn;
import org.apache.commons.lang3.StringUtils;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.DataFormat;
import org.apache.poi.ss.usermodel.FillPatternType;
import org.apache.poi.ss.usermodel.HorizontalAlignment;
import org.apache.poi.ss.usermodel.IndexedColors;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This report lists genes that have modifying changes across all genomes in the special subset.
 * The output is in Excel format.  The columns of the sheet are as follows:
 *
 * 		fig ID			feature ID
 * 		start_loc		start location
 * 		stop_loc		stop location
 * 		strand			strand containing gene
 * 		gene_name		common gene name
 * 		length			DNA length of feature
 * 		groups			comma-delimited list of groups (modulons, regulons, subsystems)
 *		all				X if changed in all strains
 *
 * @author Bruce Parrello
 *
 */
public abstract class MajorSnipReporter extends SnipReporter {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(MajorSnipReporter.class);
    /** current workbook */
    private Workbook workbook;
    /** current worksheet */
    private Sheet worksheet;
    /** style for integer cells */
    private CellStyle numStyle;
    /** style for header cells */
    private CellStyle headStyle;
    /** style for flag cells */
    private CellStyle flagStyle;
    /** style for group cells */
    private CellStyle groupStyle;
    /** next row number */
    private int rowNum;
    /** current row */
    private org.apache.poi.ss.usermodel.Row ssRow;
    /** set of IDs for special genomes */
    private Set<String> specials;
    /** set of IDs for all genomes */
    private Set<String> all;
    /** set of genomes with protein changes in the current alignment */
    private Set<String> changed;
    /** current feature */
    private Feature feat;
    /** ID of the base genome */
    private String baseId;
    /** controlling processor */
    private SnipReporter.IParms processor;
    /** saved output stream */
    private OutputStream outStream;
    /** gene name match pattern */
    private static final Pattern GENE_NAME = Pattern.compile("[a-z]{3}(?:[A-Z])?");
    /** function column width */
    private static int TEXT_WIDTH = 30 * 256;

    /**
     * Set up the Excel output and get the special subset IDs from the command processor.
     *
     * @param output		output stream
     * @param processor		controlling processor
     */
    public MajorSnipReporter(OutputStream output, IParms processor) {
        super(output, processor);
        // Denote we do not yet know the base genome.
        this.baseId = null;
        // Get the main output worksheet.
        this.workbook = new XSSFWorkbook();
        this.worksheet = this.workbook.createSheet("Changes");
        // Get a data formatter.
        DataFormat format = this.workbook.createDataFormat();
        short fmt = format.getFormat("###,##0");
        // Create the header style.
        this.headStyle = this.workbook.createCellStyle();
        this.headStyle.setFillForegroundColor(IndexedColors.GREY_25_PERCENT.getIndex());
        this.headStyle.setFillPattern(FillPatternType.SOLID_FOREGROUND);
        // Create the number style.
        this.numStyle = this.workbook.createCellStyle();
        this.numStyle.setDataFormat(fmt);
        this.numStyle.setAlignment(HorizontalAlignment.RIGHT);
        // Create the flag style.
        this.flagStyle = this.workbook.createCellStyle();
        this.flagStyle.setAlignment(HorizontalAlignment.CENTER);
        // Create the group style.
        this.groupStyle = this.workbook.createCellStyle();
        this.groupStyle.setWrapText(true);
        // Initialize the row number.
        this.rowNum = 0;
        // Get the special genome IDs.
        this.specials = processor.getSpecial();
        // Initialize the set of all genome IDs to empty.
        this.all = new TreeSet<String>();
        // Save the command processor and the output stream.
        this.processor = processor;
        this.outStream = output;
    }

    @Override
    protected void registerGenome(Genome genome) {
        if (this.baseId == null) {
            // Here we have the base genome, so save its ID.
            this.baseId = genome.getId();
        } else {
            // Here we have a mutant genome, so add its ID to the all-genome set.
            this.all.add(genome.getId());
        }
    }

    @Override
    protected void openReport(List<GenomeLabel> genomeLabels) {
        // Fill in the headings.
        this.addRow();
        this.setTextCell(0, "fig_id", this.headStyle);
        this.setTextCell(1, "start_loc", this.headStyle);
        this.setTextCell(2, "stop_loc", this.headStyle);
        this.setTextCell(3, "strand", this.headStyle);
        this.setTextCell(4, "gene_name", this.headStyle);
        this.setTextCell(5, "length", this.headStyle);
        this.setTextCell(6, "function", this.headStyle);
        this.setTextCell(7, "groups", this.headStyle);
        this.setTextCell(8, "all", this.headStyle);
        // Fix the text column widths.
        this.worksheet.setColumnWidth(6, TEXT_WIDTH);
        this.worksheet.setColumnWidth(7, TEXT_WIDTH);
    }

    @Override
    protected void openAlignment(String title, RegionList regions, Feature feat) {
        // Denote that no genomes have significant changes in this alignment.
        this.changed = new TreeSet<String>();
        // Save the base-genome feature.
        this.feat = feat;
        ExtendedProteinRegion baseRegion = regions.get(feat.getId());
        // Loop through the regions, finding changes.
        for (ExtendedProteinRegion region : regions) {
            if (region != baseRegion) {
                Feature feat2 = region.getFeature();
                if (this.testRegions(baseRegion, region)) {
                    // Here we have a significant change.
                    this.changed.add(feat2.getParent().getId());
                }
            }
        }
    }

    @Override
    protected void processSnips(SnipColumn snipCol) {
    }

    @Override
    protected void closeAlignment() {
        // Output this feature if every genome in the subset had a significant change.
        if (this.changed.containsAll(this.specials)) {
            // Determine whether or not it had significant changes in all the mutants.
            String flag = (this.changed.containsAll(this.all) ? "X" : "");
            // Get the feature location.
            Location loc = this.feat.getLocation();
            // Find the gene name.
            String gene = "";
            for (String alias : this.feat.getAliases()) {
                if (GENE_NAME.matcher(alias).matches())
                    gene = alias;
            }
            // Form the group list.
            String groups = StringUtils.join(this.processor.getGroups(this.feat.getId()), " | ");
            Collection<String> subs = this.feat.getSubsystems();
            if (subs.size() > 0) {
                if (groups != null)
                    groups += " | ";
                groups += StringUtils.join(subs, " | ");
            } else if (groups == null)
                groups = "";
            // Now create the row.
            this.addRow();
            this.setTextCell(0, this.feat.getId(), null);
            this.setNumCell(1, loc.getBegin());
            this.setNumCell(2, loc.getEnd());
            this.setTextCell(3, loc.getStrand(), this.flagStyle);
            this.setTextCell(4, gene, null);
            this.setNumCell(5, loc.getLength());
            this.setTextCell(6, this.feat.getPegFunction(), this.groupStyle);
            this.setTextCell(7, groups, this.groupStyle);
            this.setTextCell(8, flag, this.flagStyle);
        }

    }

    @Override
    protected void closeReport() {
        // Fix the column widths.
        for (int i = 0; i < 6; i++)
            this.worksheet.autoSizeColumn(i);
        // Write the spreadsheet.
        try {
            this.workbook.write(this.outStream);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    /**
     * Add a new row to the spreadsheet.
     */
    protected void addRow() {
        this.ssRow = this.worksheet.createRow(this.rowNum);
        this.rowNum++;
    }

    /**
     * Create a text cell in the current row.
     *
     * @param i				column number of new cell
     * @param string		content of the cell
     * @param style			style to give the cell, or NULL to use the default
     *
     * @return the created cell
     */
    protected Cell setTextCell(int i, String string, CellStyle style) {
        Cell retVal = this.ssRow.createCell(i);
        retVal.setCellValue(string);
        if (style != null)
            retVal.setCellStyle(style);
        return retVal;
    }

    /**
     * Create a number cell in the current row.
     *
     * @param i				column number of new cell
     * @param num			content of the cell
     *
     * @return the created cell
     */
    protected Cell setNumCell(int i, int num) {
        Cell retVal = this.ssRow.createCell(i);
        retVal.setCellValue((double) num);
        retVal.setCellStyle(this.numStyle);
        return retVal;
    }

    /**
     * Compare two regions to detect a change.
     *
     * @param region1	base region
     * @param region2	test region
     *
     * @return TRUE if there is a reportable change, else FALSE
     */
    protected abstract boolean testRegions(ExtendedProteinRegion region1, ExtendedProteinRegion region2);

    /**
     * @return TRUE if the two sequence strings represent a significant change, else FALSE.
     *
     * @param seq1		base sequence
     * @param seq2		other sequence
     */
    public boolean isSignificant(String seq1, String seq2) {
        return ! StringUtils.endsWith(seq1, seq2) && ! StringUtils.endsWith(seq2, seq1);
    }

    /**
     * This nested class produces a major-change report focused on protein changes.
     */
    public static class Protein extends MajorSnipReporter {

        public Protein(OutputStream output, IParms processor) {
            super(output, processor);
        }

        @Override
        protected boolean testRegions(ExtendedProteinRegion region1, ExtendedProteinRegion region2) {
            String prot1 = region1.getProteinTranslation();
            String prot2 = region2.getProteinTranslation();
            return super.isSignificant(prot1, prot2);
        }

    }

    /**
     * This nested class produces a major-change report focused on upstream changes.
     */
    public static class Upstream extends MajorSnipReporter {

        public Upstream(OutputStream output, IParms processor) {
            super(output, processor);
        }

        @Override
        protected boolean testRegions(ExtendedProteinRegion region1, ExtendedProteinRegion region2) {
            String dna1 = region1.getUpstreamDna();
            String dna2 = region2.getUpstreamDna();
            return super.isSignificant(dna1, dna2);
        }

    }




}
