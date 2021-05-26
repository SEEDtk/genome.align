/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.util.List;

import org.theseed.genome.Genome;
import org.theseed.sequence.Sequence;

/**
 * This is the base class for a multiple-alignment report.
 *
 * @author Bruce Parrello
 *
 */
public abstract class MultiAlignReporter implements AutoCloseable {

    // FIELDS
    /** output file */
    private File outFile;
    /** current output writer, or NULL if none */
    private PrintWriter writer;

    /**
     * Construct a report object for a specified output stream.
     *
     * @param outFile		output file to receive the report
     */
    public MultiAlignReporter(File outFile) {
        this.outFile = outFile;
        this.writer = null;
    }

    /**
     * Initialize the report.
     *
     * @param genome	base genome
     * @param altBases 	array of IDs for additional base genomes
     */
    public abstract void openReport(Genome genome, String[] altBases);

    /**
     * Register a genome used in the report (optional).
     *
     * @param genome	genome to register
     */
    public void registerGenome(Genome genome) { }

    /**
     * Output a particular alignment.
     *
     * @param fid		base genome feature ID
     * @param title		name of the alignment
     * @param alignment	aligned sequences
     */
    public abstract void writeAlignment(String fid, String title, List<Sequence> alignment);

    /**
     * Finish the report.
     */
    public abstract void closeReport();

    /**
     * Write a formatted output line.
     */
    protected void print(String format, Object... args) {
        this.getWriter().format(format, args);
        this.getWriter().println();
    }

    /**
     * Write an unformatted output line.
     */
    protected void println(String line) {
        this.getWriter().println(line);
    }

    /**
     * Write a blank output line.
     */
    protected void println() {
        this.getWriter().println();
    }

    /**
     * @return the writer object
     */
    protected PrintWriter getWriter() {
        if (this.writer == null) {
            try {
                this.writer = new PrintWriter(this.outFile);
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
        return writer;
    }

    @Override
    public void close() {
        // Close according to the type of output file.
        if (this.writer != null) {
               this.writer.close();
        }
    }

    /**
     * Enum for report types
     */
    public static enum Type {
            TEXT, HTML, SNIPS, INDELS;

        public MultiAlignReporter create(File outFile) {
            MultiAlignReporter retVal = null;
            switch (this) {
            case TEXT :
                retVal = new TextMultiAlignReporter(outFile);
                break;
            case HTML :
                retVal = new HtmlMultiAlignReporter(outFile);
                break;
            case SNIPS :
                retVal = new SnipMultiAlignReporter(outFile);
                break;
            case INDELS :
                retVal = new IndelMultiAlignReporter(outFile);
                break;
            }
            return retVal;
        }
    }

    /**
     * @return the output file
     */
    public File getOutFile() {
        return this.outFile;
    }

}
