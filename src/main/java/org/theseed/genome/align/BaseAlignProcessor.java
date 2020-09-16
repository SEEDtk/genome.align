/**
 *
 */
package org.theseed.genome.align;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.sequence.DnaKmers;
import org.theseed.utils.BaseProcessor;

/**
 * @author Bruce Parrello
 *
 */
public abstract class BaseAlignProcessor extends BaseProcessor {

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BaseAlignProcessor.class);
    /** working directory */
    @Option(name = "--workDir", metaVar = "Temp", usage = "working directory for temporary files")
    private File workDir;
    /** minimum kmer distance for a feature to be aligned */
    @Option(name = "-m", aliases = {
            "--maxDist" }, metaVar = "0.5", usage = "maximum kmer distance for a sequence to be placed in an alignment")
    private double maxDist;
    /** kmer size for computing distances */
    @Option(name = "-K", metaVar = "15", usage = "kmer size for computing sequence distances")
    private int kmerSize;

    @Override
    protected void setDefaults() {
        this.workDir = new File(System.getProperty("user.dir"), "Temp");
        this.kmerSize = DnaKmers.kmerSize();
        this.maxDist = 0.6;
        setProcessDefaults();
    }

    /**
     * Set the defaults for the parameters local to the processor.
     */
    protected abstract void setProcessDefaults();

    @Override
    protected boolean validateParms() throws IOException {
        // Verify the work directory.
        if (! this.workDir.exists()) {
            log.info("Creating work directory {}.", this.workDir);
            FileUtils.forceMkdir(this.workDir);
            this.workDir.deleteOnExit();
        } else if (! this.workDir.isDirectory())
            throw new FileNotFoundException("Work directory " + this.workDir + " not found or invalid.");
        // Validate the maximum distance.
        if (this.maxDist <= 0.0)
            throw new IllegalArgumentException("Maximum distance " + this.maxDist + " must be greater than 0.");
        // Process the kmer size.
        if (this.kmerSize < 3 || this.kmerSize > 100)
            throw new IllegalArgumentException("Kmer size " + this.kmerSize + " is out of range.  Must be >=3 and <= 100.");
        DnaKmers.setKmerSize(this.kmerSize);
        this.validateProcessParms();
        return true;
    }

    /**
     * Validate the parameters local to the processor.
     *
     * @throws IOException
     */
    protected abstract void validateProcessParms() throws IOException;

    /**
     * @return the working directory
     */
    public File getWorkDir() {
        return this.workDir;
    }

    /**
     * @return the maximum acceptable kmer distance
     */
    public double getMaxDist() {
        return this.maxDist;
    }

}
