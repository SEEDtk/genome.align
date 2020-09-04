/**
 *
 */
package org.theseed.genome.align;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.proteins.Function;
import org.theseed.proteins.FunctionMap;
import org.theseed.reports.MultiAlignReporter;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.clustal.ClustalPipeline;
import org.theseed.utils.BaseProcessor;

/**
 * This command creates alignments for a list of genomes.  The DNA sequences will be organized by functional assignment and then
 * aligned.  Only functions that occur at least three times will participate in the alignment, and only protein-encoding genes
 * will be considered.  At least one sequence must be in the first genome, as well.
 *
 * The positional parameters are the name of the file containing the first genome, then the names of the files containing the other
 * genomes.  The alignments will be written to the standard output.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	show more detailed progress messages
 * -m	maximum kmer distance for a feature to be aligned
 * -K	kmer size for computing distances
 *
 * --format		output report format
 * --workDir	working directory name; the default is "Temp" in the current directory
 *
 * @author Bruce Parrello
 *
 */
public class GtoAlignProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GtoAlignProcessor.class);
    /** role map of functional assignments */
    private FunctionMap functionMap;
    /** hash of function IDs to sequence lists */
    private Map<String, SequenceList> sequenceMap;
    /** base genome */
    private Genome baseGenome;

    // COMMAND-LINE OPTIONS

    /** output format */
    @Option(name = "--format", usage = "output format")
    private MultiAlignReporter.Type outputFormat;

    /** working directory */
    @Option(name = "--workDir", metaVar = "Temp", usage = "working directory for temporary files")
    private File workDir;

    /** minimum kmer distance for a feature to be aligned */
    @Option(name = "-m", aliases = { "--maxDist" }, metaVar = "0.5", usage = "maximum kmer distance for a sequence to be placed in an alignment")
    private double maxDist;

    /** kmer size for computing distances */
    @Option(name = "-K", metaVar = "15", usage = "kmer size for computing sequence distances")
    private int kmerSize;

    /** first genome */
    @Argument(index = 0, metaVar = "baseGTO", usage = "primary genome file", required = true)
    private File gtoBaseFile;

    /** comparative genomes */
    @Argument(index = 1, metaVar = "gto2 gto3 ...", usage = "genomes to align to primary", multiValued = true, required = true)
    private List<File> gtoFiles;

    @Override
    protected void setDefaults() {
        this.workDir = new File(System.getProperty("user.dir"), "Temp");
        this.kmerSize = DnaKmers.kmerSize();
        this.maxDist = 0.6;
        this.outputFormat = MultiAlignReporter.Type.TEXT;
    }

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
        // Verify the genome files.
        if (! this.gtoBaseFile.canRead())
            throw new FileNotFoundException("Genome file " + this.gtoBaseFile + " not found or unreadable.");
        for (File gtoFile : this.gtoFiles) {
            if (! gtoFile.canRead())
                throw new FileNotFoundException("Genome file " + gtoFile + " not found or unreadable.");
        }
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Create the temporary FASTA file.
        File tempFile = File.createTempFile("align", ".fna", this.workDir);
        tempFile.deleteOnExit();
        // Create the function maps.
        log.info("Initializing function maps.");
        this.functionMap = new FunctionMap();
        this.sequenceMap = new HashMap<String, SequenceList>(3000);
        // Process the base genome to compute the functions of interest.
        this.processBase(this.gtoBaseFile);
        // Loop through the other genomes.
        for (File gtoFile : this.gtoFiles) {
            Genome genome = new Genome(gtoFile);
            log.info("Scanning genome {}.", genome);
            // Loop through the genome's pegs.
            int kept = 0;
            for (Feature feat : genome.getPegs()) {
                // Get the function and look for a sequence list.
                String function = feat.getFunction();
                Function fun = this.functionMap.getByName(function);
                if (fun != null) {
                    SequenceList seqs = this.sequenceMap.get(fun.getId());
                    if (seqs != null) {
                        if (this.addFeature(seqs, feat))
                            kept++;
                    }
                }
            }
            log.info("{} sequences kept from {}.", kept, genome);
        }
        // Process the alignments.
        try (MultiAlignReporter reporter = this.outputFormat.create(System.out)) {
            // Initialize the output report.
            reporter.openReport(this.baseGenome);
            int alignCount = 0;
            // Loop through the alignments
            for (Map.Entry<String, SequenceList> alignRequest : this.sequenceMap.entrySet()) {
                // Verify that this alignment is big enough and has variations.
                SequenceList seqs = alignRequest.getValue();
                if (seqs.size() >= 3 && seqs.getMaxDist() > 0.0) {
                    String function = this.functionMap.getName(alignRequest.getKey());
                    log.info("Processing alignment for {}.", function);
                    // Save the sequences.
                    seqs.save(tempFile);
                    // Process the alignment.
                    ClustalPipeline aligner = new ClustalPipeline(tempFile);
                    List<Sequence> alignment = aligner.run();
                    reporter.writeAlignment(function, alignment);
                    alignCount++;
                }
            }
            // Finish the report.
            reporter.closeReport();
            log.info("{} alignments output.", alignCount);
        }

    }

    /**
     * Process the base genome.  This includes isolating all the non-hypothetical functions and storing the
     * DNA sequences in the sequence lists.
     *
     * @param gtoFile	file containing the base genome
     *
     * @throws IOException
     */
    private void processBase(File gtoFile) throws IOException {
        // This will count the total features kept.
        int kept = 0;
        // Read in the genome.
        Genome genome = new Genome(gtoFile);
        log.info("Scanning features from {}.", genome);
        this.baseGenome = genome;
        // Loop through the pegs.
        for (Feature feat : genome.getPegs()) {
            String function = feat.getFunction();
            if (! function.contentEquals("hypothetical protein")) {
                String funId = this.functionMap.findOrInsert(function).getId();
                SequenceList seqList = this.sequenceMap.get(funId);
                if (seqList != null) {
                    // Here we have an existing function, so add this feature to the sequence list.
                    if (this.addFeature(seqList, feat)) kept++;
                } else {
                    // Here we need to create a new sequence list.
                    Location loc = feat.getLocation();
                    seqList = new SequenceList(feat.getId(), loc.toString(), genome.getDna(loc));
                    this.sequenceMap.put(funId, seqList);
                    kept++;
                }
            }
        }
        log.info("{} functions found in {}, {} features kept.", this.sequenceMap.size(), genome, kept);
    }

    /**
     * Add a new feature to a sequence list.
     *
     * @param seqList	target sequence list
     * @param feat		feature whose sequence should be added
     *
     * @return TRUE if the new feature was eligible, else FALSE
     */
    private boolean addFeature(SequenceList seqList, Feature feat) {
        Location loc = feat.getLocation();
        String dna = feat.getParent().getDna(loc);
        boolean retVal = seqList.add(feat.getId(), loc.toString(), dna, this.maxDist);
        return retVal;
    }
}
