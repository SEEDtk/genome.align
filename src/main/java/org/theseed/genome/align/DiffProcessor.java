/**
 *
 */
package org.theseed.genome.align;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.FeatureFilter;
import org.theseed.proteins.Function;
import org.theseed.proteins.FunctionMap;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command measures how many useful snip changes there are between close genomes.  The proteins in a genome are aligned by
 * role name.  A useful snip change occurs if two non-hypothetical proteins have different protein sequences and neither is a
 * substring of the other.  (This last restriction is to remove differences caused by truncation or a change in the start
 * call.  For each genome other than the base, the count of useful snip changes is output.
 *
 * The positional parameters are the file name of the base genome followed by the file names of the other genomes to look at.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 *
 * --filter		feature filter types (multiple)
 *
 * @author Bruce Parrello
 *
 */
public class DiffProcessor extends BaseProcessor implements FeatureFilter.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(DiffProcessor.class);
    /** function ID map */
    private FunctionMap funMap;
    /** directory of base genome feature proteins by function ID */
    private Map<String, List<Feature>> featureMap;
    /** list of feature filters */
    private List<FeatureFilter> filters;

    // COMMAND-LINE OPTIONS

    /** alternate base genomes */
    @Option(name = "--alt", usage = "alternate base genome")
    private List<File> altGenomeFiles;

    /** filtering schemes for features */
    @Option(name = "--filter", usage = "type of feature filtering")
    private List<FeatureFilter.Type> filterTypes;

    /** feature-filter file for LIST filtering */
    @Option(name = "--fidFile", metaVar = "fidsToKeep.tbl", usage = "file containing list of acceptable features for LIST filter")
    private File fidFile;

    /** base genome file */
    @Argument(index = 0, metaVar = "base.gto", usage = "base genome file", required = true)
    private File baseFile;

    /** genome files to examine for snip changes */
    @Argument(index = 1, metaVar = "test1.gto test2.gto ...", usage = "GTOs to test for snip changes", required = true)
    private List<File> testFiles;

    @Override
    protected void setDefaults() {
        this.altGenomeFiles = new ArrayList<File>();
        this.filterTypes = new ArrayList<FeatureFilter.Type>();
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Verify that the alternate genomes exist.
        for (File altFile : this.altGenomeFiles) {
            if (! altFile.canRead())
                throw new FileNotFoundException("Alternate base genome file " + altFile + " not found or unreadable.");
        }
        // Verify that the base genome exists.
        if (! this.baseFile.canRead())
            throw new FileNotFoundException("Base genome file " + this.baseFile + " not found or unreadable.");
        // Verify that the test genomes exist.
        for (File testFile : this.testFiles) {
            if (! testFile.canRead())
                throw new FileNotFoundException("Test genome file " + testFile + " not found or unreadable.");
        }
        // Build the filters.
        this.filters = new ArrayList<FeatureFilter>(this.filterTypes.size());
        for (FeatureFilter.Type filterType : this.filterTypes)
            this.filters.add(filterType.create(this));
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Parse the base genome and create the function map.
        this.parseBaseGenome();
        // Parse the alternate genomes as well.
        for (File altFile : this.altGenomeFiles)
            this.parseAltGenome(altFile);
        // Write the output header.
        System.out.println("genome_id\tgenome_name\tchanges");
        // Loop through the test genomes.
        for (File testFile : this.testFiles) {
            log.info("Loading test genome from {}.", testFile);
            Genome genome = new Genome(testFile);
            int count = 0;
            int changes = 0;
            for (Feature feat : genome.getPegs()) {
                // Get the function's feature list.
                List<Feature> featList = this.getFeatureList(feat);
                if (featList != null) {
                    count++;
                    // This feature is a change if it does not match ANY of the proteins in the feature list.
                    String protein = feat.getProteinTranslation();
                    boolean changed = featList.stream().noneMatch(x -> this.protMatch(x, protein));
                    if (changed) changes++;
                }
            }
            log.info("{} of {} proteins changed in {}.", changes, count, genome);
            System.out.format("%s\t%s\t%d%n", genome.getId(), genome.getName(), changes);
        }
    }

    /**
     * @return TRUE if the protein in the specified feature is different from the specified protein string
     *
     * @param feat		feature of interest
     * @param protein	protein string for comparison
     */
    private boolean protMatch(Feature feat, String protein) {
        String fProtein = feat.getProteinTranslation();
        boolean retVal = (fProtein.endsWith(protein) || protein.endsWith(fProtein));
        return retVal;
    }

    /**
     * Load the base genome and run through its features, building the function map and assigning each feature to
     * its function.
     *
     * @throws IOException
     */
    private void parseBaseGenome() throws IOException {
        log.info("Processing base genome in {}.", this.baseFile);
        Genome baseGenome = new Genome(this.baseFile);
        this.funMap = new FunctionMap();
        this.featureMap = new HashMap<String, List<Feature>>(4000);
        int count = 0;
        for (Feature feat : baseGenome.getPegs()) {
            // Apply the filters.  All must accept the feature.
            if (this.filters.stream().allMatch(x -> x.filter(feat))) {
                String pegFunction = feat.getPegFunction();
                // Skip hypothetical or missing proteins.
                if (! Feature.isHypothetical(pegFunction) && feat.getProteinTranslation() != null) {
                    // Compute the function ID and put this feature in the function -> feature map.
                    Function fun = this.funMap.findOrInsert(pegFunction);
                    List<Feature> featList = this.featureMap.computeIfAbsent(fun.getId(), x -> new ArrayList<Feature>());
                    featList.add(feat);
                    count++;
                }
            }
        }
        log.info("{} features in {} functions found in {}.", count, this.featureMap.size(), baseGenome);
    }

    /**
     * Load an alternate base genome and add its features to the feature map.  Unlike the base genome, we don't
     * create new functions this time.
     *
     * @param altFile	file containing the alternate base genome
     *
     * @throws IOException
     */
    private void parseAltGenome(File altFile) throws IOException {
        log.info("Processing alternate base genome in {}.", altFile);
        Genome altGenome = new Genome(altFile);
        int count = 0;
        for (Feature feat : altGenome.getPegs()) {
            // Get the function's feature list.
            List<Feature> featList = this.getFeatureList(feat);
            if (featList != null && feat.getProteinTranslation() != null) {
                // Here the function was found in the base genome.  Add our feature to it.
                featList.add(feat);
                count++;
            }
        }
        log.info("{} features found in {}.", count, altGenome);
    }

    /**
     * @return the feature list for the specified feature's function, or NULL if it is not in the base genome
     *
     * @param feat	feature whose function should be interrogated
     */
    private List<Feature> getFeatureList(Feature feat) {
        List<Feature> retVal = null;
        String pegFunction = feat.getPegFunction();
        Function fun = this.funMap.getByName(pegFunction);
        if (fun != null)
            retVal = this.featureMap.get(fun.getId());
        return retVal;
    }


    @Override
    public Set<String> getFeatureSet() throws ParseFailureException, IOException {
        if (this.fidFile == null)
            throw new ParseFailureException("Feature-list file required for filtering.");
        if (! this.fidFile.canRead())
            throw new FileNotFoundException("Feature-list file is not found or unreadable.");
        // Read the set of feature IDs from column 1 of the file
        var retVal = TabbedLineReader.readSet(this.fidFile, "1");
        return retVal;
    }

}
