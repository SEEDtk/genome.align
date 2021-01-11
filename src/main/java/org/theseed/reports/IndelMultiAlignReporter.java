/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.sequence.Sequence;

/**
 * This report looks for indels.  If the alignment contains more than one row per genome, the alignment is ignored.  If there is no long indel
 * sequence, it is ignored.  The alternate bases will be eliminated from the output.
 *
 * @author Bruce Parrello
 *
 */
public class IndelMultiAlignReporter extends MultiAlignReporter {

    // FIELDS
    /** minimum indel length */
    private static final String MIN_INDEL = StringUtils.repeat('-', 20);
    /** ID of the alternate base genomes */
    private Set<String> altGenomeIds;
    /** ID of the first base genome */
    private String base0GenomeId;

    public IndelMultiAlignReporter(OutputStream outStream) {
        super(outStream);
    }

    @Override
    public void openReport(Genome genome, String[] altBases) {
        // Save the base genome ID.
        this.base0GenomeId = genome.getId();
        // Save the alternate base IDs.
        this.altGenomeIds = new TreeSet<String>();
        this.altGenomeIds.addAll(Arrays.asList(altBases));
    }

    @Override
    public void writeAlignment(String title, List<Sequence> alignment) {
        // Isolate the aligned sequences (that is, the ones that aren't bases.
        List<Sequence> aligned = new ArrayList<Sequence>(alignment.size());
        // This will hold the base sequence.
        Sequence base = null;
        // This will count the number of aligned sequences with long indels.
        int indelCount = 0;
        // This will be TRUE if the base has a long indel.
        boolean indelBase = false;
        // These track if we've found more than one sequence per genome.
        boolean dups = false;
        Set<String> foundGenomes = new TreeSet<String>();
        // Loop through the sequences of the alignment.
        for (Sequence seq : alignment) {
            String genomeId = Feature.genomeOf(seq.getLabel());
            if (genomeId.contentEquals(this.base0GenomeId)) {
                if (base != null)
                    dups = true;
                else {
                    base = seq;
                    indelBase = StringUtils.contains(seq.getSequence(), MIN_INDEL);
                }
            } else if (! this.altGenomeIds.contains(genomeId)) {
                if (foundGenomes.contains(genomeId))
                    dups = true;
                else {
                    aligned.add(seq);
                    if (StringUtils.contains(seq.getSequence(), MIN_INDEL))
                        indelCount++;
                    foundGenomes.add(genomeId);
                }
            }
        }
        if (! dups && (indelBase && indelCount < aligned.size() || ! indelBase && indelCount > 0)) {
            // Here we want to print the alignment, with the base on top.
            this.println(title);
            this.println();
            showSequence(base);
            for (Sequence seq : aligned)
                showSequence(seq);
            this.println();
        }
    }

    /**
     * Display a single sequence in an alignment.
     *
     * @param seq		sequence to display
     */
    protected void showSequence(Sequence seq) {
        this.print("%s\t%s\t%s", seq.getLabel(), seq.getComment(), seq.getSequence());
    }

    @Override
    public void closeReport() {
    }

}
