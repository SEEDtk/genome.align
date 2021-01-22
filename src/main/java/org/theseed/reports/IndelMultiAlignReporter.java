/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
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
    /** map of genome IDs to genomes */
    private Map<String, Genome> genomeMap;

    public IndelMultiAlignReporter(OutputStream outStream) {
        super(outStream);
        this.genomeMap = new HashMap<String, Genome>();
    }

    @Override
    public void openReport(Genome genome, String[] altBases) {
        // Save the base genome ID.
        this.base0GenomeId = genome.getId();
        // Add it to the genome map.
        this.genomeMap.put(this.base0GenomeId, genome);
        // Save the alternate base IDs.
        this.altGenomeIds = new TreeSet<String>();
        this.altGenomeIds.addAll(Arrays.asList(altBases));
    }

    @Override
    public void registerGenome(Genome genome) {
        this.genomeMap.put(genome.getId(), genome);
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
                    // Fix the indels at either end of the sequence.
                    this.fixIndels(genomeId, seq);
                }
            } else if (! this.altGenomeIds.contains(genomeId)) {
                if (foundGenomes.contains(genomeId))
                    dups = true;
                else {
                    aligned.add(seq);
                    if (StringUtils.contains(seq.getSequence(), MIN_INDEL))
                        indelCount++;
                    foundGenomes.add(genomeId);
                    // Fix the indels at either end of the sequence.
                    this.fixIndels(genomeId, seq);
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
     * This method will convert the sequence to include the upstream and downstream regions for the indels.
     *
     * @param genomeId	ID of the target genome
     * @param seq		sequence being processed
     */
    private void fixIndels(String genomeId, Sequence seq) {
        // Get the genome containing this sequence.
        Genome genome = this.genomeMap.get(genomeId);
        // Convert the current sequence to upper case.
        String sequence = seq.getSequence().toUpperCase();
        // Get the DNA location.
        Location loc = Location.fromString(seq.getComment());
        // Check for starting and ending indels.
        int start = 0;
        while (sequence.charAt(start) == '-') start++;
        if (start > 0) {
            Location loc2 = loc.upstream(start);
            String dna = StringUtils.leftPad(genome.getDna(loc2), start, '-');
            sequence = dna.toLowerCase() + StringUtils.substring(sequence, start);
        }
        int last = sequence.length() - 1;
        int end = last;
        while (sequence.charAt(end) == '-') end--;
        if (end < last) {
            Location loc2 = loc.downstream(last - end);
            String dna = StringUtils.rightPad(genome.getDna(loc2), loc2.getLength(), '-');
            sequence = StringUtils.substring(sequence, 0, end + 1) + dna.toLowerCase();
        }
        // Update the sequence.
        seq.setSequence(sequence);
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
