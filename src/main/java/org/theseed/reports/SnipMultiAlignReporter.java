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

import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.sequence.Sequence;

/**
 * This report lists the snips from each alignment.  It should only be used in cases where the alignments are fairly close, because the size
 * of the report can otherwise get very large.
 *
 * @author Bruce Parrello
 *
 */
public class SnipMultiAlignReporter extends MultiAlignReporter {

    // FIELDS
    /** ID of the base genomes */
    private Set<String> baseGenomeIds;
    /** ID of the first base genome */
    private String base0GenomeId;

    public SnipMultiAlignReporter(OutputStream outStream) {
        super(outStream);
    }

    @Override
    public void openReport(Genome genome, String[] altBases) {
        this.base0GenomeId = genome.getId();
        // We use a tree set because we expect only one or two elements at most.
        this.baseGenomeIds = new TreeSet<String>();
        this.baseGenomeIds.add(this.base0GenomeId);
        this.baseGenomeIds.addAll(Arrays.asList(altBases));
        // Start the report.
        this.println("function\tpeg\tsnip\toriginal\tlocation");
    }

    @Override
    public void writeAlignment(String title, List<Sequence> alignment) {
        // We need to separate the sequences into base sequences and variable sequences.  The first base
        // sequence is special.
        String base0Sequence = null;
        List<String> bases = new ArrayList<String>(this.baseGenomeIds.size());
        List<Sequence> others = new ArrayList<Sequence>(alignment.size() - this.baseGenomeIds.size());
        for (int i = 0; i < alignment.size(); i++) {
            Sequence curr = alignment.get(i);
            String thisGenome = Feature.genomeOf(curr.getLabel());
            if (! this.baseGenomeIds.contains(thisGenome))
                others.add(curr);
            else {
                String seq = curr.getSequence();
                bases.add(seq);
                if (thisGenome.contentEquals(this.base0GenomeId))
                    base0Sequence = seq;
            }
        }
        // Now we process each sequence other than the base, building snips.
        int width = bases.get(0).length();
        for (Sequence curr : others) {
            // This will track our current offset in the current sequence.
            int offset = 0;
            // Get the current sequence string.
            String currSequence = curr.getSequence();
            // Get the location of the aligned sequence.
            Location loc = Location.fromString(curr.getComment());
            // Loop through the sequence.
            int p = 0;
            while (p < width) {
                char c = currSequence.charAt(p);
                if (difference(c, bases, p)) {
                    // Here we have a snip.  We need to find the end of it.
                    int offsetIn = offset;
                    int pIn = p;
                    p++; if (c != '-') offset++;
                    while (p < width && difference(currSequence.charAt(p), bases, p)) {
                        if (currSequence.charAt(p) != '-') offset++;
                        p++;
                    }
                    // Output the full snip.
                    String currSnip = currSequence.substring(pIn, p);
                    String baseSnip = base0Sequence.substring(pIn, p);
                    Location snipLoc = loc.subLocation(offsetIn, offset - offsetIn);
                    String locString;
                    if (snipLoc.getLength() > 0)
                        locString = snipLoc.toString();
                    else {
                        // Here we have a zero-length location, that is, it is entirely gap characters.  Create a bogus length-1
                        // version, then forge a fake string out of it.
                        snipLoc = loc.subLocation(offsetIn, 1);
                        locString = String.format("%s%c%d", loc.getContigId(), loc.getDir(), loc.getBegin());
                    }
                    this.print("%s\t%s\t%s\t%s\t%s", title, curr.getLabel(), currSnip, baseSnip, locString);
                } else {
                    // No snip, keep going.
                    p++; if (c != '-') offset++;
                }
            }
        }
    }

    /**
     * @return TRUE if the input character does not match any of the base sequence characters at the same position
     *
     * @param c		character to check
     * @param bases	list of sequence strings to check against
     * @param p		position to check in the sequence strings
     */
    private boolean difference(char c, List<String> bases, int p) {
        return ! bases.stream().anyMatch(x -> x.charAt(p) == c);
    }

    @Override
    public void closeReport() {
    }

}
