/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import java.util.List;

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
    /** ID of the base genome */
    private String baseGenomeId;

    public SnipMultiAlignReporter(OutputStream outStream) {
        super(outStream);
    }

    @Override
    public void openReport(Genome genome) {
        this.baseGenomeId = genome.getId();
        this.println("function\tpeg\tsnip\toriginal\tlocation");

    }

    @Override
    public void writeAlignment(String title, List<Sequence> alignment) {
        // Find the base genome's sequence.
        Sequence base = null;
        for (int i = 0; base == null && i < alignment.size(); i++) {
            if (Feature.genomeOf(alignment.get(i).getLabel()).contentEquals(this.baseGenomeId))
                base = alignment.get(i);
        }
        // Now we process each sequence other than the base, building snips.
        String baseSequence = base.getSequence();
        int width = baseSequence.length();
        for (Sequence curr : alignment) {
            if (curr != base) {
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
                    if (c != baseSequence.charAt(p)) {
                        // Here we have a snip.  We need to find the end of it.
                        int offsetIn = offset;
                        int pIn = p;
                        p++; if (c != '-') offset++;
                        while (p < width && currSequence.charAt(p) != baseSequence.charAt(p)) {
                            if (currSequence.charAt(p) != '-') offset++;
                            p++;
                        }
                        // Output the full snip.
                        String currSnip = currSequence.substring(pIn, p);
                        String baseSnip = baseSequence.substring(pIn, p);
                        Location snipLoc = loc.subLocation(offsetIn, offset - offsetIn);
                        this.print("%s\t%s\t%s\t%s\t%s", title, curr.getLabel(), currSnip, baseSnip, snipLoc.toString());
                    } else {
                        // No snip, keep going.
                        p++; if (c != '-') offset++;
                    }
                }
            }
        }
    }

    @Override
    public void closeReport() {
    }

}
