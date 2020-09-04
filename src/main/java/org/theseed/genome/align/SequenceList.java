/**
 *
 */
package org.theseed.genome.align;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;

/**
 * This class represents a list of sequences.  The first sequence is converted into DnaKmers so that it may be used to
 * determine distance.
 *
 * @author Bruce Parrello
 */
public class SequenceList {

    // FIELDS
    /** DNA kmers for base sequence */
    private DnaKmers baseKmers;
    /** maximum distance found */
    private double maxDist;
    /** list of attached sequences (including the base) */
    private List<Sequence> members;

    /**
     * Construct a sequence list from a base sequence.
     *
     * @param id		ID of the sequence
     * @param comment	sequence comment
     * @param dna		sequence DNA
     */
    public SequenceList(String id, String comment, String dna) {
        // Create the kmers for testing distance.
        this.baseKmers = new DnaKmers(dna);
        // Add the sequence to the sequence list.
        this.members = new ArrayList<Sequence>();
        this.add(id, comment, dna);
        // Denote we only have one sequence.
        this.maxDist = 0.0;
    }

    /**
     * Add a sequence to this list.
     *
     * @param id		ID of the sequence
     * @param comment	sequence comment
     * @param dna		sequence DNA
     */
    protected void add(String id, String comment, String dna) {
        Sequence seq = new Sequence(id, comment, dna);
        this.members.add(seq);
    }

    /**
     * Add a sequence to this list if it is within a specified distance.
     *
     * @param id		ID of the sequence
     * @param comment	sequence comment
     * @param dna		sequence DNA
     * @param max		maximum acceptable distance
     *
     * @return TRUE if the sequence was added, else FALSE
     */
    public boolean add(String id, String comment, String dna, double max) {
        double distance = this.distance(dna);
        boolean retVal = (distance <= max);
        if (retVal) {
            this.add(id, comment, dna);
            if (distance > this.maxDist)
                this.maxDist = distance;
        }
        return retVal;
    }

    /**
     * @return the number of sequences stored
     */
    public int size() {
        return this.members.size();
    }

    /**
     * Save the sequences to a file.
     *
     * @param fileName	name of the output file
     *
     * @throws IOException
     */
    public void save(File fileName) throws IOException {
        try (FastaOutputStream outStream = new FastaOutputStream(fileName)) {
            outStream.write(this.members);
        }
    }

    /**
     * @return the distance to a proposed sequence
     *
     * @param dna	DNA sequence to test
     */
    public double distance(String dna) {
        DnaKmers other = new DnaKmers(dna);
        return this.baseKmers.distance(other);
    }

    /**
     * @return the maxDist
     */
    public double getMaxDist() {
        return this.maxDist;
    }

}
