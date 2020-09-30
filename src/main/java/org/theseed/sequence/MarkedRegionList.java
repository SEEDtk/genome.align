/**
 *
 */
package org.theseed.sequence;

import java.util.Collection;

import org.theseed.genome.Genome;

/**
 * This is a subclass of the region list that supports a counter.
 *
 * @author Bruce Parrello
 *
 */
public class MarkedRegionList extends RegionList {

    // FIELDS
    /** serialization version ID */
    private static final long serialVersionUID = -8648986586599616750L;
    /** number of regions that have no significant changes */
    int counter;

    /**
     * Create an empty region list with the default capacity.
     */
    public MarkedRegionList() {
        super();
        this.counter = 0;
    }

    /**
     * Create an empty region list with the specified initial capacity.
     *
     * @param initialCapacity	starting capacity of the array list
     */
    public MarkedRegionList(int initialCapacity) {
        super(initialCapacity);
        this.counter = 0;
    }

    /**
     * Create a region list from the specified list collection.
     *
     * @param c		collection of regions to use for initializing the list
     */
    public MarkedRegionList(Collection<? extends ExtendedProteinRegion> c) {
        super(c);
        this.counter = 0;
    }

    /**
     * Create a list of extended regions for a genome. Each extended region will contain a protein and its upstream
     * region on the strand, restrained by a maximum.
     *
     * @param genome			genome of interest
     * @param upstreamLimit		maximum upstream value
     */
    public MarkedRegionList(Genome genome, int upstreamLimit) {
        super(genome, upstreamLimit);
        this.counter = 0;
    }

    /**
     * Increment the counter.
     */
    public void increment() {
        this.counter++;
    }

    /**
     * @return the counter value
     */
    public int getCounter() {
        return this.counter;
    }

}
