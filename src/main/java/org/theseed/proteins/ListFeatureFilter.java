/**
 *
 */
package org.theseed.proteins;

import java.io.IOException;
import java.util.Set;

import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;

/**
 * This is a feature filter that only accepts features from a specified list of feature IDs.  In the
 * context of the alignment process, the feature IDs must be from the base genome.
 *
 * @author Bruce Parrello
 *
 */
public class ListFeatureFilter extends FeatureFilter {

    // FIELDS
    /** set of acceptable feature IDs */
    private Set<String> fidSet;

    /**
     * Construct a list-based feature-filter.
     *
     * @param processor		controlling command processor
     *
     * @throws IOException
     * @throws ParseFailureException
     */
    public ListFeatureFilter(IParms processor) throws ParseFailureException, IOException {
        super(processor);
        this.fidSet = processor.getFeatureSet();
    }

    @Override
    public boolean filter(Feature feat) {
        return this.fidSet.contains(feat.getId());
    }

}
