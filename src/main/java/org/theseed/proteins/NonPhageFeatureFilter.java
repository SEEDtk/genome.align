/**
 *
 */
package org.theseed.proteins;

import org.theseed.genome.Feature;

/**
 * This feature filters out phage-related roles.
 *
 * @author Bruce Parrello
 */
public class NonPhageFeatureFilter extends FeatureFilter {

    /**
     * Construct a non-phage feature filter.
     *
     * @param processor		constructing client
     */
    public NonPhageFeatureFilter(IParms processor) {
        super(processor);
    }

    @Override
    public boolean filter(Feature feat) {
        String function = feat.getPegFunction();
        return ! function.contains("phage");
    }

}
