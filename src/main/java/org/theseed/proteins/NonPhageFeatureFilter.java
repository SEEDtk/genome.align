/**
 *
 */
package org.theseed.proteins;

import java.io.IOException;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Feature;
import org.theseed.utils.ParseFailureException;

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
     *
     * @throws IOException
     * @throws ParseFailureException
     */
    public NonPhageFeatureFilter(IParms processor) throws ParseFailureException, IOException {
        super(processor);
    }

    @Override
    public boolean filter(Feature feat) {
        String function = feat.getPegFunction();
        return ! StringUtils.containsIgnoreCase(function, "phage");
    }

}
