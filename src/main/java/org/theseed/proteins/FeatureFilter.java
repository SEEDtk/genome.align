/**
 *
 */
package org.theseed.proteins;

import org.theseed.genome.Feature;

/**
 * THis is the base class for feature filters.  These take a feature as input and either pass or fail the feature.
 *
 * @author Bruce Parrello
 */
public abstract class FeatureFilter {

    /**
     * This interface should be supported by the processing clients.  It is used to extract parameters
     * from the constructing client.
     */
    public interface IParms {

    }

    /**
     * Construct a filter.
     *
     * @param processor		constructing client, containing parameters
     */
    public FeatureFilter(IParms processor) {
    }

    /**
     * @return TRUE if the specified feature should be accepted, else FALSE
     *
     * @param feat		feature to check
     */
    public abstract boolean filter(Feature feat);


    /**
     * This enum describes the typer of filters.
     */
    public static enum Type {
        NONPHAGE;

        /**
         * @return a feature filter of the specified type.
         *
         * @param processor		constructoing client
         */
        public FeatureFilter create(IParms processor) {
            FeatureFilter retVal = null;
            switch (this) {
            case NONPHAGE:
                retVal = new NonPhageFeatureFilter(processor);
                break;
            }
            return retVal;
        }

    }

}
