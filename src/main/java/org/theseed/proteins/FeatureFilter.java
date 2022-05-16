/**
 *
 */
package org.theseed.proteins;

import java.io.IOException;
import java.util.Set;

import org.theseed.genome.Feature;
import org.theseed.utils.ParseFailureException;

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

        /**
         * @return the set of acceptable feature IDs
         *
         * @throws IOException
         * @throws ParseFailureException
         */
        public Set<String> getFeatureSet() throws ParseFailureException, IOException;

    }

    /**
     * Construct a filter.
     *
     * @param processor		constructing client, containing parameters
     *
     * @throws IOException
     * @throws ParseFailureException
     */
    public FeatureFilter(IParms processor) throws ParseFailureException, IOException {
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
        NONPHAGE, LIST;

        /**
         * @return a feature filter of the specified type.
         *
         * @param processor		constructing client
         *
         * @throws IOException
         * @throws ParseFailureException
         */
        public FeatureFilter create(IParms processor) throws ParseFailureException, IOException {
            FeatureFilter retVal = null;
            switch (this) {
            case NONPHAGE:
                retVal = new NonPhageFeatureFilter(processor);
                break;
            case LIST:
                retVal = new ListFeatureFilter(processor);
            }
            return retVal;
        }

    }

}
