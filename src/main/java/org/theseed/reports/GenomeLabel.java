/**
 *
 */
package org.theseed.reports;

import org.theseed.genome.Genome;

/**
 * This object describes the labeling for a genome.  It contains the genome ID, its preferred tooltip, and
 * its preferred display label.  The default display label is the genome ID, and the default tooltip is the
 * genome name.  These can be overridden, however, by creating a custom genome label.  Equality for a genome
 * label is based solely on the genome ID.
 *
 * @author Bruce Parrello
 *
 */
public class GenomeLabel {

    // FIELDS
    /** ID of the relevant genome */
    private String id;
    /** preferred tooltip */
    private String tooltip;
    /** preferred column header */
    private String header;

    /**
     * Construct a default genome label from a genome.
     *
     * @param genome	genome from which the label information should be built
     */
    public GenomeLabel(Genome genome) {
        this.id = genome.getId();
        this.tooltip = genome.getName();
        this.header = genome.getId();
    }

    /**
     * Construct a genome label from a genome ID, label, and tooltip.
     *
     * @param genomeId	ID of the genome of interest
     * @param tooltip	tooltip to use
     * @param header	header label to use
     */
    public GenomeLabel(String genomeId, String tooltip, String header) {
        this.id = genomeId;
        this.tooltip = tooltip;
        this.header = header;
    }

    /**
     * @return the ID of this genome
     */
    public String getId() {
        return this.id;
    }

    /**
     * @return the preferred column header for this genome
     */
    public String getHeader() {
        return this.header;
    }

    /**
     * @return the preferred tooltip for this genome
     */
    public String getToolTip() {
        return this.tooltip;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((this.id == null) ? 0 : this.id.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof GenomeLabel)) {
            return false;
        }
        GenomeLabel other = (GenomeLabel) obj;
        if (this.id == null) {
            if (other.id != null) {
                return false;
            }
        } else if (!this.id.equals(other.id)) {
            return false;
        }
        return true;
    }

}
