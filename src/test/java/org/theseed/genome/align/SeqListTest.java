package org.theseed.genome.align;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;

/**
 * Unit test for simple App.
 */
public class SeqListTest
    extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public SeqListTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( SeqListTest.class );
    }

    /**
     * Test sequence lists
     * @throws IOException
     */
    public void testSeqList() throws IOException {
        Genome gto = new Genome(new File("data", "W3110-wild.gto"));
        Feature feat = gto.getFeature("fig|316407.41.peg.703");
        Location loc = feat.getLocation();
        String baseDna = gto.getDna(loc);
        DnaKmers baseKmers = new DnaKmers(baseDna);
        SequenceList list = new SequenceList(feat.getId(), loc.toString(), baseDna);
        assertThat(list.size(), equalTo(1));
        assertThat(list.getMaxDist(), equalTo(0.0));
        gto = new Genome(new File("data", "W3110-30318.gto"));
        feat = gto.getFeature("fig|316407.119.peg.3360");
        loc = feat.getLocation();
        String dna1 = gto.getDna(loc);
        double dist = list.distance(dna1);
        DnaKmers kmers = new DnaKmers(dna1);
        assertThat(dist, equalTo(baseKmers.distance(kmers)));
        boolean added = list.add(feat.getId(), loc.toString(), dna1, 0.1);
        assertThat(added, equalTo(true));
        assertThat(list.getMaxDist(), equalTo(dist));
        gto = new Genome(new File("data", "W3110-30317.gto"));
        feat = gto.getFeature("fig|316407.118.peg.3302");
        loc = feat.getLocation();
        String dna2 = gto.getDna(loc);
        double dist2 = list.distance(dna2);
        kmers = new DnaKmers(dna2);
        assertThat(dist2, equalTo(baseKmers.distance(kmers)));
        added = list.add(feat.getId(), loc.toString(), dna2, 0.1);
        assertThat(added, equalTo(true));
        assertThat(list.getMaxDist(), equalTo(dist2));
        gto = new Genome(new File("data", "W3110-30316.gto"));
        feat = gto.getFeature("fig|316407.117.peg.3328");
        loc = feat.getLocation();
        String dna3 = gto.getDna(loc);
        dist = list.distance(dna3);
        kmers = new DnaKmers(dna3);
        assertThat(dist, equalTo(baseKmers.distance(kmers)));
        added = list.add(feat.getId(), loc.toString(), dna3, 0.1);
        assertThat(added, equalTo(false));
        assertThat(list.size(), equalTo(3));
        assertThat(list.getMaxDist(), equalTo(dist2));
        File temp = new File("data", "temp.ser");
        list.save(temp);
        try (FastaInputStream inStream = new FastaInputStream(temp)) {
            Iterator<Sequence> iter = inStream.iterator();
            Sequence seq = iter.next();
            assertThat(seq.getLabel(), equalTo("fig|316407.41.peg.703"));
            assertThat(seq.getSequence(), equalTo(baseDna));
            assertThat(iter.next().getSequence(), equalTo(dna1));
            assertThat(iter.next().getSequence(), equalTo(dna2));
            assertThat(iter.hasNext(), equalTo(false));
        }
    }
}
