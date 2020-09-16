/**
 *
 */
package org.theseed.genome.align;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import junit.framework.TestCase;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

/**
 * @author Bruce Parrello
 *
 */
public class SimpleTest extends TestCase {

    public void testJoin() {
        List<String> cols = Arrays.asList("col1", "col2", "col3");
        String joined = cols.stream().collect(Collectors.joining("\t", "function\t", ""));
        assertThat(joined, equalTo("function\tcol1\tcol2\tcol3"));
    }
}
