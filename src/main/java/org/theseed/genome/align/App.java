package org.theseed.genome.align;

import java.util.Arrays;

import org.theseed.basic.BaseProcessor;

/**
 * Commands for Alignment-related utilities.
 *
 * gtos			align a list of GTOs and output the snips
 * genomes		align the genomes in a directory and output the snips
 * splice		align a genome with a reference genome to fill in the gaps
 * diff			compute the number of proteins in each of a set of genomes that differ from a base strain
 * snipCount	count the snip changes in a groups.snips.tbl file
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        // Determine the command to process.
        switch (command) {
        case "gtos" :
            processor = new GtoAlignProcessor();
            break;
        case "genomes" :
            processor = new GenomeAlignProcessor();
            break;
        case "splice" :
            processor = new SpliceProcessor();
            break;
        case "diff" :
            processor = new DiffProcessor();
            break;
        case "snipCount" :
            processor = new SnipCountProcessor();
            break;
        default:
            throw new RuntimeException("Invalid command " + command);
        }
        // Process it.
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
