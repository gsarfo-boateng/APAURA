package org.apaura.io;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class SamTextReaderTest {

    @Test
    void cigarParserDoesNotCountSkippedIntronAsAlignedBlock() {
        SamTextReader.CigarInfo cigar = SamTextReader.parseCigar("50M1000N50M", 100);

        assertEquals(1100, cigar.referenceConsumed);
        assertEquals(100, cigar.queryConsumed);
        assertEquals(2, cigar.blocks.size());
        assertEquals(100, cigar.blocks.get(0).start0);
        assertEquals(150, cigar.blocks.get(0).end0);
        assertEquals(1150, cigar.blocks.get(1).start0);
        assertEquals(1200, cigar.blocks.get(1).end0);
    }
}
