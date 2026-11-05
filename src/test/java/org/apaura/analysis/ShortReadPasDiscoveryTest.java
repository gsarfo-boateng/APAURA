package org.apaura.analysis;

import org.apaura.model.Strand;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class ShortReadPasDiscoveryTest {

    @Test
    public void discoversReplicatedPlusStrandBreakpoint() {
        List<ShortReadPasDiscovery.SampleProfile> profiles =
                new ArrayList<ShortReadPasDiscovery.SampleProfile>();

        profiles.add(profile("S1", Strand.PLUS, 60, 40, 4));
        profiles.add(profile("S2", Strand.PLUS, 61, 42, 5));
        profiles.add(profile("S3", Strand.PLUS, 59, 38, 4));
        profiles.add(flatProfile("S4", Strand.PLUS, 40));

        ShortReadPasDiscovery.Result result =
                ShortReadPasDiscovery.discover(
                        profiles,
                        20,
                        10.0,
                        0.50,
                        2,
                        4,
                        20
                );

        assertEquals(1, result.accepted.size());

        ShortReadPasDiscovery.Candidate candidate =
                result.accepted.get(0);

        assertEquals(3, candidate.supportingSamples);
        assertTrue(candidate.position0 >= 158);
        assertTrue(candidate.position0 <= 160);
        assertTrue(candidate.coverageDrop > 0.80);
    }

    @Test
    public void handlesMinusStrandInTranscriptOrientation() {
        List<ShortReadPasDiscovery.SampleProfile> profiles =
                new ArrayList<ShortReadPasDiscovery.SampleProfile>();

        profiles.add(profile("S1", Strand.MINUS, 60, 40, 4));
        profiles.add(profile("S2", Strand.MINUS, 61, 40, 4));
        profiles.add(profile("S3", Strand.MINUS, 59, 40, 5));

        ShortReadPasDiscovery.Result result =
                ShortReadPasDiscovery.discover(
                        profiles,
                        20,
                        10.0,
                        0.50,
                        2,
                        4,
                        20
                );

        assertEquals(1, result.accepted.size());

        ShortReadPasDiscovery.Candidate candidate =
                result.accepted.get(0);

        assertEquals(3, candidate.supportingSamples);
        assertTrue(candidate.position0 >= 139);
        assertTrue(candidate.position0 <= 141);
    }

    @Test
    public void flatCoverageProducesNoCandidate() {
        List<ShortReadPasDiscovery.SampleProfile> profiles =
                new ArrayList<ShortReadPasDiscovery.SampleProfile>();

        profiles.add(flatProfile("S1", Strand.PLUS, 40));
        profiles.add(flatProfile("S2", Strand.PLUS, 42));
        profiles.add(flatProfile("S3", Strand.PLUS, 38));

        ShortReadPasDiscovery.Result result =
                ShortReadPasDiscovery.discover(
                        profiles,
                        20,
                        10.0,
                        0.50,
                        2,
                        4,
                        20
                );

        assertEquals(0, result.accepted.size());
        assertEquals(0, result.candidates.size());
    }

    @Test
    public void singleSampleArtifactIsRejected() {
        List<ShortReadPasDiscovery.SampleProfile> profiles =
                new ArrayList<ShortReadPasDiscovery.SampleProfile>();

        profiles.add(profile("S1", Strand.PLUS, 60, 40, 3));
        profiles.add(flatProfile("S2", Strand.PLUS, 40));
        profiles.add(flatProfile("S3", Strand.PLUS, 40));

        ShortReadPasDiscovery.Result result =
                ShortReadPasDiscovery.discover(
                        profiles,
                        20,
                        10.0,
                        0.50,
                        2,
                        4,
                        20
                );

        assertEquals(0, result.accepted.size());
        assertEquals(1, result.candidates.size());
        assertEquals(
                "insufficient_replicate_support",
                result.candidates.get(0).status
        );
    }

    private static ShortReadPasDiscovery.SampleProfile profile(
            String sample,
            Strand strand,
            int breakpoint,
            int upstreamDepth,
            int downstreamDepth
    ) {
        int[] transcriptOriented = new int[100];

        for (int i = 0; i < transcriptOriented.length; i++) {
            transcriptOriented[i] =
                    i < breakpoint
                            ? upstreamDepth
                            : downstreamDepth;
        }

        int[] genomic = strand == Strand.MINUS
                ? reverse(transcriptOriented)
                : transcriptOriented;

        return new ShortReadPasDiscovery.SampleProfile(
                sample,
                "chr1",
                strand,
                100,
                200,
                genomic
        );
    }

    private static ShortReadPasDiscovery.SampleProfile flatProfile(
            String sample,
            Strand strand,
            int depth
    ) {
        int[] values = new int[100];

        for (int i = 0; i < values.length; i++) {
            values[i] = depth;
        }

        return new ShortReadPasDiscovery.SampleProfile(
                sample,
                "chr1",
                strand,
                100,
                200,
                values
        );
    }

    private static int[] reverse(int[] values) {
        int[] reversed = new int[values.length];

        for (int i = 0; i < values.length; i++) {
            reversed[i] = values[values.length - 1 - i];
        }

        return reversed;
    }
}
