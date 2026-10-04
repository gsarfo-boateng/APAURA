package org.apaura.analysis;

import org.apaura.model.EndKey;
import org.apaura.model.Sample;
import org.apaura.model.SiteCluster;
import org.apaura.model.Strand;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;

class SiteClustererTest {

    @Test
    void modeUsesMostSupportedCoordinateAndKeepsSampleCounts() {
        List<Sample> samples = new ArrayList<Sample>();
        samples.add(new Sample("C1", "control", "C1.bam"));
        samples.add(new Sample("T1", "treated", "T1.bam"));

        Map<String, Map<EndKey, Integer>> counts =
                new HashMap<String, Map<EndKey, Integer>>();
        counts.put("C1", sampleCounts(100, 2, 101, 5));
        counts.put("T1", sampleCounts(100, 1, 101, 3));

        List<SiteCluster> clusters = SiteClusterer.cluster(samples, counts, 1, "mode");

        assertEquals(1, clusters.size());
        assertEquals(101, clusters.get(0).representativePos0);
        assertEquals(7, clusters.get(0).counts.get("C1"));
        assertEquals(4, clusters.get(0).counts.get("T1"));
    }

    private static Map<EndKey, Integer> sampleCounts(
            int firstPosition,
            int firstCount,
            int secondPosition,
            int secondCount
    ) {
        Map<EndKey, Integer> counts = new HashMap<EndKey, Integer>();
        counts.put(new EndKey("chr1", firstPosition, Strand.PLUS), firstCount);
        counts.put(new EndKey("chr1", secondPosition, Strand.PLUS), secondCount);
        return counts;
    }
}
