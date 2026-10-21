package org.apaura.stats;

import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertTrue;

class ModeratedDirichletMultinomialTest {

    @Test
    void ranksStrongCompositionShiftAheadOfStableEvent() {
        List<ModeratedDirichletMultinomial.Event> events =
                new ArrayList<ModeratedDirichletMultinomial.Event>();

        events.add(new ModeratedDirichletMultinomial.Event(
                "stable",
                new int[][]{
                        {80, 10, 10},
                        {78, 12, 10},
                        {82, 8, 10},
                        {80, 10, 10},
                        {78, 12, 10},
                        {82, 8, 10}
                },
                new int[]{0, 0, 0, 1, 1, 1}
        ));

        events.add(new ModeratedDirichletMultinomial.Event(
                "shifted",
                new int[][]{
                        {80, 10, 10},
                        {78, 12, 10},
                        {82, 8, 10},
                        {10, 10, 80},
                        {12, 8, 80},
                        {8, 12, 80}
                },
                new int[]{0, 0, 0, 1, 1, 1}
        ));

        ModeratedDirichletMultinomial.BatchResult batch =
                ModeratedDirichletMultinomial.test(events);

        ModeratedDirichletMultinomial.Result stable =
                batch.get("stable");
        ModeratedDirichletMultinomial.Result shifted =
                batch.get("shifted");

        assertTrue(Double.isFinite(batch.commonPrecision));
        assertTrue(batch.commonPrecision > 0.0);
        assertTrue(stable != null && shifted != null);
        assertTrue(Double.isFinite(stable.pvalue));
        assertTrue(Double.isFinite(shifted.pvalue));
        assertTrue(shifted.pvalue < stable.pvalue);
        assertTrue(shifted.pvalue < 0.05);
    }

    @Test
    void supportsTwoFeatureIpaLikeEvents() {
        List<ModeratedDirichletMultinomial.Event> events =
                new ArrayList<ModeratedDirichletMultinomial.Event>();

        events.add(new ModeratedDirichletMultinomial.Event(
                "ipa",
                new int[][]{
                        {15, 85},
                        {18, 82},
                        {12, 88},
                        {80, 20},
                        {82, 18},
                        {78, 22}
                },
                new int[]{0, 0, 0, 1, 1, 1}
        ));

        ModeratedDirichletMultinomial.Result result =
                ModeratedDirichletMultinomial.test(events).get("ipa");

        assertTrue(result != null);
        assertTrue(result.converged);
        assertTrue(Double.isFinite(result.pvalue));
        assertTrue(result.pvalue < 0.05);
    }
}
