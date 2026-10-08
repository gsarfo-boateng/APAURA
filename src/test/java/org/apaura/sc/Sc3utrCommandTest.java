package org.apaura.sc;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class Sc3utrCommandTest {

    private static final String EXPECTED_SHA256 =
            "dff45872e3bb155bf600d826514c9dd3f657a1aac73870d8c8f950c7246ead56";

    @Test
    void embeddedCoreHasExpectedChecksum() throws Exception {

        assertEquals(
                EXPECTED_SHA256,
                Sc3utrCommand.verifyCore()
        );
    }
}