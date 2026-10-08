package org.apaura.sc;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class Sc3utrCommandTest {

    @Test
    void embeddedArchitectureHasExpectedChecksum() throws Exception {

        assertEquals(
                "cda3d113e8d9a59e264daf3392fc065dff2156bda9a4ae3f91957010996c2674",
                Sc3utrCommand.verifyArchitecture()
        );
    }
}