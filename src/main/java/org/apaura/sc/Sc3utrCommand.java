package org.apaura.sc;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Adapter for the preserved sc3UTR Java 8 core.
 *
 * The embedded sc3UTR JAR is checksum-verified before execution.
 * Arguments are forwarded unchanged to preserve sc3UTR behavior.
 */
public final class Sc3utrCommand {

    public static final String VERSION = "0.9.5-java8";

    private static final String CORE_RESOURCE =
            "/org/apaura/sc/sc3utr-core-0.9.4-java8.jar";

    private static final String CORE_FILE =
            "sc3utr-core-0.9.4-java8.jar";

    private static final String CORE_SHA256 =
            "dff45872e3bb155bf600d826514c9dd3f657a1aac73870d8c8f950c7246ead56";

    private Sc3utrCommand() {
    }

    public static void runCore(String[] args) throws Exception {

        MaterializedCore core = materializeCoreJar();

        try {
            List<String> command = new ArrayList<String>();

            command.add(javaExecutable());
            command.add("-jar");
            command.add(core.jar.toString());
            command.addAll(Arrays.asList(args));

            Process process =
                    new ProcessBuilder(command)
                            .inheritIO()
                            .start();

            int status;

            try {
                status = process.waitFor();
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw e;
            }

            if (status != 0) {
                throw new IllegalStateException(
                        "sc3UTR exited with status " + status
                );
            }

        } finally {
            core.close();
        }
    }

    public static String verifyCore() throws Exception {

        MaterializedCore core = materializeCoreJar();

        try {
            return sha256(core.jar);
        } finally {
            core.close();
        }
    }

    private static MaterializedCore materializeCoreJar()
            throws IOException, NoSuchAlgorithmException {

        InputStream input =
                Sc3utrCommand.class.getResourceAsStream(CORE_RESOURCE);

        if (input == null) {
            throw new IOException(
                    "Embedded sc3UTR core resource is missing: "
                            + CORE_RESOURCE
            );
        }

        Path directory =
                Files.createTempDirectory("apaura-sc3utr-");

        Path jar =
                directory.resolve(CORE_FILE);

        try {
            Files.copy(
                    input,
                    jar,
                    StandardCopyOption.REPLACE_EXISTING
            );
        } finally {
            input.close();
        }

        String actual = sha256(jar);

        if (!CORE_SHA256.equals(actual)) {

            deleteQuietly(jar);
            deleteQuietly(directory);

            throw new IOException(
                    "Embedded sc3UTR core checksum mismatch.\n"
                            + "Expected: " + CORE_SHA256 + "\n"
                            + "Observed: " + actual
            );
        }

        return new MaterializedCore(
                directory,
                jar
        );
    }

    private static String sha256(Path path)
            throws IOException, NoSuchAlgorithmException {

        MessageDigest digest =
                MessageDigest.getInstance("SHA-256");

        InputStream input =
                Files.newInputStream(path);

        try {
            byte[] buffer = new byte[8192];
            int n;

            while ((n = input.read(buffer)) != -1) {
                digest.update(buffer, 0, n);
            }

        } finally {
            input.close();
        }

        byte[] hash = digest.digest();

        StringBuilder out =
                new StringBuilder(hash.length * 2);

        for (byte b : hash) {
            out.append(
                    String.format(
                            "%02x",
                            b & 0xff
                    )
            );
        }

        return out.toString();
    }

    private static String javaExecutable() {

        File bin =
                new File(
                        System.getProperty("java.home"),
                        "bin"
                );

        File java =
                new File(
                        bin,
                        "java"
                );

        if (java.isFile()) {
            return java.getAbsolutePath();
        }

        File javaExe =
                new File(
                        bin,
                        "java.exe"
                );

        if (javaExe.isFile()) {
            return javaExe.getAbsolutePath();
        }

        return "java";
    }

    private static void deleteQuietly(Path path) {

        try {
            Files.deleteIfExists(path);

        } catch (IOException ignored) {
            path.toFile().deleteOnExit();
        }
    }

    private static final class MaterializedCore {

        private final Path directory;
        private final Path jar;

        private MaterializedCore(
                Path directory,
                Path jar
        ) {
            this.directory = directory;
            this.jar = jar;
        }

        private void close() {
            deleteQuietly(jar);
            deleteQuietly(directory);
        }
    }
}