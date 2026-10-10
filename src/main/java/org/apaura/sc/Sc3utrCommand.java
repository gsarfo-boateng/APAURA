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
 * Internal adapter for APAURA single-cell analysis.
 *
 * The preserved single-cell core and architecture engines are embedded
 * as resources and checksum-verified before execution.
 */
public final class Sc3utrCommand {

    private static final String CORE_RESOURCE =
            "/org/apaura/sc/sc3utr-core-0.9.4-java8.jar";

    private static final String CORE_FILE =
            "sc3utr-core-0.9.4-java8.jar";

    private static final String CORE_SHA256 =
            "dff45872e3bb155bf600d826514c9dd3f657a1aac73870d8c8f950c7246ead56";

    private static final String ARCH_RESOURCE =
            "/org/apaura/sc/sc3utr-architecture-0.9.5-java8.jar";

    private static final String ARCH_FILE =
            "sc3utr-architecture-0.9.5-java8.jar";

    private static final String ARCH_SHA256 =
            "cda3d113e8d9a59e264daf3392fc065dff2156bda9a4ae3f91957010996c2674";

    private Sc3utrCommand() {
    }

    public static void runCore(String[] args) throws Exception {
        MaterializedJar core = materializeJar(
                CORE_RESOURCE,
                CORE_FILE,
                CORE_SHA256
        );

        try {
            List<String> command = new ArrayList<String>();
            command.add(javaExecutable());
            command.add("-jar");
            command.add(core.jar.toString());
            command.addAll(Arrays.asList(args));

            runProcess(command, "Single-cell core");
        } finally {
            core.close();
        }
    }

    public static void runGtfSites(String[] args) throws Exception {
        runCore(args);

        String annotation = getArg(args, "--annotation");
        String sites = getArg(args, "--out");

        if (annotation == null || sites == null) {
            return;
        }

        File sitesFile = new File(sites);

        if (!sitesFile.isFile()) {
            return;
        }

        String clusterBp = getArg(args, "--cluster-bp");

        if (clusterBp == null) {
            clusterBp = "24";
        }

        String prefix = sites.endsWith(".tsv")
                ? sites.substring(0, sites.length() - 4)
                : sites;

        runArchitecture(new String[]{
                "build",
                "--annotation", annotation,
                "--sites", sites,
                "--out-prefix", prefix,
                "--cluster-bp", clusterBp,
                "--overwrite"
        });
    }

    public static void runReference(String[] args) throws Exception {
        runCore(args);
        annotateResults(args);
    }

    public static void runCompare(String[] args) throws Exception {
        runCore(args);
        annotateResults(args);
        ScApaUsageHarmonizer.run(args);
    }

    private static void annotateResults(String[] args) throws Exception {

        String sites = getArg(args, "--sites");
        String out = getArg(args, "--out");

        if (sites == null || out == null) {
            return;
        }

        String prefix = sites.endsWith(".tsv")
                ? sites.substring(0, sites.length() - 4)
                : sites;

        File geneArchitecture =
                new File(prefix + ".gene_architecture.tsv");

        File outDir = new File(out);

        if (!geneArchitecture.isFile() || !outDir.isDirectory()) {
            return;
        }

        runArchitecture(new String[]{
                "annotate-results",
                "--gene-architecture", geneArchitecture.getPath(),
                "--results-root", outDir.getPath(),
                "--overwrite"
        });
    }

    public static String verifyCore() throws Exception {
        MaterializedJar core = materializeJar(
                CORE_RESOURCE,
                CORE_FILE,
                CORE_SHA256
        );

        try {
            return sha256(core.jar);
        } finally {
            core.close();
        }
    }

    public static String verifyArchitecture() throws Exception {
        MaterializedJar architecture = materializeJar(
                ARCH_RESOURCE,
                ARCH_FILE,
                ARCH_SHA256
        );

        try {
            return sha256(architecture.jar);
        } finally {
            architecture.close();
        }
    }

    private static void runArchitecture(String[] args) throws Exception {
        MaterializedJar architecture = materializeJar(
                ARCH_RESOURCE,
                ARCH_FILE,
                ARCH_SHA256
        );

        try {
            List<String> command = new ArrayList<String>();
            command.add(javaExecutable());
            command.add("-cp");
            command.add(architecture.jar.toString());
            command.add("io.sc3utr.architecture.ArchitectureMain");
            command.addAll(Arrays.asList(args));

            runProcess(command, "APA architecture annotation");
        } finally {
            architecture.close();
        }
    }

    private static void runProcess(List<String> command, String name)
            throws IOException, InterruptedException {

        Process process = new ProcessBuilder(command)
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
                    name + " exited with status " + status
            );
        }
    }

    private static String getArg(String[] args, String key) {
        for (int i = 0; i < args.length - 1; i++) {
            if (key.equals(args[i])) {
                return args[i + 1];
            }
        }

        return null;
    }

    private static MaterializedJar materializeJar(
            String resource,
            String fileName,
            String expectedSha256
    ) throws IOException, NoSuchAlgorithmException {

        InputStream input =
                Sc3utrCommand.class.getResourceAsStream(resource);

        if (input == null) {
            throw new IOException(
                    "Embedded resource is missing: " + resource
            );
        }

        Path directory =
                Files.createTempDirectory("apaura-sc-");

        Path jar =
                directory.resolve(fileName);

        try {
            Files.copy(
                    input,
                    jar,
                    StandardCopyOption.REPLACE_EXISTING
            );
        } finally {
            input.close();
        }

        String actualSha256 = sha256(jar);

        if (!expectedSha256.equals(actualSha256)) {
            deleteQuietly(jar);
            deleteQuietly(directory);

            throw new IOException(
                    "Checksum mismatch for " + fileName
                            + ". Expected " + expectedSha256
                            + " but found " + actualSha256
            );
        }

        return new MaterializedJar(directory, jar);
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

        StringBuilder hash = new StringBuilder();

        for (byte b : digest.digest()) {
            hash.append(String.format("%02x", b & 0xff));
        }

        return hash.toString();
    }

    private static String javaExecutable() {
        File bin = new File(
                System.getProperty("java.home"),
                "bin"
        );

        File java = new File(bin, "java");

        if (java.isFile()) {
            return java.getAbsolutePath();
        }

        File javaExe = new File(bin, "java.exe");

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

    private static final class MaterializedJar {

        private final Path directory;
        private final Path jar;

        private MaterializedJar(Path directory, Path jar) {
            this.directory = directory;
            this.jar = jar;
        }

        private void close() {
            deleteQuietly(jar);
            deleteQuietly(directory);
        }
    }
}