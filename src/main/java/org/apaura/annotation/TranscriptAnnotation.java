package org.apaura.annotation;

import org.apaura.model.Strand;
import org.apaura.model.TranscriptModel;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * Transcript-aware GTF parser used by the short-read annotated-PAS engine. 
 */
public final class TranscriptAnnotation {
    private final List<TranscriptModel> transcripts;
    
    private TranscriptAnnotation(List<TranscriptModel> transcripts) {
        this.transcripts = transcripts;
    }
    
    public List<TranscriptModel> transcripts() {
        return Collections.unmodifiableList(transcripts);
    }
    
    public static TranscriptAnnotation load(String path) throws IOException {
        Map<String, TranscriptModel> tx = new LinkedHashMap<String, TranscriptModel>();
        InputStream in = new FileInputStream(path);
        if (path.endsWith(".gz")) in = new GZIPInputStream(in);
        BufferedReader br = new BufferedReader(new InputStreamReader(in, "UTF-8"));
        String line;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty() || line.charAt(0) == '#') continue;
            String[] f = line.split("\t");
            if (f.length < 9) continue;
            String feature = f[2].toLowerCase(Locale.ROOT);
            if (!feature.equals("exon") && !feature.equals("cds")) continue;
            Map<String, String> a = parseAttributes(f[8]);
            String tid = a.get("transcript_id");
            String gid = a.get("gene_id");
            if (tid == null || gid == null) continue;
            String gname = first(a, "gene_name");
            if (gname == null) gname = gid;
            String bio = first(a, "transcript_type", "transcript_biotype", "gene_type", "gene_biotype");
            if (bio == null) bio = ".";
            String key = gid + "\u0001" + tid;
            TranscriptModel t = tx.get(key);
            if (t == null) {
                t = new TranscriptModel(tid, gid, gname, bio, f[0], Strand.parse(f[6]));
                tx.put(key, t);
            }
            int s = Integer.parseInt(f[3]) - 1;
            int e = Integer.parseInt(f[4]);
            if (feature.equals("exon")) t.addExon(s, e); else t.addCds(s, e);
        }
        br.close();
        List<TranscriptModel> out = new ArrayList<TranscriptModel>();
        for (TranscriptModel t : tx.values()) if (!t.exons.isEmpty() && t.strand != Strand.UNKNOWN) out.add(t);
        System.err.println("[APAURA] loaded transcript models=" + out.size());
        return new TranscriptAnnotation(out);
    }
    
    private static String first(Map<String, String> a, String... keys) {
        for (String k : keys) if (a.containsKey(k)) return a.get(k);
        return null;
    }
    
    private static Map<String, String> parseAttributes(String s) {
        Map<String, String> m = new HashMap<String, String>();
        for (String x : s.split(";")) {
            x = x.trim();
            if (x.isEmpty()) continue;
            int sp = x.indexOf(' ');
            int eq = x.indexOf('=');
            int p = (sp >= 0 && (eq < 0 || sp < eq)) ? sp : eq;
            if (p < 0) continue;
            String k = x.substring(0, p).trim();
            String v = x.substring(p + 1).trim();
            if (v.startsWith("\"") && v.endsWith("\"") && v.length() >= 2) v = v.substring(1, v.length() - 1);
            m.put(k, v);
        }
        return m;
    }
}