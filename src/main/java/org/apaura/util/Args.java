package org.apaura.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public final class Args {
    private final Map<String, String> values = new LinkedHashMap<String, String>();
    private final Set<String> flags = new HashSet<String>();
    private final List<String> positional = new ArrayList<String>();

    public Args(String[] args, int start) {
        for (int i = start; i < args.length; i++) {
            String a = args[i];
            if (a.startsWith("--")) {
                String key = a.substring(2);
                int eq = key.indexOf('=');
                if (eq >= 0) {
                    values.put(key.substring(0, eq), key.substring(eq + 1));
                } else if (i + 1 < args.length && !args[i + 1].startsWith("--")) {
                    values.put(key, args[++i]);
                } else {
                    flags.add(key);
                }
            } else {
                positional.add(a);
            }
        }
    }

    public boolean has(String key) {
        return values.containsKey(key) || flags.contains(key);
    }

    public String get(String key, String def) {
        return values.containsKey(key) ? values.get(key) : def;
    }

    public String require(String key) {
        String v = values.get(key);
        if (v == null || v.isEmpty()) throw new IllegalArgumentException("Missing required option --" + key);
        return v;
    }

    public int getInt(String key, int def) {
        return Integer.parseInt(get(key, Integer.toString(def)));
    }

    public double getDouble(String key, double def) {
        return Double.parseDouble(get(key, Double.toString(def)));
    }

    public boolean flag(String key) {
        return flags.contains(key) || "true".equalsIgnoreCase(values.get(key));
    }

    public List<String> positional() {
        return Collections.unmodifiableList(positional);
    }
}