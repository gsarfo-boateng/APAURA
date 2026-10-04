# APA architecture model

APAURA separates **architecture** from **direction**.

Architecture asks what terminal structures exist. Direction asks how their usage changes.

## TUTR

Multiple PAS within the same terminal-exon branch.

Typical direction labels include `shortening`, `lengthening`, and `no_directional_shift`.

## ALE

Alternative terminal-exon branches.

Direction can be summarized as proximal/distal switching where supported.

## mixed

Both TUTR and ALE architecture are present.

## none

No terminal architecture class is assigned under the current annotation/quantification rules.

## IPA

Intronic polyadenylation is tracked separately rather than forced into a terminal architecture label.

## Transcript orientation

Proximal/distal is defined in transcript 5′->3′ orientation.

Plus strand:

```text
lower genomic coordinate -> more proximal
higher genomic coordinate -> more distal
```

Minus strand:

```text
higher genomic coordinate -> more proximal
lower genomic coordinate -> more distal
```
