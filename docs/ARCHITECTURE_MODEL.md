# APA architecture model

APAURA separates **terminal architecture** from **differential usage** and from
**termination mechanism**. Architecture is a structural property of the retained
terminal PAS and their reconciled terminal exon groups; it does not depend on
treatment, statistical significance, or direction of change.

## Public terminal-architecture classes

| Architecture | Definition |
|---|---|
| `SINGLE_PAS` | one resolved terminal exon group and one retained terminal PAS |
| `TUTR` | one resolved terminal exon group and two or more retained terminal PAS |
| `ALE` | two or more resolved terminal exon groups, with one retained PAS in each group |
| `mixed` | two or more resolved terminal exon groups and at least one group with two or more retained PAS |
| `UNRESOLVED` | terminal architecture cannot be assigned completely from the retained evidence |

The decision logic is:

```text
0 terminal groups or 0 retained terminal PAS  -> UNRESOLVED
1 terminal group + 1 retained PAS             -> SINGLE_PAS
1 terminal group + >=2 retained PAS           -> TUTR
>=2 terminal groups + one PAS in each group   -> ALE
>=2 terminal groups + tandem PAS in any group -> mixed
```

If any retained terminal PAS cannot be assigned to a reconciled terminal exon
group, the gene-level terminal architecture is `UNRESOLVED`. Resolved local
structure-specific metrics can still be retained for the portions of the gene
that are structurally interpretable.

## Terminal exon groups

GTF transcript models can contain redundant or slightly shifted terminal exon
annotations. APAURA therefore reconciles compatible terminal annotations before
classifying architecture. The reconciled unit is described publicly as a
**terminal exon group**.

Some output-column names retain historical `branch` terminology
(`terminal_branch_rank`, `n_terminal_branches`, `n_tutr_branches`) for v0.10.1
compatibility. Those columns refer to reconciled terminal exon groups rather than
to an independent biological ontology.

## TUTR

TUTR represents multiple retained PAS within the same resolved terminal exon
group.

Direction labels include:

```text
shortening
lengthening
no_directional_shift
```

A TUTR change is not determined by `delta_wpi` alone. Distributional movement
can be supported by PAS-usage redistribution even when opposing shifts make the
net weighted-position change small.

## ALE

ALE represents alternative resolved terminal exon groups. Where supported,
direction is summarized as proximal or distal switching in transcript
orientation.

## mixed

`mixed` indicates that both ALE structure and within-group tandem PAS are
present in the same gene.

## IPA

Intronic polyadenylation (IPA) is a **termination mechanism**, not a terminal
architecture class. IPA can co-occur with `SINGLE_PAS`, `TUTR`, `ALE`, `mixed`,
or a globally `UNRESOLVED` terminal architecture.

A mechanism is reported as changed only when its effect threshold is met and
significant site/event support belongs to that mechanism. Gene-wide significance
from another mechanism is not sufficient.

## Transcript orientation

Proximal/distal is defined in transcript 5'->3' orientation.

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

## Single-cell naming

The single-cell architecture files retain uppercase `MIXED` for compatibility
with the integrated sc3UTR output. Bulk/direct/short-read APAURA output uses
lowercase `mixed`. The biological meaning is the same.
