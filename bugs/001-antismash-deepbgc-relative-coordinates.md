# Bug 001: antiSMASH and deepBGC BGC coordinates were relative to the BGC, not the contig

| | |
|---|---|
| **Status** | Fixed |
| **Script** | `bin/utils.py` (`antismash_annot_parser`, `deepbgc_annot_parser`) |
| **Lines modified** | `bin/utils.py` 311–316 (antiSMASH) and 380–385 (deepBGC) |
| **Fixed in commit** | `fd50a17` (fd50a17e08c9e6259e0f9e67a192b949423f15b2), 2026-10-14 |
| **Present since** | the legacy code: `old/bin/src/bgc_annot/utilities.py` lines 330 and 416; still in `89b50d4` |
| **Affects** | modules 1 and 2 (`annot_metadata.tsv`), module 4 (dereplication), `MODULE5_PRECOMPUTED` (CDS selection by coordinates) |

## Problem

Both parsers took a BGC's `start` and `end` from a feature inside its GBK (`cand_cluster` for
antiSMASH, `cluster` for deepBGC). Those GBKs hold only the BGC's sequence, so the coordinates
were relative to the BGC and `start` was always about 0. GECCO's parser reads its
`*.clusters.tsv`, which is in contig coordinates. Module 4 therefore compared antiSMASH and
deepBGC positions with GECCO positions measured from a different origin.

Example (`OceanDNA-b45222_redu`): deepBGC `OceanDNA-b45222_00387` lies at 2–3137, exactly like
GECCO's prediction, but was recorded as 0–3135, so module 4 treated the two as staggered
instead of identical. deepBGC `OceanDNA-b45222_00076` (257–5379) was recorded as 0–5122.

## Fix

The parsers now add the BGC's start on the contig (zero based) to the feature coordinates:

- **antiSMASH:** `Orig. start` from the GBK's `antiSMASH-Data` structured comment.
- **deepBGC:** the `_<start>-<end>` suffix of the record's accession.

## Verification

- The deepBGC rows match deepBGC's own `.bgc.tsv` (`nucl_start`/`nucl_end`) for all three test MAGs.
- The antiSMASH parser returns 1000–30137 for a test GBK whose `Orig. start` was set to 1000.
- After a full pipeline run, deepBGC and GECCO agree exactly on `OceanDNA-b45222_00387`
  (2–3137).
- Changed `bgc_id`s on the test data: `b45222_00076` (0–5122 → 257–5379), `b45222_00083`
  (0–5266 → 1–5267), `b45222_00387` (0–3135 → 2–3137) and `b9899_00139` (0–3422 → 12–3434).
  The set of dereplicated BGCs and the GCF assignment did not change.
- All antiSMASH test regions start at position 0, so their outputs are unchanged.
