# Bug 003: `-resume` ignored changed parameters

| | |
|---|---|
| **Status** | Fixed |
| **Scripts** | `sedna.nf` and the six process files `modules/sedna-module-{1..6}.nf` |
| **Lines modified** | `sedna.nf` 88–110 (new `moduleN_opts` maps), 138, 142, 146, 165, 191 and 199–200 (opts passed to MODULE1–6). `modules/sedna-module-1.nf` 18–21, 31–34. `modules/sedna-module-2.nf` 18–21, 31–34. `modules/sedna-module-3.nf` 18–21, 31–32. `modules/sedna-module-4.nf` 24–27, 38. `modules/sedna-module-5.nf` 20–23, 32–35 (and the `def cut_ga` line after 28 removed). `modules/sedna-module-6.nf` 21–24, 33–35. In each process file: the new `val(opts)` input, and `${params.x}` replaced by `${opts.x}` in the script |
| **Fixed in commit** | `174e6f2` (174e6f2ec5f5663e398eedec6dad918446bfc821), 2026-10-14 |
| **Present since** | the Nextflow port, `ae09487`; still in `731dcce` |
| **Affects** | every `-resume` run after a parameter change (Nextflow 26.04) |

## Problem

Each process read its settings as `params.*` inside its `script:` block. Examples:
`--evalue_thres`, `--overlap_thres`, `--birch_threshold`, the deepBGC/GECCO thresholds, and the
database paths.

Nextflow decides whether a task can be reused from its task hash. A `-dump-hashes` run shows that,
in Nextflow 26.04, the hash covers:
- the process script text, with `${params.x}` not yet filled in;
- the inputs;
- the container;
- the `bin/` script.

The param values themselves are not included. So after a parameter change, `-resume` found every
task cached and silently returned results computed with the old value.

Example: a run with `--evalue_thres 0.01` and `-resume` reused all the module 5 tasks of the
previous run, which had used 0.001.

## Fix

`sedna.nf` builds one map per module (`module1_opts` … `module6_opts`) holding the params that
module's script uses. Each process receives it as a `val(opts)` input and reads `${opts.x}`.
Inputs are part of the hash, so a changed value now re-runs that module and everything
downstream.

Not changed:
- `publishDir` and `containerOptions` still read `params`, because they do not change results.
- `--derep_version` was already passed as an input path (`731dcce`).

## Verification

Six resumed runs on the three test MAGs, each changing one parameter:

| Run | Re-ran | Cached |
|---|---|---|
| baseline | all 16 tasks | none |
| identical re-run | none | all 16 |
| `--evalue_thres 0.01` | MODULE5 ×3, MODULE6 | MODULE1–4 |
| `--overlap_thres 0.5` | MODULE4–6 | MODULE1–3 |
| `--birch_threshold 0.3` | MODULE6 | MODULE1–5 |
| `--gecco_score_thres 0.9` | MODULE3–6 | MODULE1–2 |

With default parameters, the clustering is identical to the run before the fix.
