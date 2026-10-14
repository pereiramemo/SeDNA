# Bug 002: dereplication re-added BGCs that had already lost a comparison

| | |
|---|---|
| **Status** | Fixed |
| **Script** | `bin/utils.py` (`dereplicate_shared_contigs`, called by `recursive_dereplication` in `bin/sedna-module-4.py`) |
| **Lines modified** | `bin/utils.py` 706–710 (new `lost_bgcs` set); one `lost_bgcs.add(...)` line per overlap case at 757, 771, 785, 799, 827, 861, 875, 889, 903, 932, 966; 971–972 (case 13 deletes `id2`); 974–975 (losers removed after all pairs) |
| **Fixed in commit** | `9fd90a8` (9fd90a8baf58b315c54f69b1db34a5d81180584c), 2026-10-14 |
| **Present since** | the legacy code: `old/bin/src/bgc_derep/dsc.py`; still in `8942cf4` |
| **Affects** | module 4: the dereplicated catalog (and so modules 5 and 6) |

## Problem

`dereplicate_shared_contigs` compares every BGC of one table with every BGC of the other on each
shared contig, one pair at a time. When two BGCs overlap, the shorter one is removed from the
catalog. When two BGCs are disjoint (cases 6 and 12), **both** are added. A BGC that had already
lost to a longer one in an earlier pair was therefore put back whenever a later pair found it
disjoint from some other BGC. The result depended on the order of the pairs.

Any contig with two or more predictions from the same tool was exposed, even with only two tools.
Example: toolA has BGCs at 1000–9000 and 12500–15000, and toolB has BGCs at 2000–8000 and
22000–25000, all on one contig.
1. toolA 1000–9000 vs toolB 2000–8000 is case 8 (nested), so toolB 2000–8000 is removed.
2. toolA 12500–15000 vs toolB 2000–8000 is case 6 (disjoint), so toolB 2000–8000 is added back.

The catalog ended with both toolA 1000–9000 and the toolB 2000–8000 nested inside it.

Case 13 (identical coordinates) also deleted `dereplicated_bgcs[tool2][id1]` instead of `[id2]`
(docs §11.7).

## Fix

Each overlap case now records its losing BGC in `lost_bgcs`. Once all pairs of the contig round
have been compared, every loser is removed from the catalog.

- **Tie on length** (cases 4, 5, 10, 11, 13): tool 1's BGC is kept, as before.
- **Unchanged:** the 13 cases, longest-wins and `--overlap_thres`, as well as the audit tables
  (`overlapped/`, `partially_overlapped/`, `non_overlapped/`).
- **Case 13** now deletes `id2`.

## Verification

Synthetic fixtures with fake tools on the test contigs (`tests/output/ntool/`, git-ignored):

- **Two tools** (the example above): the catalog keeps toolA 1000–9000, toolA 12500–15000 and
  toolB 22000–25000. toolB 2000–8000 is no longer kept.
- **Five tools:** an all-overlap group of four tools now gives one BGC (toolA 1000–9000) instead of
  five, and a chain of four tools gives one (toolC 12000–26000).
- **The three test MAGs:** no contig has two predictions from the same tool, so their results do
  not change.
