# QC-DB conventions

Context for working in this directory. See `project_brief.md` for the
pipeline's purpose and design; this file covers operational conventions
that aren't obvious from the code alone.

## Slurm: no time limits

Don't add `#SBATCH --time=` to any script in `Scripts/`. The `amd`
partition has `MaxTime=UNLIMITED`, so omitting `--time=` entirely is
sufficient — no need for a large placeholder value.

**Why:** a serial 12h-time-limited run of `01-Check_pooled_counts.sh` was
killed by Slurm mid-taxon and lost *all* progress, because it only wrote
output via one atomic `mv` at the very end. The fix was two-part: (1) drop
the time limit, and (2) redesign `01` to be array-batched with one output
file per taxon (see below), so a kill only costs the one in-flight task,
not the whole run.

**Tradeoff to know:** with no time limit, a job stuck on a dead socket or
hung process now runs forever instead of surfacing via `TIMEOUT`. If a
stage seems stuck, check `squeue`'s elapsed time by hand rather than
assuming Slurm will flag it. This is part of why array-batching (many
short-lived tasks) is preferred over long serial jobs here — one hung task
out of many is far easier to spot and reap than one hung monolithic job.

All jobs use `--nodelist=pichu` — that's where the conda envs, `/bank`
mounts, and prior downloads (BUSCO lineage, `LUCA.h5`) live. Other nodes were checked:
`guanine`, `lapras`, `caffeine` are DOWN; `dicer` and `polygon2` lack `/bank` mounts; `pikachu` lacks python dependency modules. Thus, `pichu` is the sole functional node for QC-DB.

## Dependency chains: afterok vs afterany

Chain steps that follow a **Slurm array** must use
`--dependency=afterany:<arrayjobid>`, not `afterok`. `afterok` on an array
requires *every* task to exit 0 — with ~70 taxa hitting NCBI eutils, at
least one transient failure is expected, and `afterok` would leave the
whole chain permanently stuck (`DependencyNeverSatisfied`, unrecoverable —
this happened once already). The step consuming the array's output (e.g.
`tmp-collapse-notify.sbatch`) does its own pass/fail judgement instead,
based on how much actual output was produced.

Ordinary single-script steps (not arrays) downstream of that can safely
use `afterok` — their exit code is trustworthy.

## Conda envs

- `busco` — **also** has `agat`, `gffread`, and `seqkit`, not just BUSCO (non-obvious
  from the name).
- `omark` — omamer/omark, plus `pytables` (the `tables` module) for
  reading `.h5` files. **Not** `h5py` — it isn't installed here.
- `NCBI-download` — `datasets`, `blastdbcmd`, `epost`, `efetch`.
- `mail` — used only for sending notifications via `mail.py`.

## Mail notifications

`lib_mail.sh` (sourced by the `tmp-*.sbatch` chain scripts) provides
`send_mail` and `fail`. `fail` logs, sends a "FAILED" mail, and exits 1 —
use it instead of a bare `exit 1` in any chain script so failures surface
by mail, not just by someone noticing `squeue` is empty.

## Cross-taxon comparability: unified lineages/DBs

BUSCO always runs against `eukaryota_odb10` and OMArk always against
`LUCA.h5` — deliberately not taxon-specific, so scores are comparable
across the full taxonomic range rather than optimized per-clade.

`Asgard_Archaea` is excluded from all downstream analysis — skip it in any
new script that iterates over taxa.

## Per-taxon files, not one big buffered write

Two independent examples of the same pattern, both existing because a
single script processing ~70 taxa serially and writing output only once
at the end is fragile (loses everything on any interruption):

- **`01-Check_pooled_counts.sh`**: array-batched, one Slurm array task per
  `species_info.tsv` row (`SLURM_ARRAY_TASK_ID` = row index). Each task
  writes `Results/.pooled_counts/<NNN>.tsv` (zero-padded row index, one
  15-column line, no header). It's idempotent: if that file already
  exists with a non-`NA` result, the task exits immediately. Resubmitting
  the same array script is how `NA` rows (transient IPG failures) get
  retried — there's no separate retry script.
  `01-Collapse-pooled-counts.sh` concatenates all present per-row files
  (in row order) into `species_info_pooled_counts.tsv`, reporting how many
  rows are missing/still-NA without treating that as fatal itself.
- **`03-BUSCO-and-OMArk.sh`**: per-taxon, writes
  `Results/.done/<Taxonomy>_<taxid>.done` as the very last action, only
  reached if every step above it succeeded (`set -euo pipefail`). Lets you
  check progress across all taxa (`ls Results/.done/ | wc -l`) or gate a
  rerun without opening every taxon's own results folder.

Concurrency for the `01` array is capped at `%3`, matching NCBI eutils'
3 req/s per-IP limit (no API key is configured — getting one would allow
raising both the per-request rate and this concurrency cap safely).

## `tmp-*` files

Everything prefixed `tmp-` in `Scripts/` is one-off chain orchestration,
not part of the permanent pipeline — safe to delete once a chain run has
gone through. Don't hardcode Slurm job IDs as defaults inside these
scripts; job IDs are stale within a day and a stale default silently
starts depending on the wrong (or a long-finished) job.
