#!/bin/bash
# lib_mail.sh
#
# Shared mail helper for the temporary QC-DB Slurm job chain (tmp-*.sbatch).
# Reuses the pattern from TENT/TENT2/Scripts/03-BLASTp-search.sh.

MAIL_TO="${PIPELINE_MAIL_TO:-siwawa@snu.ac.kr}"
MAIL_SCRIPT="/rna/liha/tools/send_mail/mail.py"

send_mail() {
    local subject="$1" body="$2"
    [[ -n "$MAIL_TO" && -s "$MAIL_SCRIPT" ]] || return 0

    local mail_status
    set +e
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate mail
    python3 "$MAIL_SCRIPT" --to "$MAIL_TO" --subject "$subject" --body "$body"
    mail_status=$?
    conda deactivate
    set -e

    [[ "$mail_status" -eq 0 ]] || echo "[WARNING] mail send failed (exit ${mail_status}): ${subject}" >&2
}

# fail <message>: log, mail a FAILED notice, and exit 1. Any downstream
# chain step (submitted with --dependency=afterok) will simply never run.
fail() {
    echo "[ERROR] $1" >&2
    send_mail "QC-DB pipeline FAILED" "$1"
    exit 1
}
