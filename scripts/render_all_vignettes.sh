#!/bin/bash
# Render all vignettes in vignettes/ to gfm; log pass/fail and elapsed time.
set -u
cd "$(dirname "$0")/.." || exit 1
LOG=/tmp/vign_run_log.csv
FAIL=/tmp/vign_failures.txt
FAIL_LOGS=/tmp/vign_failure_logs
mkdir -p "$FAIL_LOGS"
echo "vignette,status,elapsed_seconds" > "$LOG"
: > "$FAIL"

for d in $(ls -d vignettes/[0-9]*/ | sort -V); do
    name=$(basename "$d")
    qmd="$d${name}.qmd"
    if [[ ! -f "$qmd" ]]; then
        echo "$name,missing_qmd,0" >> "$LOG"
        echo "$name MISSING $qmd" >> "$FAIL"
        continue
    fi
    echo "[$(date +%H:%M:%S)] rendering $name ..."
    t0=$(date +%s)
    if (cd "$d" && quarto render "${name}.qmd" --to gfm) > "$FAIL_LOGS/${name}.log" 2>&1; then
        t1=$(date +%s)
        echo "$name,ok,$((t1-t0))" >> "$LOG"
        rm -f "$FAIL_LOGS/${name}.log"
    else
        t1=$(date +%s)
        echo "$name,fail,$((t1-t0))" >> "$LOG"
        echo "$name FAILED in $((t1-t0))s" >> "$FAIL"
    fi
done
echo "DONE. log=$LOG fail=$FAIL fail_logs=$FAIL_LOGS"
