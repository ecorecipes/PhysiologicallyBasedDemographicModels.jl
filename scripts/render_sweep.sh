#!/usr/bin/env bash
# Re-render every numbered vignette to GitHub-flavoured Markdown and record
# pass/fail + elapsed time in a report at vignettes/_reports/render-sweep.md.
#
# Goal: catch drift in any vignette not touched recently. Uses --to gfm for
# speed (no PDF rendering) and does not pre-purge _freeze/ so Quarto can
# reuse cached cells where the source has not changed.
#
# Usage:
#   scripts/render_sweep.sh               # render all
#   scripts/render_sweep.sh 35 40 41      # render a subset by prefix number

set -u

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VIGDIR="$ROOT/vignettes"
REPORT_DIR="$VIGDIR/_reports"
REPORT="$REPORT_DIR/render-sweep.md"
LOGDIR="$REPORT_DIR/render-logs"

mkdir -p "$REPORT_DIR" "$LOGDIR"

if [ "$#" -gt 0 ]; then
    DIRS=()
    for num in "$@"; do
        for d in "$VIGDIR"/"${num}"_*/; do
            [ -d "$d" ] && DIRS+=("$d")
        done
    done
else
    DIRS=("$VIGDIR"/[0-9][0-9]_*/)
fi

total=${#DIRS[@]}
echo "[sweep] rendering $total vignettes"

pass=0
fail=0
results=()
start_all=$(date +%s)

for i in "${!DIRS[@]}"; do
    d="${DIRS[$i]}"
    name="$(basename "$d")"
    qmd="$d$name.qmd"
    if [ ! -f "$qmd" ]; then
        echo "[sweep] $((i+1))/$total $name — SKIP (no qmd)"
        results+=("SKIP|$name|0|no .qmd file")
        continue
    fi
    log="$LOGDIR/$name.log"
    t0=$(date +%s)
    printf "[sweep] %2d/%d %-40s " "$((i+1))" "$total" "$name"
    (cd "$d" && quarto render "$name.qmd" --to gfm --execute) \
        >"$log" 2>&1
    status=$?
    t1=$(date +%s)
    elapsed=$((t1 - t0))
    if [ "$status" -eq 0 ]; then
        pass=$((pass + 1))
        printf "PASS (%ds)\n" "$elapsed"
        results+=("PASS|$name|$elapsed|")
    else
        fail=$((fail + 1))
        err_tail=$(tail -5 "$log" | tr '\n' ' ' | tr '|' '/' | cut -c1-240)
        printf "FAIL (%ds)\n" "$elapsed"
        results+=("FAIL|$name|$elapsed|$err_tail")
    fi
done

total_elapsed=$(( $(date +%s) - start_all ))

{
    echo "# Vignette render sweep"
    echo
    echo "- Total vignettes: $total"
    echo "- Passed: $pass"
    echo "- Failed: $fail"
    echo "- Skipped: $((total - pass - fail))"
    echo "- Elapsed: ${total_elapsed}s"
    echo "- Command: \`quarto render <name>.qmd --to gfm --execute\`"
    echo
    echo "## Results"
    echo
    echo "| Status | Vignette | Elapsed (s) | Notes |"
    echo "|--------|----------|-------------|-------|"
    for r in "${results[@]}"; do
        IFS='|' read -r status name elapsed note <<< "$r"
        printf "| %s | %s | %s | %s |\n" "$status" "$name" "$elapsed" "$note"
    done
} > "$REPORT"

echo "[sweep] wrote $REPORT"
echo "[sweep] $pass passed, $fail failed, ${total_elapsed}s"

if [ "$fail" -gt 0 ]; then
    exit 1
fi
exit 0
