#!/usr/bin/env bash
# Pack 10,000 sample dirs into 10 tar.gz shards (1000 samples each).
# Uses pigz for parallel compression.
set -euo pipefail

SRC=${SRC:-/root/55miner/traindatasets_second}
OUT=${OUT:-/root/55miner/traindatasets_second_shards}
SHARD_SIZE=${SHARD_SIZE:-1000}

mkdir -p "$OUT"
cd "$SRC"

# Determine total samples
total=$(ls -d sample_* | wc -l)
echo "Sharding $total samples into chunks of $SHARD_SIZE → $OUT"

shard_idx=1
for start in $(seq 1 "$SHARD_SIZE" "$total"); do
    end=$(( start + SHARD_SIZE - 1 ))
    [[ $end -gt $total ]] && end=$total

    # Build list of sample dirs in this shard range
    list=$(mktemp)
    for i in $(seq -f "sample_%05g" "$start" "$end"); do
        echo "$i" >> "$list"
    done

    shard_name=$(printf "samples_%05d-%05d.tar.gz" "$start" "$end")
    out_path="$OUT/$shard_name"

    if [[ -f "$out_path" ]]; then
        echo "  [$shard_idx] $shard_name already exists, skipping"
    else
        echo "  [$shard_idx] packing $shard_name ($((end-start+1)) samples)"
        tar --use-compress-program="pigz -p 4" -cf "$out_path" -T "$list"
    fi
    rm -f "$list"
    shard_idx=$((shard_idx + 1))
done

echo
echo "=== shard summary ==="
ls -lh "$OUT" | tail -20
echo
echo "total: $(du -sh "$OUT" | cut -f1)"
