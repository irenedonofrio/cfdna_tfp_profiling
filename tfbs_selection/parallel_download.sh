#!/usr/bin/env bash
# Parallel, resumable download of one large file using HTTP range requests.
# usage: URL=... SIZE=... N=8 OUT=file.zip bash parallel_download.sh
set -u
: "${URL:?set URL}" "${SIZE:?set SIZE (bytes)}"
N=${N:-8}
OUT=${OUT:-Homo_sapiens_meta_clusters.zip}

# The server must honour ranges, otherwise every piece would get the whole file.
code=$(curl -s -L -r 0-0 -o /dev/null -w '%{http_code}' "$URL")
if [ "$code" != "206" ]; then
  echo "server answered $code to a range request (need 206): parallel download not possible"; exit 1
fi

CH=$(( (SIZE + N - 1) / N ))

get_part() {   # $1 = piece index, $2 = first byte, $3 = last byte
  local f="$OUT.part$1" want=$(( $3 - $2 + 1 )) have
  while :; do
    if [ -f "$f" ]; then have=$(( $(wc -c < "$f") )); else have=0; fi
    if [ "$have" -gt "$want" ]; then          # e.g. a reused download that ran past this piece
      head -c "$want" "$f" > "$f.tmp" && mv "$f.tmp" "$f"; have=$want
    fi
    [ "$have" -ge "$want" ] && break
    curl -s -L --speed-limit 1000 --speed-time 60 -r "$(( $2 + have ))-$3" "$URL" >> "$f" \
      || sleep 10
  done
}

i=0
while [ $i -lt $N ]; do
  s=$(( i * CH )); e=$(( s + CH - 1 )); [ $e -ge $SIZE ] && e=$(( SIZE - 1 ))
  get_part $i $s $e &
  i=$(( i + 1 ))
done

# progress every 30 s until all pieces are done
while [ -n "$(jobs -r)" ]; do
  sleep 30
  tot=0; for f in "$OUT".part*; do [ -f "$f" ] && tot=$(( tot + $(wc -c < "$f") )); done
  echo "$(date +%H:%M)  $(( tot / 1048576 )) / $(( SIZE / 1048576 )) MiB"
done
wait

: > "$OUT"
i=0; while [ $i -lt $N ]; do cat "$OUT.part$i" >> "$OUT"; i=$(( i + 1 )); done
got=$(( $(wc -c < "$OUT") ))
if [ "$got" -eq "$SIZE" ]; then
  echo "done: $OUT is $got bytes (expected $SIZE). You can delete $OUT.part*"
else
  echo "SIZE MISMATCH: got $got, expected $SIZE -- keep the parts and tell Claude"; exit 1
fi
