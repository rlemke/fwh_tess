
mkdir -p ./tl
for s in $(seq 1 96); do
  S=$(printf "%03d" $s)
  f=./tl/S$S.txt
  [ -s $f ] && continue
  for v in 1 2 3; do
    curl -sfL -o $f https://tess.mit.edu/public/target_lists/2m/all_targets_S${S}_v$v.txt && break
  done
done
ls ./tl | wc -l
