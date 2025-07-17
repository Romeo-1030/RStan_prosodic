#!/bin/bash

N_INTERVALS=3
OVERLAP=0.4
OUT_DIR=../data/cluster_result

SEEDS=(642 647 672 604 645 614 680 681 657 629)

for SEED in "${SEEDS[@]}"
do
  echo "Running Mapper with seed: $SEED"
  python run_mapper.py --seed $SEED --n_intervals $N_INTERVALS --overlap $OVERLAP --out_dir $OUT_DIR/seed${SEED}
done

