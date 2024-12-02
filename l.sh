#!/usr/bin/env bash

for it in $(seq 1 30); do
  prev=$((it - 1))
  time poetry run python mfl/main.py -- \
    -p ./populations/"$(printf "%03d" "$prev")".csv \
    -d ./populations \
    -m ./meta.json \
    -r "$it" \
    -s 100 \
    --max-added 10 \
    --max-replaced 10 \
    --max-removed 10 \
    --max-crossed 40 \
    --bank-size 100

  if test $? -ne 0; then
    exit 1
  fi
done
