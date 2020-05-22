#!/bin/bash

for i in {0..21}; do
    INDEX=$(printf "%02d" $i)
    IP=$(expr $i + 1)
    TIME=$(head -n $IP times.txt | tail -n 1)
    touch -t 010100$TIME slides-$INDEX.png
done
