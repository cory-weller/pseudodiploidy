#!/usr/bin/env bash

# get_read_set_counts

input=${1}

awk 'NR % 4 - 2 == 0' ${input} | \
    awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | \
    sort -nr