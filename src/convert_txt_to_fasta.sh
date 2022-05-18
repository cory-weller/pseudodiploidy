#!/usr/bin/env bash

awk '{{print NR,$0}}' ${1} |sed 's/^/>/g' | tr " " "\n" | fold