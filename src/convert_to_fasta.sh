#!/usr/bin/env bash

input=${1}

awk '{ print ">"NR"-"$1"\n"$2 }' ${input} | fold -w 80