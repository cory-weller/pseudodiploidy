#!/usr/bin/env bash

singularity exec -B ~/.Xauthority:/home/.Xauthority -H ${PWD}:/home src/R.sif R
