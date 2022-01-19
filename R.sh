#!/usr/bin/env bash

singularity exec --bind ${PWD} src/R.sif R
