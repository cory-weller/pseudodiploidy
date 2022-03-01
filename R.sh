#!/usr/bin/env bash

singularity exec --bind ${PWD} src/pseudodiploidy.sif R
