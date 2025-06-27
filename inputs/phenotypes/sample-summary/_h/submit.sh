#!/bin/bash

head -n 1 ../../_m/phenotypes-AA.tsv | tr '\t' '\n' | nl > variables.tsv
