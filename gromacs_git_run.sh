#!/bin/bash -e
THREADS=24
ALIGNMENTS=0 32 64

cd /home/gad4877/proj/bigfoot
git clone https://github.com/gadube/gromacs_example.git
cd gromacs_example
./run.sh
