#!/bin/bash

for f in $( ls submit_*.sh ); do

    sbatch $f

done
