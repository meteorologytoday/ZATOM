#!/bin/bash

N=2

(
for f in $( ls -r submit_*.sh ); do

    #((i=i%N)); ((i++==0)) && wait
    sbatch $f

done
)


wait

echo "Done"
