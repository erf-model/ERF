#!/bin/bash

studydir='convergence'

run_case()
{
    Nz="$1"
    name="$2"
    ./ekman_spiral inputs_ex amr.n_cell="4 4 $Nz" &> tmp.log
    lastline=`tail -n 1 tmp.log`
    if [[ "$lastline" != "AMReX"*"finalized" ]]; then
        echo "Case $name failed"
        exit 1
    else
        rm tmp.log
    fi
    mkdir -p "$studydir/$name"
    mv plt* $studydir/$name/
    echo "Case $name complete"
}

rm -rf $studydir

#geometry.prob_lo     =  0.  0.    0.
#geometry.prob_hi     = 50. 50. 5000.
run_case 200 'dz25.0'
run_case 400 'dz12.5'
run_case 800 'dz6.25'
run_case 1600 'dz3.125'
run_case 3200 'dz1.5625'
run_case 6400 'dz0.78125'

python check_convergence.py

