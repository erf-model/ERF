#!/bin/bash

studydir='convergence'

run_case()
{
    Nz="$1"
    name="$2"
    dt="$3"
    max_step="$4"
    ./ekman_spiral inputs_ex amr.n_cell="4 4 $Nz" erf.fixed_dt="$dt" max_step = "$max_step" &> tmp.log
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

run_case 200 'dz25.0'     0.5            10
run_case 400 'dz12.5'     0.25           20
run_case 800 'dz6.25'     0.125          40
run_case 1600 'dz3.125'   0.0625         80
run_case 3200 'dz1.5625'  0.03125        160
run_case 6400 'dz0.78125' 0.015625       320

python check_convergence.py

