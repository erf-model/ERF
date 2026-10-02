#!/bin/bash
# Run both cases of the comparison and compare them (see README.md).
#
#   ./run_comparison.sh /path/to/erf_exec [nranks]
#
# Needs an executable built with -DERF_ENABLE_NOAHMP=ON (which needs a parallel NetCDF),
# ncgen on the PATH, and Python 3 with numpy and matplotlib. Each case runs in its own
# directory under ./runs; compare.py writes comparison.csv and comparison.png here.
set -euo pipefail

exe=${1:?usage: $0 /path/to/erf_exec [nranks]}
nranks=${2:-4}
here=$(cd "$(dirname "$0")" && pwd)
erf_home=$(cd "$here/../../../.." && pwd)
table="$erf_home/Submodules/Noah-MP/parameters/NoahmpTable.TBL"
[ -f "$table" ] || { echo "Noah-MP's parameter table is missing: $table (initialise Submodules/Noah-MP)"; exit 1; }

mkdir -p "$here/runs"
for case in noahmp force_restore; do
    dir="$here/runs/$case"
    rm -rf "$dir"; mkdir -p "$dir"
    cp "$here/inputs_common" "$here/inputs_$case" "$here/input_sounding" "$here/namelist.erf" \
       "$here/make_land_files.py" "$table" "$dir/"
    (cd "$dir" && python3 make_land_files.py > /dev/null)
    echo "running $case on $nranks ranks in $dir"
    (cd "$dir" && mpiexec -n "$nranks" "$exe" "inputs_$case" > run.log 2>&1) || {
        echo "the $case run failed; see $dir/run.log"; exit 1; }
done
python3 "$here/compare.py" --noahmp "$here/runs/noahmp" --force-restore "$here/runs/force_restore" \
        --out "$here"
