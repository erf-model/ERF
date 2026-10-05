#!/bin/bash
# Run both cases of the comparison and compare them (see README.md).
#
#   ./run_comparison.sh /path/to/erf_exec [nranks] [--barren | --barren-moist]
#
# The flags may come anywhere among the arguments; nranks defaults to 4.
#
# --barren runs the dry-surface variant: bare land (Noah-MP land-use 16, no vegetation) on the
# same soil at its wilting point, where both surfaces evaporate little. --barren-moist runs
# the same bare land at 0.25 m3/m3, where the bare-soil resistances set the evaporation.
#
# Needs an executable built with -DERF_ENABLE_NOAHMP=ON (which needs a parallel NetCDF),
# ncgen on the PATH, and Python 3 with numpy, matplotlib and yt. Each case runs in its own
# directory under ./runs (./runs_barren, ./runs_barren_moist); compare.py writes comparison.csv and
# comparison.png there.
set -euo pipefail

usage="usage: $0 /path/to/erf_exec [nranks] [--barren | --barren-moist] (flags may come anywhere)"
# The variant is a flag in any position; the rest are the executable and the rank count.
barren=""
args=()
for arg in "$@"; do
    case "$arg" in
        --barren|--barren-moist)
            [ -z "$barren" ] || { echo "give one variant; $usage"; exit 1; }
            barren="$arg" ;;
        -*) echo "unknown option $arg; $usage"; exit 1 ;;
        *) args+=("$arg") ;;
    esac
done
[ "${#args[@]}" -ge 1 ] && [ "${#args[@]}" -le 2 ] || { echo "$usage"; exit 1; }
exe=${args[0]}
nranks=${args[1]:-4}
case "$nranks" in
    ''|*[!0-9]*|0) echo "nranks must be a positive integer, got '$nranks'; $usage"; exit 1 ;;
esac
here=$(cd "$(dirname "$0")" && pwd)
erf_home=$(cd "$here/../../../.." && pwd)
table="$erf_home/Submodules/Noah-MP/parameters/NoahmpTable.TBL"
[ -f "$table" ] || { echo "Noah-MP's parameter table is missing: $table (initialise Submodules/Noah-MP)"; exit 1; }

land_option=""
force_restore_deck=inputs_force_restore
out="$here/runs"
if [ "$barren" = "--barren" ]; then
    land_option="--barren"
    force_restore_deck=inputs_force_restore_barren
    out="$here/runs_barren"
elif [ "$barren" = "--barren-moist" ]; then
    land_option="--barren --smois 0.25"
    force_restore_deck=inputs_force_restore_barren_moist
    out="$here/runs_barren_moist"
fi

mkdir -p "$out"
for case in noahmp force_restore; do
    dir="$out/$case"
    rm -rf "$dir"; mkdir -p "$dir"
    cp "$here/inputs_common" "$here/inputs_noahmp" "$here/inputs_force_restore_common" \
       "$here/inputs_force_restore" "$here/inputs_force_restore_barren" \
       "$here/inputs_force_restore_barren_moist" "$here/input_sounding" \
       "$here/namelist.erf" "$here/make_land_files.py" "$table" "$dir/"
    (cd "$dir" && python3 make_land_files.py $land_option > /dev/null)
    deck=inputs_noahmp
    if [ "$case" = force_restore ]; then deck=$force_restore_deck; fi
    echo "running $case ($deck) on $nranks ranks in $dir"
    (cd "$dir" && mpiexec -n "$nranks" "$exe" "$deck" > run.log 2>&1) || {
        echo "the $case run failed; see $dir/run.log"; exit 1; }
done
python3 "$here/compare.py" --noahmp "$out/noahmp" --force-restore "$out/force_restore" --out "$out"
