#!/bin/bash
#
# This script runs a sweep of canonical ABL cases with different MOST settings for manual
# regression testing
#
BUILDDIR='.'
exe=$BUILDDIR/ERF3d.gnu.TEST.MPI.ex
comp=''

comp='/projects/erf/equon/ERF/MyBuild_serial/install/bin/amrex_fcompare'
refdir='0_init'

#==============================================================================
maxstep=10
common_params="max_step=$maxstep erf.plot_int_1=$maxstep erf.check_int=-1"

declare -A landsea_params
landsea_params['land']="erf.is_land=1 erf.most.z0=0.16"
landsea_params['sea0']="erf.is_land=0 erf.most.z0=1e-4 erf.most.roughness_type_sea=charnock"
landsea_params['sea1']="erf.is_land=0 erf.most.z0=1e-4 erf.most.roughness_type_sea=modified_charnock"
landsea_params['sea2']="erf.is_land=0 erf.most.z0=1e-4 erf.most.roughness_type_sea=donelan"
landsea_params['sea3']="erf.is_land=0 erf.most.z0=1e-4 erf.most.roughness_type_sea=coare3.0"

declare -A surf_params
surf_params['adiabatic']="erf.most.surf_temp_flux=0"
surf_params['heatflux']="erf.most.surf_temp_flux=0.2"
surf_params['coolflux']="erf.most.surf_temp_flux=-0.05"
surf_params['heatrate']="erf.most.surf_temp=300 erf.most.surf_heating_rate=1.0"
surf_params['coolrate']="erf.most.surf_temp=300 erf.most.surf_heating_rate=-0.25"
#==============================================================================


set -e
(cd $BUILDDIR && make -j64 COMP=gnu USE_FFT=TRUE TINY_PROFILE=FALSE)
set +e


for landsea in "${!landsea_params[@]}"; do
    for surf in "${!surf_params[@]}"; do
        desc="${landsea}_${surf}"
        allparams="${landsea_params["$landsea"]} ${surf_params["$surf"]} ${common_params} erf.plot_file_1=$desc"
        echo '=============================================================================='
        echo "Running $desc ($allparams)..."
        srun -n ${SLURM_NTASKS} $exe inputs_anel_most $allparams $* &> log.$desc

        #tail log.$desc

        if [ -n "$comp" ]; then
            pltfile=`printf "${desc}%05d" $maxstep`
            $comp $pltfile $refdir/$pltfile 2>&1 | tee log.fcompare.$desc
        fi

        sleep 1  # workaround for slurm scheduler hangup?
    done
done
