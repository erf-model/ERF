!
! Flush the Fortran runtime's standard output and standard error units.
!
! Called from the exit trap in main.cpp, which ends the process with std::_Exit and so
! skips the rest of the normal exit sequence -- including the Fortran runtime's own
! cleanup, which is what writes out its unit buffers. Without this, a Fortran library
! that WRITEs its diagnostics and then STOPs (Noah-MP's energy-budget check, for
! example) loses those lines whenever standard output is redirected to a file. Compiled
! only when ERF builds Fortran (ERF_HAS_FORTRAN); main.cpp has a no-op otherwise.
!
subroutine erf_flush_fortran_units () bind(C, name="erf_flush_fortran_units")
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none
  flush(output_unit)
  flush(error_unit)
end subroutine erf_flush_fortran_units
