!=================================================================================================================
! WDM6 ISO Helper Module
! Helper functions for WDM6 microphysics using iso_c_binding
! Based on ERF_module_mp_wsm6_isohelper.F90
!=================================================================================================================

module wdm6_isohelper
  use iso_c_binding, only: c_double
  implicit none
  private
  public :: wdm6_rgmma

  integer, parameter :: kind_phys = c_double

contains

  !=================================================================================================================
  ! Gamma function using Weierstrass infinite product form
  ! Used for computing slope parameter coefficients
  !=================================================================================================================
  function wdm6_rgmma(x) result(rgmma_result)
    implicit none
    real(kind=kind_phys), intent(in) :: x
    real(kind=kind_phys) :: rgmma_result
    real(kind=kind_phys), parameter :: euler = 0.577215664901532_kind_phys
    real(kind=kind_phys) :: y
    integer :: i

    if (x == 1.0_kind_phys) then
      rgmma_result = 0.0_kind_phys
      return
    endif

    rgmma_result = x * exp(euler * x)
    do i = 1, 10000
      y = real(i, kind=kind_phys)
      rgmma_result = rgmma_result * (1.0_kind_phys + x/y) * exp(-x/y)
    enddo
    rgmma_result = 1.0_kind_phys / rgmma_result

  end function wdm6_rgmma

end module wdm6_isohelper
