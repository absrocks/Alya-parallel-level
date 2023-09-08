subroutine hlm_mlsint(tpoin,nn,ii,Asx,Asy,Asz,Psis)

  !-------------------------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_mlsint.f90
  ! NAME
  !    hlm_mlsint
  ! DESCRIPTION
  !    This routine performs the Moving Least-Squares Interpolation of 
  !    FE-computed scalar and vector secondary potentials at a test point
  !    (site) that moves over the domain. 
  ! INPUT ARGUMENTS
  !    TPOIN ... Cartesian coordinates of a test point (xt, yt, zt)
  !    NN ...... Number of closest mesh nodes to a test point, N 
  !    II ...... Number of a test point        
  ! OUTPUT ARGUMENTS
  !    Asx ..... Coefficients for interpolation of x component of vector secondary potential
  !    Asy ..... Coefficients for interpolation of y component of vector secondary potential
  !    Asz ..... Coefficients for interpolation of z component of vector secondary potential
  !    Psis .... Coefficients for interpolation of scalar secondary potential 
  ! USES
  !    hlm_closnod
  !    hlm_wls
  ! USED BY
  !    hlm_fivecs
  !-------------------------------------------------------------------------------------------

  use def_master
  use def_domain
  use def_helmoz
  use def_parame

  implicit none

  real(rp),    intent(in) :: tpoin(3)
  integer(ip), intent(in) :: nn,ii
  complex(rp), intent(out):: Asx(4),Asy(4),Asz(4),Psis(4)

  integer(ip) :: jj,kk
  integer(ip) :: clnod1(nn)
  real(rp)    :: clnod2(nn)         

  kk = nn * (ii - 1_ip)
  do jj = 1,nn
      clnod1(jj) = clnod1_hlm(kk+jj)
      clnod2(jj) = clnod2_hlm(kk+jj)
  enddo

  !WLS interpolation of vector and scalar secondary potentials
  call hlm_wls(tpoin,nn,ii,clnod1,clnod2,Asx,Asy,Asz,Psis)

end subroutine hlm_mlsint
