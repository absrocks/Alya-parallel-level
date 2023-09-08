
 subroutine hlm_dirbcsdiffAssem(indvars)

  !---------------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_dirbcsdiff.f90
  ! NAME 
  !    hlm_dirbcsdiff
  ! DESCRIPTION
  !    This routine imposes Dirichlet boundary conditions on the system equations.
  ! USES
  ! USED BY
  !    hlm_matrix
  !---------------------------------------------------------------------------------

  use def_parame
  use def_master
  use def_helmoz
  use def_domain

  implicit none

  integer(ip), intent(in)  :: indvars
  integer(ip) :: ipoin,poin,ii,ih

  ! Modify the matrix columns
  do ii = 1,nzdom
      if (kfl_fixno_hlm(1,c_dom(ii)) == 1_ip) then       !If the element is in the 'boundary' column
          ih = 16_ip * (ii - 1_ip)
          damatx(ih+1)  = (0.0_rp,0.0_rp)
          !amatx(ih+2)  = (0.0_rp,0.0_rp)
          !amatx(ih+3)  = (0.0_rp,0.0_rp)
          damatx(ih+4)  = (0.0_rp,0.0_rp)
          !amatx(ih+5)  = (0.0_rp,0.0_rp)
          damatx(ih+6)  = (0.0_rp,0.0_rp)
          !amatx(ih+7)  = (0.0_rp,0.0_rp)
          damatx(ih+8)  = (0.0_rp,0.0_rp)
          !amatx(ih+9)  = (0.0_rp,0.0_rp)
          !amatx(ih+10) = (0.0_rp,0.0_rp)
          damatx(ih+11) = (0.0_rp,0.0_rp)
          damatx(ih+12) = (0.0_rp,0.0_rp)
          damatx(ih+13) = (0.0_rp,0.0_rp)
          damatx(ih+14) = (0.0_rp,0.0_rp)
          damatx(ih+15) = (0.0_rp,0.0_rp)
          damatx(ih+16) = (0.0_rp,0.0_rp)
          
      endif
  enddo

  ! Modify the matrix rows and RHS
  do ipoin = 1,npoin
      if (kfl_fixno_hlm(1,ipoin) == 1_ip) then
          poin = 4_ip * (ipoin - 1_ip)
          !drhsi(poin+1,indvars) = (0.0_rp,0.0_rp)
          !drhsi(poin+2,indvars) = (0.0_rp,0.0_rp)
          !drhsi(poin+3,indvars) = (0.0_rp,0.0_rp)
          !drhsi(poin+4,indvars) = (0.0_rp,0.0_rp)
          drhsix(poin+1) = (0.0_rp,0.0_rp)
          drhsix(poin+2) = (0.0_rp,0.0_rp)
          drhsix(poin+3) = (0.0_rp,0.0_rp)
          drhsix(poin+4) = (0.0_rp,0.0_rp)

          do ii = r_dom(ipoin),r_dom(ipoin+1)-1_ip
              if (c_dom(ii) == ipoin) then
                  ih = 16_ip * (ii - 1_ip)
                  damatx(ih+1)  = (0.0_rp,0.0_rp)    !<-- differentiated
                  !amatx(ih+2)  = (0.0_rp,0.0_rp)
                  !amatx(ih+3)  = (0.0_rp,0.0_rp)
                  damatx(ih+4)  = (0.0_rp,0.0_rp)
                  !amatx(ih+5)  = (0.0_rp,0.0_rp)
                  damatx(ih+6)  = (0.0_rp,0.0_rp)    !<-- differentiated 
                  !amatx(ih+7)  = (0.0_rp,0.0_rp)
                  damatx(ih+8)  = (0.0_rp,0.0_rp)
                  !amatx(ih+9)  = (0.0_rp,0.0_rp)
                  !amatx(ih+10) = (0.0_rp,0.0_rp)
                  damatx(ih+11) = (0.0_rp,0.0_rp)    !<-- differentiated  
                  damatx(ih+12) = (0.0_rp,0.0_rp)
                  damatx(ih+13) = (0.0_rp,0.0_rp)
                  damatx(ih+14) = (0.0_rp,0.0_rp)
                  damatx(ih+15) = (0.0_rp,0.0_rp)
                  damatx(ih+16) = (0.0_rp,0.0_rp)    !<-- differentiated   
              else 
                  ih = 16_ip * (ii - 1_ip)
                  damatx(ih+1)  = (0.0_rp,0.0_rp)
                  !amatx(ih+2)  = (0.0_rp,0.0_rp)
                  !amatx(ih+3)  = (0.0_rp,0.0_rp)
                  damatx(ih+4)  = (0.0_rp,0.0_rp)
                  !amatx(ih+5)  = (0.0_rp,0.0_rp)
                  damatx(ih+6)  = (0.0_rp,0.0_rp)
                  !amatx(ih+7)  = (0.0_rp,0.0_rp)
                  damatx(ih+8)  = (0.0_rp,0.0_rp)
                  !amatx(ih+9)  = (0.0_rp,0.0_rp)
                  !amatx(ih+10) = (0.0_rp,0.0_rp)
                  damatx(ih+11) = (0.0_rp,0.0_rp)
                  damatx(ih+12) = (0.0_rp,0.0_rp)
                  damatx(ih+13) = (0.0_rp,0.0_rp)
                  damatx(ih+14) = (0.0_rp,0.0_rp)
                  damatx(ih+15) = (0.0_rp,0.0_rp)
                  damatx(ih+16) = (0.0_rp,0.0_rp)
              endif
          enddo
      endif
  enddo


end subroutine hlm_dirbcsdiffAssem
