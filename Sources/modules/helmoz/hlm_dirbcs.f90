 subroutine hlm_dirbcs()

  !---------------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_dirbcs.f90
  ! NAME 
  !    hlm_dirbcs
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

  integer(ip) :: ipoin,poin,ii,ih

  ! Modify the matrix columns
  do ii = 1,nzdom
      if (kfl_fixno_hlm(1,c_dom(ii)) == 1_ip) then       !If the element is in the 'boundary' column
          ih = 16_ip * (ii - 1_ip)
          amatx(ih+1)  = (0.0_rp,0.0_rp)
          !amatx(ih+2)  = (0.0_rp,0.0_rp)
          !amatx(ih+3)  = (0.0_rp,0.0_rp)
          amatx(ih+4)  = (0.0_rp,0.0_rp)
          !amatx(ih+5)  = (0.0_rp,0.0_rp)
          amatx(ih+6)  = (0.0_rp,0.0_rp)
          !amatx(ih+7)  = (0.0_rp,0.0_rp)
          amatx(ih+8)  = (0.0_rp,0.0_rp)
          !amatx(ih+9)  = (0.0_rp,0.0_rp)
          !amatx(ih+10) = (0.0_rp,0.0_rp)
          amatx(ih+11) = (0.0_rp,0.0_rp)
          amatx(ih+12) = (0.0_rp,0.0_rp)
          amatx(ih+13) = (0.0_rp,0.0_rp)
          amatx(ih+14) = (0.0_rp,0.0_rp)
          amatx(ih+15) = (0.0_rp,0.0_rp)
          amatx(ih+16) = (0.0_rp,0.0_rp)
          
      endif
  enddo

  ! Modify the matrix rows and RHS
  do ipoin = 1,npoin
      if (kfl_fixno_hlm(1,ipoin) == 1_ip) then
          poin = 4_ip * (ipoin - 1_ip)
          rhsix(poin+1) = (0.0_rp,0.0_rp)
          rhsix(poin+2) = (0.0_rp,0.0_rp)
          rhsix(poin+3) = (0.0_rp,0.0_rp)
          rhsix(poin+4) = (0.0_rp,0.0_rp)

          do ii = r_dom(ipoin),r_dom(ipoin+1)-1_ip
              if (c_dom(ii) == ipoin) then	
                  ih = 16_ip * (ii - 1_ip)
                  amatx(ih+1)  = (1.0_rp,0.0_rp)
                  !amatx(ih+2)  = (0.0_rp,0.0_rp)
                  !amatx(ih+3)  = (0.0_rp,0.0_rp)
                  amatx(ih+4)  = (0.0_rp,0.0_rp)
                  !amatx(ih+5)  = (0.0_rp,0.0_rp)
                  amatx(ih+6)  = (1.0_rp,0.0_rp)
                  !amatx(ih+7)  = (0.0_rp,0.0_rp)
                  amatx(ih+8)  = (0.0_rp,0.0_rp)
                  !amatx(ih+9)  = (0.0_rp,0.0_rp)
                  !amatx(ih+10) = (0.0_rp,0.0_rp)
                  amatx(ih+11) = (1.0_rp,0.0_rp)
                  amatx(ih+12) = (0.0_rp,0.0_rp)
                  amatx(ih+13) = (0.0_rp,0.0_rp)
                  amatx(ih+14) = (0.0_rp,0.0_rp)
                  amatx(ih+15) = (0.0_rp,0.0_rp)
                  amatx(ih+16) = (1.0_rp,0.0_rp)
              else 
                  ih = 16_ip * (ii - 1_ip)
                  amatx(ih+1)  = (0.0_rp,0.0_rp)
                  !amatx(ih+2)  = (0.0_rp,0.0_rp)
                  !amatx(ih+3)  = (0.0_rp,0.0_rp)
                  amatx(ih+4)  = (0.0_rp,0.0_rp)
                  !amatx(ih+5)  = (0.0_rp,0.0_rp)
                  amatx(ih+6)  = (0.0_rp,0.0_rp)
                  !amatx(ih+7)  = (0.0_rp,0.0_rp)
                  amatx(ih+8)  = (0.0_rp,0.0_rp)
                  !amatx(ih+9)  = (0.0_rp,0.0_rp)
                  !amatx(ih+10) = (0.0_rp,0.0_rp)
                  amatx(ih+11) = (0.0_rp,0.0_rp)
                  amatx(ih+12) = (0.0_rp,0.0_rp)
                  amatx(ih+13) = (0.0_rp,0.0_rp)
                  amatx(ih+14) = (0.0_rp,0.0_rp)
                  amatx(ih+15) = (0.0_rp,0.0_rp)
                  amatx(ih+16) = (0.0_rp,0.0_rp)
              endif
          enddo
      endif
  enddo

end subroutine hlm_dirbcs
