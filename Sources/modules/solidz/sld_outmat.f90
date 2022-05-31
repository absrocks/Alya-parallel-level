!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_outmat.f90
!> @author  Gerard Guillamet
!> @date    July, 2018
!>          - Subroutine written
!> @brief   Post-process output global system matrix
!> @details
!>          \verbatim
!>          The global system matrix can be postprocess in the following
!>          formats:
!>            - Postcript format (.ps)
!>            - Data File (.dat)
!>
!>          \endverbatim
!>
!> @}
!------------------------------------------------------------------------

subroutine sld_outmat()

  use def_kintyp,     only : ip, rp
  use def_master,     only : ISEQUEN
  use def_master,     only : itinn, ittim
  use def_master,     only : modul, title, namod
  use def_master,     only : solve
  use def_master,     only : amatr
  use def_domain,     only : npoin, c_dom, r_dom
  use mod_matrix,     only : matrix_output_dense_format
  use def_solidz,     only : kfl_psmat_sld, psmat_sld, lun_psmat_sld

  implicit none

  !
  ! Matrix output
  !
  if( ISEQUEN ) then
     if ( ittim == int(psmat_sld(1),ip) .and. itinn(modul) == int(psmat_sld(2),ip) ) then

        ! Matrix in Postcript format
        call pspltm(&
             npoin,npoin,solve(1)%ndofn,0_ip,c_dom,r_dom,amatr,&
             trim(title)//': '//namod(modul),0_ip,18.0_rp,'cm',&
             0_ip,0_ip,2_ip,lun_psmat_sld)
        call sld_openfi(2_ip)

        ! Matrix output file
        call matrix_output_dense_format(npoin,solve(1)%ndofn,r_dom,c_dom,amatr)
        kfl_psmat_sld = 0_ip

     end if
  end if

end subroutine sld_outmat
