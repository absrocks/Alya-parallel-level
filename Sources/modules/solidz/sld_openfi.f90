!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_openfi.f90
!> @author  Gerard Guillamet
!> @date    July, 2018
!> @brief   Open file names
!> @details
!>
!>          This subroutine gets ALL the file names and open them to be
!>          used by the module in two possible ways:\n
!>
!>          (a) Recalling them from the environment, when Alya is
!>              launched encapsulated in a shell script, or\n
!>          (b) Composing the names out of the problem name which is
!>              given as argument when the binary file Alya is
!>              launched "naked".\n
!>
!> @}
!------------------------------------------------------------------------

subroutine sld_openfi(itask)

use def_kintyp,  only : ip, rp
use def_master,  only : INOTSLAVE, modul, namda, exmod
use def_master,  only : kfl_naked, kfl_rstar
use mod_iofile,  only : iofile
use def_solidz,  only : kfl_psmat_sld, lun_psmat_sld
use def_solidz,  only : kfl_volca_sld, lun_carcy_res_sld, lun_carcy_cvg_sld

implicit none

integer(ip),   intent(in) :: itask
character(150)            :: fil_psmat_sld
character(150)            :: fil_carcy_res_sld
character(150)            :: fil_carcy_cvg_sld
character(7)              :: statu
character(11)             :: forma
character(6)              :: posit

    if ( INOTSLAVE ) then
        !
        ! Define unit opening option if this is a restart run
        !
        if(kfl_rstar == 2_ip) then
        statu='old'
        forma='formatted'
        posit='append'
        else
        statu='unknown'
        forma='formatted'
        posit='asis'
        end if

        select case (itask)

        case (1_ip)
            !
            ! Open files needed occasionally
            !
            if (kfl_naked==0) then
                call GET_ENVIRONMENT_VARIABLE('FOR114',fil_psmat_sld)
                call GET_ENVIRONMENT_VARIABLE('FOR115',fil_carcy_res_sld)
                call GET_ENVIRONMENT_VARIABLE('FOR116',fil_carcy_cvg_sld)
            else
                fil_psmat_sld = adjustl(trim(namda))//'-matrix.'//exmod(modul)//'.ps'
                fil_carcy_res_sld = adjustl(trim(namda))//'-cardiac-cycle.'//exmod(modul)//'.res'
                fil_carcy_cvg_sld = adjustl(trim(namda))//'-cardiac-cycle.'//exmod(modul)//'.cvg'
            end if
            !
            ! Matrix profile
            !
            if( kfl_psmat_sld > 0 ) then
                call iofile(0_ip,lun_psmat_sld,fil_psmat_sld,'SOLIDZ MATRIX',statu,forma,posit)
            end if
            !
            ! Cardia cycle
            !
            if( kfl_volca_sld > 0_ip ) then
                call iofile(0_ip,lun_carcy_res_sld,fil_carcy_res_sld,'SOLIDZ CARDIAC CYCLE',statu,forma,posit)
                call iofile(0_ip,lun_carcy_cvg_sld,fil_carcy_cvg_sld,'SOLIDZ CARDIAC CYCLE',statu,forma,posit)
            end if

        case (2_ip)
            !
            ! Close files
            !
            call iofile(2_ip,lun_psmat_sld,' ','SOLIDZ MATRIX PROFILE')

        end select

    end if

end subroutine sld_openfi

