!-----------------------------------------------------------------------
!> @addtogroup Optsol
!> @{
!> @file    opt_reapro.f90
!> @author  Oscar Peredo
!> @date    October 20, 2013
!> @brief   Read problem data
!> @details Read problem data
!> @} 
!-----------------------------------------------------------------------
subroutine opt_reapro()
  use      def_parame
  use      def_master
  use      def_optsol
  use      def_inpout
  use mod_ecoute, only :  ecoute
  implicit none
  character(10) :: messa

  !kfl_goopt=1 is setted in kernel/master/inirun.f90
  kfl_first_opt=1
  kfl_curstp_opt=1
  kfl_curlin_opt=1
  kfl_curres_opt=0
  costf_prev=1e+100_rp
  stepj=1.0_rp
  kfl_steplength_opt=1.0_rp

  if(INOTSLAVE) then
  !
  ! Reach the section
  !
  servi=ID_OPTSOL 
  messa='opt_reapro'

  rewind(lisda)

  do while(words(1)/='RUNDA')
     call ecoute(messa)
  end do
  do while(words(1)/='ENDRU')
     call ecoute(messa)
  end do
  do while(words(1)/='PROBL')
     call ecoute(messa)
  end do
  !
  ! Read data
  !      
  do while(words(1)/='ENDPR')
     call ecoute(messa)
     if(words(1)=='OPTSO') then 
        if(exists('ON   ')) then
           kfl_servi(servi)=1
           do while(words(1)/='ENDOP')
              call ecoute(messa)
              if(words(1)=='ADJOI') then
                if(words(2)=='YES  ') then
                  kfl_method_opt=1
                end if
              end if
              if(words(1)=='CONVE') then
                kfl_conver_opt=param(1)
              end if
              if(words(1)=='BARRIE') then
                kfl_conbar_opt=param(1)
              end if
              if(words(1)=='MAXST') then
                kfl_maxstp_opt=int(param(1))
              end if
              if(words(1)=='MAXLI') then
                kfl_maxlin_opt=int(param(1))
              end if
              if(words(1)=='MAXBAR') then
                kfl_maxbar_opt=int(param(1))
              end if
              if(words(1)=='MAXRES') then
                kfl_maxres_opt=int(param(1))
              end if
              if(words(1)=='SCHEM') then
                kfl_scheme_opt=int(param(1))
              end if
              if(words(1)=='NUMDE') then
                kfl_ndvars_opt=int(param(1))
              end if
              if(words(1)=='SCALE') then
                kfl_scale_opt=param(1)
              end if
              if(words(1)=='STEPL') then
                kfl_steplength_opt=param(1)
              end if

           end do
        else
           kfl_servi(servi)=0
        end if
     end if
  end do
  end if
end subroutine opt_reapro
