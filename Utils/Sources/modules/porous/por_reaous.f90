!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_reanut.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Read output strategy
!> @details Read output strategy
!> @} 
!------------------------------------------------------------------------
subroutine por_reaous()
  use def_parame
  use def_inpout
  use def_master
  use def_porous
  use def_domain
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: dummi

  if( INOTSLAVE ) then
     !
     ! Initializations
     !
     !
     ! Reach the section
     !
     call ecoute('por_reaous')
     do while(words(1)/='OUTPU')
        call ecoute('por_reaous')
     end do
     !
     ! Begin to read data.
     !
     do while(words(1)/='ENDOU')
        call ecoute('por_reaous')

        call posdef(2_ip,dummi)

        if(words(1)=='OUTPU') then
           !
           ! Output
           !
!           if(words(2)=='ERROR') then
              !
              ! Exact solution
              !
!              kfl_exacs_por=getint('SOLUT',1_ip,'#Exact solution')
!              expar_por=param(4:3+nexap_por)

!           else if(words(2)=='MATRI') then
              !
              ! Matrix profile
              !
!              kfl_psmat_por=getint('ITERA',1_ip,'#Iteration to postprocess matrix') 
 
!           else if(words(2)=='SOLVE') then
              !
              ! Solver convergence
              !
!              solve(1)%kfl_cvgso=1

!           end if

        end if

     end do

  end if

end subroutine por_reaous
    
