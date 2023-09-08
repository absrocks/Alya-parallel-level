subroutine rad_reaous()

  !-----------------------------------------------------------------------
  !
  ! This routine reads the output strategy for the radiation module
  ! 
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_radiat
  use def_domain
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: dummi

  if( INOTSLAVE ) then
     !
     ! Initializations
     !
     kfl_splot_rad          = 0           ! Surface plot
     kfl_psmat_rad          = 0           ! Postscript file of the matrix
     kfl_exacs_rad          = 0           ! Exact solution
     npp_bound_rad          = 0           ! Postprocess boundary conditions

     !
     ! Reach the section
     !
     call ecoute('rad_reaous')
     do while(words(1)/='OUTPU')
        call ecoute('rad_reaous')
     end do
     !
     ! Begin to read data.
     !
     do while(words(1)/='ENDOU')
        call ecoute('rad_reaous')

        call posdef(2_ip,dummi)

        if(words(1)=='OUTPU') then
           !
           ! Output
           !
           if(words(2)=='ERROR') then
              !
              ! Exact solution
              !
              kfl_exacs_rad=getint('SOLUT',1_ip,'#Exact solution')
              expar_rad=param(4:3+nexap_rad)

           else if(words(2)=='MATRI') then
              !
              ! Matrix profile
              !
              kfl_psmat_rad=getint('ITERA',1_ip,'#Iteration to postprocess matrix') 
 
           else if(words(2)=='SOLVE') then
              !
              ! Solver convergence
              !
              solve(1)%kfl_cvgso=1

           end if

        else if(words(1)=='SURFA') then
           !
           ! Draw 3D plots
           !
           kfl_splot_rad = 1

        end if

     end do

  end if

end subroutine rad_reaous
    
