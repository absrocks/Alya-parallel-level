subroutine rad_outinf()

  !-----------------------------------------------------------------------
  !
  ! This routine writes on the radiation files
  !
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_radiat
  use mod_outfor, only : outfor
  implicit none
  character(60) :: equat
  integer(ip)   :: ierhs,imate
  character(2)  :: wcpcv

  if(IMASTER) then
!!$     !
!!$     ! Write information in Result file
!!$     !
!!$     if(kfl_rstar/=2) then
!!$        equat=''
!!$        if(kfl_regim_rad==1) then
!!$           wcpcv='Cv'
!!$        else
!!$           wcpcv='Cp'
!!$        end if
!!$        if(kfl_timei_rad==1) equat=trim(equat)//'rho*'//wcpcv//'*dT/dt '
!!$        if(kfl_advec_rad>=1) equat=trim(equat)//' + rho*'//wcpcv//'*u.grad(T) '
!!$
!!$        if(react_rad>zetem)  equat=trim(equat)//' + s*T'
!!$        if(kfl_regim_rad==1) equat=trim(equat)//' + [rho*R*div(u)]*T'
!!$        equat=trim(equat)//'= '
!!$        ierhs=0
!!$        if(kfl_sourc_rad==1) then
!!$           equat=trim(equat)//' Q '
!!$           ierhs=1
!!$        end if
!!$        if(kfl_exacs_rad/=0) then
!!$           equat=trim(equat)//' Q_exact '
!!$           ierhs=1           
!!$        end if
!!$        if(kfl_regim_rad==3) then
!!$           equat=trim(equat)//' + dp0/dt '
!!$           ierhs=1                      
!!$        end if
!!$        if(ierhs==0) equat=trim(equat)//' 0 '
!!$
!!$        call outfor(25_ip,momod(modul)%lun_outpu,'DIFFERENTIAL EQUATION')
!!$        write(momod(modul)%lun_outpu,110) trim(equat),nmate
!!$        write(momod(modul)%lun_outpu,111) 'rho [ M / L^3 ]     = ',(densi_rad(1,imate),imate=1,nmate)
!!$        write(momod(modul)%lun_outpu,111) 'Cp  [ L^2 / T^3 K ] = ',(sphea_rad(1,imate),imate=1,nmate)
!!$        write(momod(modul)%lun_outpu,111) 'k   [ M L / T^3 K ] = ',(tcond_rad(1,imate),imate=1,nmate)
!!$        if(kfl_regim_rad==1) then
!!$           write(momod(modul)%lun_outpu,111) 'R   [ L^2 / T^2 K ] = ',gasco_rad
!!$           write(momod(modul)%lun_outpu,111) 'Cv  [ L^2 / T^3 K ] = ',(sphea_rad(1,imate)-gasco_rad,imate=1,nmate)  
!!$        end if
!!$        if(react_rad>zetem) then
!!$           write(momod(modul)%lun_outpu,111) 's   [ M / L T^3 K ] = ',react_rad
!!$        end if        
!!$        if(kfl_cotur_rad/=0) then
!!$           write(momod(modul)%lun_outpu,112) 'Prt                 = ',prtur_rad,' (kt=rho*cp*nut/Prt)'
!!$        end if
!!$        write(momod(modul)%lun_outpu,*) 
!!$     end if
     print *, 'RADIAT Does not write information yet.'

  end if
  !
  ! Formats
  !
110 format(/,&
       & 10x,a,//,&
       & 10x,'Physical properties: ',i2,' material(s)',/,&
       & 10x,'------------------------------------------------------------ ')
111 format(&
       & 10x,a,10(e12.6,1x))
112 format(&
       & 10x,a,e12.6,a)

end subroutine rad_outinf

