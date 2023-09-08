subroutine nsa_relvar(ktask,densx,presx,cdenx,enthx,enerx,umomx,velox,krest,ipoin)
!-----------------------------------------------------------------------
!****f* Nastal/nsa_setvar
! NAME 
!    nsa_setvar
! DESCRIPTION
!    This routine computes derived fields out of the primitive variables
!    and state equations
! USED BY
!    nsa_iniunk
!***
!-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_parame

  use      def_nastal

  implicit none
  integer(ip) :: ktask,idime,itenr,ntenr,kflnr,krest,ipoin
  real(rp)    :: densx,presx,cdenx,enthx,enerx,spien,wlore,vsqua,xauxi,&
       umomx(ndime),velox(ndime),resnr,tolnr,oldpx,newpx,usqua,pmini,pmaxi,auxi1,fpval,dfdpv

  krest= 0
  if (kfl_relat_nsa == 0)  return

  if (ktask==1) then
     !
     ! initial condition: conservative (relativistic) from primitive
     !
     vsqua= velox(1)*velox(1)+velox(2)*velox(2)
     if (ndime==3) vsqua= vsqua+velox(3)*velox(3)
!!!     spien= (enerx - 0.5 * vsqua) / densx         ! specific internal energy
     spien= presx / densx / (adgam_nsa - 1.0_rp)     ! specific internal energy
     wlore= sqrt(1.0_rp - vsqua)                  
     wlore= 1.0_rp / wlore                        ! lorentz factor
     cdenx= densx                                 ! move classical density from densi to cdens
     densx= cdenx * wlore                         ! rest mass density
     enthx= 1.0_rp + spien + presx / cdenx        ! enthalpy
     enerx= densx * wlore * enthx - presx - densx ! (rest mass) energy density
     umomx(1)= densx * wlore * enthx * velox(1)   ! (rest mass) momentum 
     umomx(2)= densx * wlore * enthx * velox(2)   
     if (ndime==3) umomx(3)= densx * wlore * enthx * velox(3)

     if (ipoin == 2000) then
        write(6,*) 'izquierda'
     else if (ipoin == 2001) then
        write(6,*) 'derecha'
     end if

  else if (ktask==2) then
     !
     ! primitive from conservative (relativistic)
     !
     tolnr= 1.0e-8
     ntenr= 300

     usqua= umomx(1)*umomx(1)+umomx(2)*umomx(2)
     if (ndime==3) usqua= vsqua+umomx(3)*umomx(3)
     usqua= sqrt(usqua)
     
     pmini= usqua-enerx-densx
     pmaxi= enerx*(adgam_nsa - 1.0_rp)

!     if (pmini > pmaxi) then
!        write(6,*) 'chau'
!        stop
!     end if

     if (pmini > presx) presx= pmini

     kflnr= 0
     itenr= 0

     if (ipoin == 2000) then
        write(6,*) 'izquierda'
        write(6,*) presx,densx,enerx
        write(6,*) velox
        write(6,*) umomx
     else if (ipoin == 2001) then
        write(6,*) 'derecha'
     end if

     do while (kflnr == 0)
        itenr= itenr+1
        xauxi= 1.0_rp / (enerx + densx + presx)
        velox(1)= umomx(1)*xauxi
        velox(2)= umomx(2)*xauxi
        if (ndime == 3) velox(3)= umomx(3)*xauxi
        vsqua= velox(1)*velox(1)+velox(2)*velox(2)
        if (ndime==3) vsqua= vsqua+velox(3)*velox(3)
        
        wlore= sqrt(1.0_rp - vsqua)                  
        wlore= 1.0_rp / wlore                        ! lorentz factor        
        auxi1= enerx + densx * (1.0_rp - wlore) + presx
        fpval= (adgam_nsa - 1.0_rp) * auxi1 / wlore / wlore - adgam_nsa * presx
        dfdpv= (adgam_nsa - 1.0_rp) * auxi1 * xauxi / wlore / wlore - 1.0_rp

        oldpx= presx 
        newpx= presx - fpval / dfdpv
        
        presx= pmini
        if (newpx > pmini) presx= newpx

        resnr= abs(1.0_rp - presx/oldpx) 

        if (itenr==ntenr) then
           write(6,*) 'cuidado' 
           write(6,*) presx,densx,enerx
           write(6,*) umomx(1:2),resnr
           krest= 1
           stop
           kflnr= 1
        end if

        if (resnr < tolnr) then
           kflnr= 1
        end if

     end do

     xauxi= 1.0_rp / (enerx + densx + presx)
     velox(1)= umomx(1)*xauxi
     velox(2)= umomx(2)*xauxi
     if (ndime == 3) velox(3)= umomx(3)*xauxi
     vsqua= velox(1)*velox(1)+velox(2)*velox(2)
     if (ndime==3) vsqua= vsqua+velox(3)*velox(3)     
     wlore= sqrt(1.0_rp - vsqua)                  
     wlore= 1.0_rp / wlore                        ! lorentz factor
     cdenx= densx / wlore             
     enthx= 1.0_rp + presx * adgam_nsa / cdenx / (adgam_nsa - 1.0_rp)

     if (ipoin == 2000) then
        write(6,*) 'izquierda'
        write(6,*) presx,cdenx
        write(6,*) velox
     else if (ipoin == 2001) then
        write(6,*) 'derecha'
     end if

  end if


end subroutine nsa_relvar
