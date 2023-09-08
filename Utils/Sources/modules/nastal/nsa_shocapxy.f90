!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_shocapxy.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Compute shock capturing terms
!> @details Compute shock capturing terms
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_shocapxy(itask,&
     kapsh,graun,&
     resid,sound,cartd,chale,xvelo,velmo,taupa,difeq,&
     zesho,qufac,shmet,shote)
  use      def_kintyp
  use      def_domain, only:  ndime
  use      def_nastal, only:  shock_nsa,kfl_shock_nsa,ndofn_nsa
  implicit none
  integer(ip)           :: itask,inode,ipoin,idome,igaus,idime,jdime,idofn,jdofn

  real(rp), intent(in)  :: graun(ndofn_nsa,ndime),resid(ndofn_nsa),cartd(ndime),chale(2),xvelo(ndime)
  real(rp)              :: qufac,taupa(ndofn_nsa),difeq(ndofn_nsa),zesho,velmo,sound
  real(rp), intent(out) :: shmet(ndime,ndime,ndofn_nsa),shote(ndofn_nsa)

  integer(ip)           :: kapsh(ndofn_nsa),napsh
  real(rp)              :: &
       grau2(ndofn_nsa),grmod(ndofn_nsa),remod(ndofn_nsa),&
       ficve(ndofn_nsa),shpec,difau,shfau,zesh2,hconv,xmach
  real(rp)              :: epsdc,epssu,epssl,epsde,fiso(ndofn_nsa),faniso(ndofn_nsa),uvel,vvel,wvel,xvel2


! !$OMP

  if (itask == 1) then
     
     shmet = 0.0_rp

     if (kfl_shock_nsa(1) == 0) return

     hconv = chale(1)

     kapsh = 1             ! default: apply shock capturing
     xvel2 = velmo*velmo   ! velocity squared
     zesh2 = zesho*zesho
!     if (xvel2 < zesho) kapsh = 0            ! velocity is the first threshold  
     if (shock_nsa < zesho) then
        kapsh = 0
        return
     end if
     
     napsh = 0
     do idofn = 1,ndofn_nsa
        grau2(idofn)= 0.0_rp
        do idime = 1,ndime                                     
           grau2(idofn) = grau2(idofn) + graun(idofn,idime)*graun(idofn,idime)            !evaluate grad modul
        end do
        grmod(idofn) = sqrt(grau2(idofn))
        remod(idofn) =  abs(resid(idofn))
        if (grau2(idofn) < zesho) then
           kapsh(idofn) = 0
           ficve(idofn) = 0.0_rp
        else
           ficve(idofn) = remod(idofn) / grmod(idofn)
        end if
        xmach= velmo/sound
        if (ficve(idofn) < zesho) kapsh(idofn)= 0

        if (ficve(idofn) > zesho) then
!           write (6,200) idofn,ficve(idofn),sound,velmo
200        format (i4,10(2x,f8.4))
        end if
        napsh = napsh + kapsh(idofn)
     end do
     
     if (napsh == 0) return       ! no one needs shock capturing, return
     
     uvel= xvelo(1)
     vvel= xvelo(2)
     wvel= xvelo(ndime)

     do idofn=1,ndofn_nsa
        shfau = shock_nsa
        if (difeq(idofn) > zesh2) then
           shpec = ficve(idofn) * hconv / 2.0_rp / difeq(idofn)           
           if (shpec > 0.0_rp) shfau = shock_nsa - 1.0 / shpec
           if (shfau < 0.0_rp) shfau = 0.0_rp
        end if
        if (idofn == ndime+1) shfau = shock_nsa  ! no diffusion for continuity
        epsdc= 0.5_rp * shfau * hconv * ficve(idofn)      
        fiso(idofn)    = epsdc
        
        epssu = taupa(idofn) * xvel2

        epssl = epsdc - epssu
        if (epssl < 0.0_rp) then
           epssl= 0.0_rp
        end if
        if (xvel2 > 0.0_rp) then
           faniso(idofn)= (epssl - epsdc)/xvel2     
        else
           faniso(idofn)= 0.0_rp
        end if

!!!!! ojo!! probar esto:
!!!!  hago isotropo porque como no se calcular la verdadera contribucion de supg, directamente la quito
!!!!  y esto parece que va mejor al menos en el shock tube
!        if (idofn == ndime+1) then
!           faniso(idofn)  = 0.0_rp      ! isotropic, no convection for continuity eq
!        end if


        ! isotropic shock capturing
        if (kfl_shock_nsa(idofn) == 2) faniso(idofn)  = 0.0_rp      

        shmet(1,1,idofn) =  faniso(idofn)*uvel*uvel + fiso(idofn)
        shmet(2,2,idofn) =  faniso(idofn)*vvel*vvel + fiso(idofn)        
        shmet(1,2,idofn) =  faniso(idofn)*uvel*vvel
        shmet(2,1,idofn) =  shmet(1,2,idofn)

        if (ndime == 3) then           
           shmet(3,3,idofn)= faniso(idofn)*wvel*wvel + fiso(idofn)
           shmet(1,3,idofn)= faniso(idofn)*uvel*wvel
           shmet(2,3,idofn)= faniso(idofn)*vvel*wvel
           shmet(3,1,idofn)= shmet(1,3,idofn)
           shmet(3,2,idofn)= shmet(2,3,idofn)
        end if

     end do

  else if (itask == 2) then

     if (kfl_shock_nsa(1) == 0) return

     if (shock_nsa < zesho) then
        kapsh = 0
        return
     end if

     do idofn= 1,ndofn_nsa
        jdofn= idofn
        if (kapsh(idofn) > 0) then
           do jdime= 1,ndime
              do idime= 1,ndime
                 shote(idofn)= shote(idofn) + (cartd(idime)*shmet(jdime,idime,jdofn))*graun(idofn,jdime)
              end do
           end do
! estoy harto de este menos!!!
!!!           shote(idofn) = -shote(idofn)
        end if
     end do

  end if

end subroutine nsa_shocapxy
