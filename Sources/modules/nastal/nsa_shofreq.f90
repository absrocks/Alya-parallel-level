subroutine nsa_shofreq(kapsh,graun,resid,chale,xvelo,velmo,zesho,mfreq)
!-----------------------------------------------------------------------------------------------
!
! 
!
!
!
!
!-----------------------------------------------------------------------------------------------
  use      def_kintyp
  use      def_domain, only:  ndime
  use      def_nastal, only:  kfl_shock_nsa,ndofn_nsa,mfreq_nsa,frmax_nsa
  implicit none
  integer(ip)           :: kfcdr,inode,ipoin,idome,igaus,itidi,idime,jdime,idofn,itime

  real(rp), intent(in)  :: graun(ndofn_nsa,ndime),resid(ndofn_nsa),chale(2),xvelo(ndime)
  real(rp), intent(in)  :: zesho,velmo

  integer(ip)           :: kapsh(ndofn_nsa),napsh,mfreq
  real(rp)              :: &
       grau2(ndofn_nsa),grmod(ndofn_nsa),remod(ndofn_nsa),&
       ficve(ndofn_nsa),shpec,difau,shfau,zesh2,hconv
  real(rp)              :: epsdc,epssu,epssl,epsde,fiso(ndofn_nsa),faniso(ndofn_nsa),uvel,vvel,wvel,xvel2

     hconv = chale(1)

     mfreq = mfreq_nsa

     kapsh = 1             ! default: apply shock capturing
     xvel2 = velmo*velmo   ! velocity squared
     zesh2 = zesho*zesho

!     if (xvel2 < zesho) kapsh = 0            ! velocity is the first threshold  
     
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
        if (ficve(idofn) < zesho) kapsh(idofn)= 0
        napsh = napsh + kapsh(idofn)
     end do

     if (napsh .ne. 0) mfreq = frmax_nsa

end subroutine nsa_shofreq
