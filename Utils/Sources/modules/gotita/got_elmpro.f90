subroutine got_elmpro(&
     pnode,pgaus,igaui,igauf,eldif,gpsha,gpvdr,gpvel,&
     gpvno,gpcdr,chale,gppor,gpdif)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmpro
  ! NAME 
  !    got_elmpro
  ! DESCRIPTION
  !    Compute the porosity-like and diffusion terms:
  !
  !                  CD*ReD             rho_a*d*||u_a-u|| 
  !    GPPOR = sig = ------*alpha, ReD= -----------------
  !                   24 K                    mu_a
  !           
  !                                     +-      +-    alpha       -+ C2-+
  !                                     | - log | ------------ + 1 |    |         
  !                                     +-      +- C1*alpha_inf   -+   -+
  !    GPDIF = k = alpha_inf*u_inf*C3*exp
  !                            
  !                         
  !    k is added so that the diffusion term appears when alpha tends
  !    to zero, that is when the momentum equation starts to be singular:
  ! 
  !    4*k/h^2 + alpha*( 2*u/h + sig )
  !    (C*alpha_max+alpha)*( 2*u/h + sig )
  !  
  ! USES  
  !
  ! USED BY
  !    got_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  use def_gotita, only     :  densi_got,deair_got,muair_got,veair_got,&
       &                      ddrop_got,kfact_got,kfl_diffu_got,&
       &                      diffu_got,almax_got,vemax_got,dimin_got,&
       &                      dimax_got
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus,igaui,igauf
  real(rp),    intent(in)  :: eldif(pnode),gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpvdr(ndime,pgaus),gpvel(ndime,pgaus)
  real(rp),    intent(in)  :: gpvno(pgaus),gpcdr(pgaus),chale(2)
  real(rp),    intent(out) :: gppor(pgaus),gpdif(pgaus)
  real(rp)                 :: vdiff(3),vnorm,ReD,CDReD,gpven,fact1
  integer(ip)              :: igaus,idime,inode
  !
  ! Porosity sig
  !
  do igaus=igaui,igauf
     !
     ! ReD= rho_a*d*||u_a-u||/mu_a
     !
     do idime=1,ndime
        vdiff(idime)=gpvel(idime,igaus)-gpvdr(idime,igaus)
     end do
     call vecnor(vdiff,ndime,vnorm,2_ip)
     Red=deair_got*ddrop_got*vnorm/muair_got
     call got_dragco(ReD,CDReD)
     !
     ! Porosity: GPPOR=sig
     !
     gppor(igaus)=CDReD/(24.0_rp*kfact_got)
  end do
  !
  ! Diffusion: GPDIF=k
  !
  if(kfl_diffu_got==0) then

     do igaus=igaui,igauf
        gpdif(igaus)=0.0_rp
     end do

  else if(kfl_diffu_got==1) then
     !
     ! GPDIF computed at Gauss points
     !
     do igaus=igaui,igauf
        call got_diffun(gpcdr(igaus),chale,gpdif(igaus))
        if(gpdif(igaus)>dimax_got) dimax_got=gpdif(igaus)
        if(gpdif(igaus)<dimin_got) dimin_got=gpdif(igaus)
     end do

  else if(kfl_diffu_got==2) then
     !
     ! Smoothed GPDIF: computed in subroutine got_diffus()
     !
     do igaus=igaui,igauf
        gpdif(igaus)=0.0_rp
        do inode=1,pnode
           gpdif(igaus)=gpdif(igaus)+gpsha(inode,igaus)*eldif(inode)
        end do
     end do
  end if

end subroutine got_elmpro
