subroutine got_exacso(itask,gpcod,gpvel,gprhs)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_exacso
  ! NAME 
  !    got_exacso
  ! DESCRIPTION
  !    This routine:
  !    ITASK=1 ... Returns exact solution
  !    ITASK=2 ... Returns force vector 
  ! USES
  ! USED BY
  !    got_elmexa
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_gotita, only       :  kfl_exacs_got,exvdr_got,excdr_got,&
       &                        exgvd_got,exgcd_got,kfact_got,&
       &                        deair_got,ddrop_got,veair_got,&
       &                        muair_got,kfl_coupl_got
  implicit none
  integer(ip), intent(in)    :: itask
  real(rp),    intent(in)    :: gpcod(*),gpvel(*)
  real(rp),    intent(inout) :: gprhs(*)
  real(rp)                   :: u,v,a,dudx,dudy,dvdx,dvdy,dadx,dady
  real(rp)                   :: x,y,ua,va,uua(3),unorm,CDReD,ReD,gppor
  !
  ! Initializations
  ! 
  x    = gpcod(1)
  y    = gpcod(2)
  a    = 0.0_rp
  u    = 0.0_rp
  v    = 0.0_rp
  dudx = 0.0_rp 
  dudy = 0.0_rp 
  dvdx = 0.0_rp 
  dvdy = 0.0_rp 
  dadx = 0.0_rp 
  dady = 0.0_rp

  if (kfl_exacs_got==1) then
     !
     ! Constant solution
     !
     u    = 1.0_rp
     v    = 2.0_rp
     a    = 3.0_rp
     dudx = 0.0_rp 
     dudy = 0.0_rp 
     dvdx = 0.0_rp 
     dvdy = 0.0_rp 
     dadx = 0.0_rp 
     dady = 0.0_rp

  else if (kfl_exacs_got==2) then
     !
     ! Linear solution
     !
     u    = x+2.0_rp*y
     v    = 0.0_rp
     a    = 1.0_rp+x+0.5_rp*y
     dudx = 1.0_rp
     dudy = 2.0_rp 
     dvdx = 0.0_rp
     dvdy = 0.0_rp
     dadx = 1.0_rp 
     dady = 0.5_rp
  end if

  if(itask==1) then
     !
     ! Exact solution and gradients
     !
     exvdr_got(1)   = u
     exvdr_got(2)   = v
     excdr_got      = a
     exgvd_got(1,1) = dudx
     exgvd_got(2,1) = dudy
     exgvd_got(1,2) = dvdx
     exgvd_got(2,2) = dvdy
     exgcd_got(1)   = dadx
     exgcd_got(2)   = dady

  else if(itask==2) then
     !
     ! Force term GPRHS
     !
     ua       = gpvel(1)
     va       = gpvel(2)
     uua(1)   = u-ua
     uua(2)   = v-va
     call vecnor(uua,ndime,unorm,2_ip)
     ReD      = deair_got*ddrop_got*veair_got*unorm/muair_got
     call got_dragco(ReD,CDReD)
     gppor    = CDReD/(24.0_rp*kfact_got)
     if(kfl_coupl_got==1) then
        gprhs(1) = a*u*dudx + a*v*dudy + gppor*a*(u-ua)
        gprhs(2) = a*u*dvdx + a*v*dvdy + gppor*a*(v-va)
     else
        gprhs(1) = u*dudx + v*dudy + gppor*(u-ua)
        gprhs(2) = u*dvdx + v*dvdy + gppor*(v-va)
     end if
     gprhs(3) = dadx*u   + dady*v   + a*(dudx+dvdy)

  end if

end subroutine got_exacso
