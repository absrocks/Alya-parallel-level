subroutine got_elmrhs(&
     pgaus,ndofn,gppor,gpcdr,gpvdr,gpgcd,gpdiv,gpugu,gpvel,gprhs)
  !----------------------------------------------------------------------
  !****f* Gotita/got_elmrhs
  ! NAME 
  !    got_elmrhs
  ! DESCRIPTION
  !    This routine computes the right-hand side of the momentum 
  !    and continuity equations, including:
  !    Time derivative ....... alpha_l*u^{n-1}/(theta*dt)
  ! INPUT
  ! OUTPUT
  !    GPRHS ................. Right hand side
  ! USES
  ! USED BY
  !    got_elmope 
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  use def_gotita, only     :  kfl_timei_got,dtinv_got,ndofn_got,&
       &                      penal_got,kfl_forme_got,kfl_exacs_got,&
       &                      kfl_probl_got,kfl_coupl_got,kfl_timem_got,&
       &                      kfl_timec_got,almax_got,kfl_linea_got
  implicit none
  integer(ip), intent(in)  :: pgaus,ndofn
  real(rp),    intent(in)  :: gppor(pgaus)
  real(rp),    intent(in)  :: gpcdr(pgaus,2),gpvdr(ndime,pgaus,2)
  real(rp),    intent(in)  :: gpgcd(ndime,pgaus),gpdiv(pgaus)
  real(rp),    intent(in)  :: gpugu(ndime,pgaus),gpvel(ndime,pgaus)
  real(rp),    intent(out) :: gprhs(ndofn,pgaus)
  integer(ip)              :: igaus,idime
  real(rp)                 :: fact1

  !----------------------------------------------------------------------
  !
  ! BOTH EQUATIONS
  ! 
  !----------------------------------------------------------------------

  if(kfl_probl_got==1) then
     !
     ! Initialization
     !
     if(kfl_exacs_got==0) then
        do igaus=1,pgaus
           do idime=1,ndofn_got(3)
              gprhs(idime,igaus)=0.0_rp
           end do
        end do
     end if
     !
     ! Time derivatives 
     !
     if(kfl_timem_got==1) then
        if(kfl_coupl_got==1) then
           !
           ! Momentum equation (coupled eqns): (alpha+penal)/(theta*dt)*u^{n-1}
           !
           do igaus=1,pgaus
              fact1=dtinv_got*(gpcdr(igaus,1)+penal_got)
              do idime=1,ndime
                 gprhs(idime,igaus)=gprhs(idime,igaus)&
                      +fact1*gpvdr(idime,igaus,2)
              end do
           end do
           if(kfl_linea_got==2) then
              do igaus=1,pgaus
                 fact1=gpcdr(igaus,1)+penal_got
                 do idime=1,ndime
                    gprhs(idime,igaus)=gprhs(idime,igaus)&
                         +fact1*gpugu(idime,igaus)
                 end do
              end do
           end if
        else
           !
           ! Momentum equation (uncoupled eqns): 1/(theta*dt)*u^{n-1} 
           !
           do igaus=1,pgaus
              fact1=dtinv_got*almax_got
              do idime=1,ndime
                 gprhs(idime,igaus)=gprhs(idime,igaus)&
                      +fact1*gpvdr(idime,igaus,2)
              end do
           end do
           if(kfl_linea_got==2) then
              do igaus=1,pgaus
                 do idime=1,ndime
                    gprhs(idime,igaus)=gprhs(idime,igaus)&
                         +gpugu(idime,igaus)*almax_got
                 end do
              end do
           end if
        end if
     end if

     if(kfl_timec_got==1) then
        !
        ! Continuity equation: 1/(theta*dt)*alpha^{n-1}
        !
        do igaus=1,pgaus
           gprhs(ndofn_got(3),igaus)=gprhs(ndofn_got(3),igaus)&
                +dtinv_got*gpcdr(igaus,2)
        end do
     end if

     if(kfl_coupl_got==1) then
        !
        ! Linearization momentum (coupled eqns): alpha*[(u.grad)u+sig*u]+eps*sig*ua
        !
        do igaus=1,pgaus
           fact1=gppor(igaus)*penal_got
           do idime=1,ndime
              !gprhs(idime,igaus)=gprhs(idime,igaus)&
              !     +gpcdr(igaus,1)&
              !     *(gpugu(idime,igaus)+gppor(igaus)*gpvdr(idime,igaus,1))&
              !     +fact1*gpvel(idime,igaus)
           end do
        end do
        do igaus=1,pgaus
           fact1=gpcdr(igaus,1)+penal_got
           !fact1=penal_got
           do idime=1,ndime
              gprhs(idime,igaus)=gprhs(idime,igaus)&
                   +fact1*gppor(igaus)*gpvel(idime,igaus)
           end do
        end do
     else
        !
        ! RHS momentum (uncoupled eqns): sig*ua
        !      
        do igaus=1,pgaus
           fact1=almax_got*gppor(igaus)
           do idime=1,ndime
              gprhs(idime,igaus)=gprhs(idime,igaus)&
                   +fact1*gpvel(idime,igaus)
           end do
        end do
     end if
     !
     ! Linearization continuity: div(alpha*u)
     !
     if(kfl_coupl_got==1) then
        do igaus=1,pgaus
           !gprhs(ndofn_got(3),igaus)=gprhs(ndofn_got(3),igaus)&    ! alpha*div(u)
           !     +gpcdr(igaus,1)*gpdiv(igaus)
           !do idime=1,ndime
           !   gprhs(ndofn_got(3),igaus)=gprhs(ndofn_got(3),igaus)& ! grad(alpha).u
           !        +gpgcd(idime,igaus)*gpvdr(idime,igaus,1)
           !end do
        end do
     end if
     !
     ! Conservative form: u*(continuity residual) 
     !
     if(kfl_forme_got==2.and.kfl_coupl_got==1) then
        do igaus=1,pgaus 
           do idime=1,ndime
              gprhs(idime,igaus)=gprhs(idime,igaus)&
                   +gprhs(ndofn_got(3),igaus)*gpvdr(idime,igaus,1)
           end do
        end do
     end if

     !----------------------------------------------------------------------
     !
     ! MOMENTUM EQUATION
     ! 
     !----------------------------------------------------------------------

  else if(kfl_probl_got==2) then
     !
     ! Initialization
     !
     if(kfl_exacs_got==0) then
        do igaus=1,pgaus
           do idime=1,ndofn
              gprhs(idime,igaus)=0.0_rp
           end do
        end do
     end if
     !
     ! Time derivative: 1/(theta*dt)*u^{n-1}
     !
     if(kfl_timem_got==1) then
        do igaus=1,pgaus
           do idime=1,ndime
              gprhs(idime,igaus)=gprhs(idime,igaus)&
                   +dtinv_got*gpvdr(idime,igaus,2)
           end do
        end do
     end if
     !
     ! Source term: sig*ua
     !
     do igaus=1,pgaus
        do idime=1,ndime
           gprhs(idime,igaus)=gprhs(idime,igaus)&
                +gppor(igaus)*gpvel(idime,igaus)
        end do
     end do
     !
     ! Newton-Raphson: (u.grad)u
     !
     if(kfl_linea_got==2) then
        do igaus=1,pgaus
           do idime=1,ndime
              gprhs(idime,igaus)=gprhs(idime,igaus)&
                   +gpugu(idime,igaus)
           end do
        end do
     end if

     !----------------------------------------------------------------------
     !
     ! CONTINUITY EQUATION
     ! 
     !----------------------------------------------------------------------

  else if(kfl_probl_got==3) then

     if(kfl_timec_got==1) then
        !
        ! Time derivatives: 1/(theta*dt)*alpha^{n-1}
        !
        do igaus=1,pgaus
           gprhs(1,igaus)=dtinv_got*gpcdr(igaus,2)
        end do
     else
        do igaus=1,pgaus
           gprhs(1,igaus)=0.0_rp
        end do
     end if

  end if

end subroutine got_elmrhs
