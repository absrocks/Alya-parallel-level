subroutine chm_boumat(&
     iclas,pnode,pnodb,pgaub,pgaus,lboel,gbsha,&
     gbsur,shaga,gbdif,elmat,elrhs)
  !------------------------------------------------------------------------
  !****f* temper/chm_boumat
  ! NAME 
  !    chm_boumat
  ! DESCRIPTION
  !    Assemble boundary contribution
  ! USES
  ! USED BY
  !    chm_bouope
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_chemic, only       :  nclas_chm,panat_chm,boltz_chm,&
       &                        equil_chm
  implicit none
  integer(ip), intent(in)    :: iclas,pnode,pnodb,pgaub,pgaus
  integer(ip), intent(in)    :: lboel(pnodb)
  real(rp),    intent(in)    :: gbsha(pnodb,pgaub)
  real(rp),    intent(in)    :: gbsur(pgaub)
  real(rp),    intent(in)    :: shaga(pgaus,pnode)
  real(rp),    intent(in)    :: gbdif(pgaub)
  real(rp),    intent(inout) :: elmat(pnode,pnode)
  real(rp),    intent(inout) :: elrhs(pnode)
  integer(ip)                :: inodb,jnodb,inode,jnode,igaub
  real(rp)                   :: xmuit,xmrhs,xmmat,ovkT
  real(rp)                   :: C0eq,Eform,Ceq,para1,para2,para3
  real(rp)                   :: k,lambda,T
  !
  ! Initialization
  !
  call chm_usrtem(T)
  ovkT = 1.0_rp/(boltz_chm*T) 

  do inode=1,pnode
     elrhs(inode)=0.0_rp
     do jnode=1,pnode
        elmat(jnode,inode)=0.0_rp
     end do
  end do

  do igaub=1,pgaub
     !
     ! Interpolate extrapolated diffusion
     !
     !k=0.0_rp
     !do igaus=1,pgaus
     !   do inodb=1,pnodb
     !      tmatr=shaga(igaus,lboel(inodb))*gbsha(inodb,igaub)
     !      k=k+gpdif(igaus)*tmatr
     !   end do
     !end do
     !
     ! Robin condition:
     !
     ! k*grad(C).n = P3 + P1 * ( C - P2 )
     !             = k/lambda * [ C -  C0eq * exp(-E_form/kT) ] + q
     !
     k      = gbdif(igaub)                       ! k
     C0eq   = equil_chm(1,iclas)                 ! C0eq
     Eform  = equil_chm(2,iclas)                 ! E_form: formation energy
     Ceq    = C0eq*exp(-Eform*ovkT)              ! Ceq = C0eq * exp(-E_form/kT)
     lambda = panat_chm(1,iclas)                 ! lambda

     para1  = k/lambda                           ! P1  = k/lambda
     para2  = Ceq                                ! P2  = Ceq
     para3  = panat_chm(2,iclas)                 ! P3  = q
     !
     ! (...) + ( -P1*C, v ) = (...) + ( P3 - P1*P2 , v )
     !
     xmrhs  = (para3+para1*para2)*gbsur(igaub)
     xmmat  = (para1)*gbsur(igaub)

     do inodb=1,pnodb  
        inode=lboel(inodb)
        elrhs(inode)=elrhs(inode)+gbsha(inodb,igaub)*xmrhs
        xmuit=xmmat*gbsha(inodb,igaub)
        do jnodb=1,pnodb
           jnode=lboel(jnodb)
           elmat(jnode,inode)=elmat(jnode,inode)+xmuit*gbsha(jnodb,igaub)
        end do
     end do

  end do

end subroutine chm_boumat
