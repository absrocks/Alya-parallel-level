subroutine qua_elmres(&
     pnode,pgaus,plapl,gpden,gpgrd,gprea,gpadv,gpcon,&
     gpcod,gpsha,gpcar,gplap,rtemp)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_elmres
  ! NAME
  !   qua_elmres
  ! DESCRIPTION
  !    Compute residual of the Schrodinger equation
  ! OUTPUT 
  !    GPRHS
  ! USES
  ! USED BY
  !    qua_elmope 
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  mnode,ndime,kfl_naxis
  use def_quanty

  implicit none
  integer(ip), intent(in)  :: pnode,pgaus,plapl
  real(rp),    intent(in)  :: gpden(pgaus),gpgrd(ndime,pgaus)
  real(rp),    intent(in)  :: gprea(pgaus),gpadv(pnode,pgaus)
  real(rp),    intent(in)  :: gpcon(pgaus),gpcod(ndime,pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gplap(pnode,pgaus)
  real(rp),    intent(out) :: rtemp(pnode,pgaus)
  integer(ip)              :: igaus,inode,idime
  real(rp)                 :: fact1,factt
  !
  ! Initialization: s*Nj
  !
  do igaus=1,pgaus                
     do inode=1,pnode
        rtemp(inode,igaus)=gprea(igaus)*gpsha(inode,igaus)
     end do
  end do
  !
  ! Time integration: rho*cp/(dt*theta)*Nj
  !
 ! if(kfl_timei_qua==1) then
 !    if(kfl_tisch_qua==1) then
 !       factt=dtinv_qua                 ! Trapezoidal rule
 !    else
 !       factt=dtinv_qua*pabdf_qua(1)    ! BDF scheme
 !    end if
 !    do igaus=1,pgaus
 !       fact1=factt*gpden(igaus)
 !       do inode=1,pnode
 !          rtemp(inode,igaus)=rtemp(inode,igaus)+fact1*gpsha(inode,igaus)
 !       end do
 !    end do
 ! end if
  !
  ! Axisymmetric flow: -k/r*dNj/dr
  !
  !if(kfl_naxis==1) then
  !   do igaus=1,pgaus
  !      fact1=-gpcon(igaus)/gpcod(1,igaus)
  !      do inode=1,pnode
  !         rtemp(inode,igaus)=rtemp(inode,igaus)&
  !              +fact1*gpcar(1,inode,igaus)
  !      end do
  !   end do     
  !end if
  !
  ! Diffusion term
  !
!  if(kfl_condu_qua==1.and.kfl_taust_qua/=0) then

!     if(plapl==1) then
        !
        ! Laplacian: -k*lapl(Nj) 
        !
        do igaus=1,pgaus
           do inode=1,pnode
              rtemp(inode,igaus)=rtemp(inode,igaus)-gplap(inode,igaus)
           end do
        end do
!     end if


!  end if

end subroutine qua_elmres
