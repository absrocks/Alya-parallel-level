subroutine elmchl(&
     tragl,hleng,elcod,elvel,chave,chale,pnode,porde,&
     hnatu,kfl_advec,kfl_ellen)
  !-----------------------------------------------------------------------
  !****f* Domain/elmchl
  ! NAME
  !   elmchl
  ! DESCRIPTION
  !   This routine computes the characteristic element lengths CHALE 
  !   according to a given strategy. CHALE is divided by two for
  !   quadratic elements:
  !   KFL_ELLEN = 0 ... CHALE(1) = Minimum element length
  !                 ... CHALE(2) = Minimum element length
  !   KFL_ELLEN = 1 ... CHALE(1) = Maximum element length
  !                 ... CHALE(2) = Maximum element length
  !   KFL_ELLEN = 2 ... CHALE(1) = Average element length
  !                 ... CHALE(2) = Average element length
  !   KFL_ELLEN = 3 ... IF KFL_ADVEC = 1:
  !                     CHALE(1) = Flow direction
  !                     CHALE(2) = Flow direction
  !                     ELSE IF KFL_ADVEC =0:
  !                     CHALE(1) = Minimum element length
  !                     CHALE(2) = Minimum element length
  !   KFL_ELLEN = 4 ... CHALE(1) = Approx. diameter=sqrt(hmin*hmax)
  !                 ... CHALE(2) = Approx. diameter=sqrt(hmin*hmax)
  !   KFL_ELLEN = 5 ... CHALE(1) = Length in flow direction
  !                 ... CHALE(2) = Minimum element kength
  ! OUTPUT
  !   CHALE
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  implicit none
  integer(ip), intent(in)  :: pnode,porde,kfl_advec,kfl_ellen
  real(rp),    intent(in)  :: hnatu
  real(rp),    intent(out) :: chale(2)
  real(rp),    intent(in)  :: tragl(ndime,ndime),hleng(ndime)
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  real(rp),    intent(in)  :: elvel(ndime,pnode)
  real(rp),    intent(out) :: chave(ndime,2)
  integer(ip)              :: idime,inode
  real(rp)                 :: elno1,elno2

  if(kfl_ellen==0) then 
     !
     ! Minimum element length
     !
     chale(1)=hleng(ndime) 
     chale(2)=chale(1)

  else if(kfl_ellen==1) then   
     !
     ! Maximum element length
     !     
     chale(1)=hleng(1) 
     chale(2)=chale(1)

  else if(kfl_ellen==2) then 
     !
     ! Average length
     !
     chale(1)=0.0_rp
     do idime=1,ndime
        chale(1)=chale(1)+hleng(idime)
     end do
     chale(1)=chale(1)/real(ndime,rp) 
     chale(2)=chale(1)

  else if(kfl_ellen==3) then 
     !
     ! Length in flow direction
     !
     if(kfl_advec/=0) then 
        !
        ! Characteristic element velocity (average)
        !
        chave=0.0_rp
        do idime=1,ndime
           do inode=1,pnode
              chave(idime,1)=chave(idime,1)+elvel(idime,inode)
           end do
           chave(idime,1)=chave(idime,1)/real(pnode,rp)
        end do
        !
        ! Characteristic element length u^l = J^(-t) u^g
        !
        call mbvab1(chave(1,2),tragl,chave(1,1),ndime,ndime,elno2,elno1)
        if(elno2>1.0e-16_rp.and.elno1>1.0e-16_rp) then
           chale(1)=hnatu*elno1/elno2
        else
           chale(1)=hleng(ndime)
        end if
        chale(2)=chale(1)
        chale(2)=hleng(ndime)
        if (ndime ==3 ) then
           chale(2)=(hleng(ndime)*hleng(2)*hleng(1))**(1.0_rp/3.0_rp)
        else if (ndime==2) then
           chale(2)=sqrt(hleng(2)*hleng(1))
        end if
     else
        chale(1)=hleng(ndime)       
        chale(2)=chale(1)
     end if

  else if(kfl_ellen==4) then 
     !
     ! sqrt(hmin*hmax)
     !
     chale(1)=sqrt(hleng(1)*hleng(ndime))     
     chale(2)=chale(1)

  else if(kfl_ellen==5) then 
     !
     ! Along velocity direction
     !
     call velchl(pnode,elcod,elvel,chale,hleng)

  else if(kfl_ellen==6) then 
     !
     ! Mixed element length - hmin for tau1, hmax for tau2 - here we only obtain the values for tau1 - tau2 directly in nsi_elmsgs
     !
     chale(1)=hleng(ndime) 
     chale(2)=chale(1)

  end if
  !
  ! Divide h by 2 for quadratic elements and 3 for cubic elements
  !
  chale(1) = chale(1)/real(porde,rp)
  chale(2) = chale(2)/real(porde,rp)

end subroutine elmchl

