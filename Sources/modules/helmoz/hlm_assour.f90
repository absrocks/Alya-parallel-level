!subroutine hlm_assour(pnode,gpmas,gpmav,gpels,gpelv,gpvol,shapf,elrhs)
  !------------------------------------------------------------------------
  !****f* Helmoz/hlm_assour
  ! NAME 
  !    hlm_assour
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    hlm_elmope
  !***
  !------------------------------------------------------------------------
  !use def_kintyp, only     :  ip,rp
  !use def_domain, only     :  ndime
  !use def_helmoz, only     :  nequs_hlm,kfl_model_hlm
  !use def_master, only     :  kfl_paral

  !implicit none

  !integer(ip), intent(in)  :: pnode
  !complex(rp), intent(in)  :: gpmas
  !complex(rp), intent(in)  :: gpmav(ndime)
  !complex(rp), intent(in)  :: gpels
  !complex(rp), intent(in)  :: gpelv(ndime)
  !real(rp),    intent(in)  :: gpvol,shapf(pnode)
  !complex(rp), intent(out) :: elrhs(pnode,nequs_hlm)
  !integer(ip)              :: inode,idime
  !integer(ip)              :: ndof2,ndof3,ndof4
  !real(rp)                 :: xfact

  !ndof2 = ndime + 1
  !ndof3 = ndime + 1
  !ndof4 = 2 * ndime + 2

  !if( kfl_model_hlm == 1 ) then
     !
     ! Pure magnetic
     !
     !do inode = 1,pnode            
        !xfact = gpvol * shapf(inode)
        !elrhs(inode,ndof2) = elrhs(inode,ndof2) + gpmas * xfact
        !do idime = 1,ndime
           !elrhs(inode,idime) = elrhs(inode,idime) + gpmav(idime) * xfact
        !end do
     !end do
  !else if( kfl_model_hlm == 2 ) then
     !
     ! Pure electric
     !
     !do inode = 1,pnode            
        !xfact = gpvol * shapf(inode)
        !elrhs(inode,ndof2) = elrhs(inode,ndof2) + gpels * xfact
        !do idime = 1,ndime
           !elrhs(inode,idime) = elrhs(inode,idime) + gpelv(idime) * xfact
        !end do
     !end do

  !else if( kfl_model_hlm == 3 ) then
     !
     ! Magnetic-electric problem
     !     
     !do inode = 1,pnode            
        !xfact = gpvol * shapf(inode)
        !elrhs(inode,ndof2) = elrhs(inode,ndof2) + gpmas * xfact
        !elrhs(inode,ndof4) = elrhs(inode,ndof4) + gpels * xfact
        !do idime = 1,ndime
           !elrhs(inode,idime)       = elrhs(inode,idime)       + gpmav(idime) * xfact
           !elrhs(inode,ndof3+idime) = elrhs(inode,ndof3+idime) + gpelv(idime) * xfact
        !end do
     !end do

  !end if

!end subroutine hlm_assour
