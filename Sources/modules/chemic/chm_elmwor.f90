subroutine chm_elmwor(pnode,iclas,elcon,elco1)
  use def_kintyp, only     :  ip,rp
  use def_chemic, only     :  nclas_chm,kfl_timei_chm,&
       &                      kfl_tisch_chm,kfl_tiacc_chm,&
       &                      nspec_chm,ncomp_chm
  implicit none
  integer(ip), intent(in)  :: pnode,iclas
  real(rp),    intent(in)  :: elcon(pnode,nspec_chm,ncomp_chm)
  real(rp),    intent(out) :: elco1(pnode,ncomp_chm)
  integer(ip)              :: inode,itime

  do inode=1,pnode
     elco1(inode,1)=elcon(inode,iclas,1)
  end do

  if(kfl_timei_chm==1) then
     do inode=1,pnode
        elco1(inode,2)=elcon(inode,iclas,2)
     end do
     if(kfl_tisch_chm==2) then
        do itime=3,1+kfl_tiacc_chm
           do inode=1,pnode
              elco1(inode,itime) = elcon(inode,iclas,itime)
           end do
        end do
     end if
  end if

end subroutine chm_elmwor
