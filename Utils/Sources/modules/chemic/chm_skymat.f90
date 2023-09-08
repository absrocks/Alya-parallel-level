subroutine chm_skymat()
  !------------------------------------------------------------------------
  !****f* partis/chm_skymat
  ! NAME 
  !    chm_skymat
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_memchk
  implicit none
  integer(ip)          :: iodes,ispec,jodes,iarea,ireac,nspe1,iode1
  integer(4)           :: istat
  logical(lg), pointer :: touch(:)

  allocate(touch(nodes_chm),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'TOUCH','chm_skymat',touch)

  allocate(iskyl_chm(nodes_chm+1),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'ISKYL_CHM','chm_skymat',iskyl_chm)
  allocate(idiag_chm(nodes_chm),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'IDIAG_CHM','chm_skymat',idiag_chm)

  nskyl_chm    = 1
  iskyl_chm(1) = 1

  do iodes = 1,nodes_chm

     ispec = iodes+nclas_chm
     do jodes = 1,nodes_chm
        touch(jodes)=.false.
     end do

     do iarea =  iarea_chm(ispec),iarea_chm(ispec+1)-1
        ireac =  jarea_chm(iarea)
        nspe1 =  size(lreac_chm(ireac)%l)
        do ispec = 1,nspe1
           jodes = lreac_chm(ireac)%l(ispec)-nclas_chm
           if(jodes>0)&
                touch(jodes)=.true.
        end do
     end do

     iode1 =  1e6
     do jodes = 1,nodes_chm
        if(touch(jodes)) then
           if(jodes<iode1) iode1=jodes
        end if
     end do
     nskyl_chm          = nskyl_chm + (iodes - iode1)*2 + 1
     iskyl_chm(iodes+1) = nskyl_chm
     idiag_chm(iodes)   = nskyl_chm - (iodes - iode1)   - 1

  end do
  nskyl_chm = nskyl_chm -1

  allocate(amatr_chm(nskyl_chm),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'AMATR_CHM','chm_skymat',amatr_chm)
  allocate(rhsid_chm(nodes_chm),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'RHSID_CHM','chm_skymat',rhsid_chm)
  
  call memchk(two,istat,mem_modul(1:2,modul),'TOUCH','chm_skymat',touch)
  deallocate(touch,stat=istat)
  if(istat/=0) call memerr(two,'TOUCH_CHM','chm_skymat',0_ip)

end subroutine chm_skymat
       
