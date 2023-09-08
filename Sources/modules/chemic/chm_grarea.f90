subroutine chm_grarea()
  !------------------------------------------------------------------------
  !****f* partis/chm_elmadr
  ! NAME 
  !    chm_elmadr
  ! DESCRIPTION
  !    Elemental operations
  ! USES
  ! USED BY
  !    chm_matrix
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_memchk
  implicit none
  integer(ip)          :: ireac,nlelp,ispec,jspec,mepoi_chm
  integer(4)           :: istat
  integer(ip), pointer :: nepoi_chm(:)
  logical(lg), pointer :: lepoi_chm(:)

  allocate(nepoi_chm(nspec_chm),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'NEPOI_CHM','connpo',nepoi_chm)
  allocate(lepoi_chm(nspec_chm),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'LEPOI_CHM','connpo',lepoi_chm)
  !
  ! NEPOI_CHM(ISPEC) = number of reaction in which ISPEC appears
  !
  do ireac = 1,nreac_chm
     do jspec = 1,nspec_chm
        lepoi_chm(jspec) = .true.
     end do
     do ispec = 1,size(lreac_chm(ireac)%l)
        jspec = abs(lreac_chm(ireac)%l(ispec))
        if( jspec /= 0 ) then
           if( lepoi_chm(jspec) ) then
              nepoi_chm(jspec) = nepoi_chm(jspec) + 1
              lepoi_chm(jspec) = .false.
           end if
        end if
     end do
  end do
  !
  ! Allocate memory for IAREA and compute it
  !
  allocate(iarea_chm(nspec_chm+1),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'IAREA_CHM','connpo',iarea_chm)
  iarea_chm(1) = 1
  do ispec = 1,nspec_chm
     iarea_chm(ispec+1) = iarea_chm(ispec) + nepoi_chm(ispec)
  end do
  !
  ! Allocate memory for LELPO and construct the list
  !
  nlelp = iarea_chm(nspec_chm+1)
  allocate(jarea_chm(nlelp),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'JAREA_CHM','connpo',jarea_chm)

  do ireac = 1,nreac_chm
     do jspec = 1,nspec_chm
        lepoi_chm(jspec) = .true.
     end do
     do ispec = 1,size(lreac_chm(ireac)%l)
        jspec = abs(lreac_chm(ireac)%l(ispec))
        if( jspec /= 0 ) then
           if( lepoi_chm(jspec) ) then
              lepoi_chm(jspec)            = .false.
              jarea_chm(iarea_chm(jspec)) = ireac
              iarea_chm(jspec)            = iarea_chm(jspec)+1
           end if
        end if
     end do
  end do
  !
  ! Recompute IAREA_CHM and maximum number of element neighbors MEPOI_CHM
  !
  iarea_chm(1) =  1
  mepoi_chm    = -1
  do ispec = 1,nspec_chm
     iarea_chm(ispec+1) = iarea_chm(ispec) + nepoi_chm(ispec)
     mepoi_chm = max(mepoi_chm,nepoi_chm(ispec))
  end do

  call memchk(two,istat,mem_modul(1:2,modul),'LEPOI_CHM','chm_membcs',lepoi_chm)
  deallocate(lepoi_chm,stat=istat)
  if(istat/=0) call memerr(two,'LEPOI_CHM','nsi_membcs',0_ip)

  call memchk(two,istat,mem_modul(1:2,modul),'NEPOI_CHM','chm_membcs',nepoi_chm)
  deallocate(nepoi_chm,stat=istat)
  if(istat/=0) call memerr(two,'NEPOI_CHM','nsi_membcs',0_ip)

end subroutine chm_grarea
