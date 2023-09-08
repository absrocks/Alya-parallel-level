subroutine got_memall()
  !-----------------------------------------------------------------------
  !****f* Gotita/got_memall
  ! NAME 
  !    got_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    NS equations
  ! USES
  !    memchk
  !    mediso
  ! USED BY
  !    got_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver
  use def_gotita
  use mod_memchk
  implicit none
  integer(ip) :: lodof,ncsgs,ielem,pgaus,pelty,nlpns,ivari,nvari
  integer(4)  :: istat
  !
  ! Problem unknowns VDROP, CDROP, VESGS and solver initialization
  !
  if( INOTMASTER ) then
     !
     ! Allocate memory for droplet velocity and droplet concentration: VDROP, PRESS
     ! vdrop(:,:,1)=u^{n,i}
     ! vdrop(:,:,2)=u^{n,i-1}
     ! vdrop(:,:,3)=u^{n-1}
     ! vdrop(:,:,4)=u^{n-2}
     ! vdrop(:,:,5)=u^{n-3}, etc.
     !
     if(kfl_probl_got<=2) then
        allocate(vdrop(ndime,npoin,ncomp_got),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'VELOC','got_memall',vdrop)
     end if
     if(kfl_probl_got/=2) then
        allocate(cdrop(npoin,ncomp_got),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'PRESS','got_memall',cdrop)
     end if
     !
     ! Subgrid scale: VESGS
     !
     if(kfl_sgsco_got==1.or.kfl_sgsti_got==1) then
        allocate(vesgs(nelem),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'VESGS','got_memall',vesgs)
        ncsgs=min(2,2*kfl_sgsti_got+kfl_sgsco_got)
        do ielem=1,nelem
           pelty=ltype(ielem)
           pgaus=ngaus(pelty)
           allocate(vesgs(ielem)%a(ndime,pgaus,ncsgs),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'VESGS','got_memall',vesgs(ielem)%a)
        end do
        if(kfl_sgsco_got==1) then
           allocate(itsta_got(misgs_got),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'ITSTA_GOT','got_memall',itsta_got)
           allocate(resis_got(2,misgs_got),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'RESIS_GOT','got_memall',resis_got)
        end if
     end if
     !
     ! Residual
     !
     if(postp(1) % npp_stepi (9)>0.and.kfl_probl_got/=2) then
        allocate(cdold_got(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'CDOLD_GOT','got_memall',cdold_got)
     end if
     if(postp(1) % npp_stepi (10)>0.and.kfl_probl_got<=2) then
        allocate(vdold_got(ndime,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'VDOLD_GOT','got_memall',vdold_got)
     end if
     if(kfl_diffu_got==2.and.kfl_probl_got<=2) then
        allocate(diffm_got(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'DIFFM_GOT','got_memall',diffm_got)
     end if
     !
     ! Elemental arrays
     !
     allocate(elcod_got(ndime,mnode),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ELCOD_GOT','got_memall',elcod_got)
     allocate(elvel_got(ndime,mnode),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ELVEL_GOT','got_memall',elvel_got)
     allocate(elvdr_got(ndime,mnode,ncomp_got),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ELVDR_GOT','got_memall',elvdr_got)
     allocate(elcdr_got(mnode,ncomp_got),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ELCDR_GOT','got_memall',elcdr_got) 
     allocate(eldif_got(mnode),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ELDIF_GOT','got_memall',eldif_got) 

     solve(2) = solve(1)
     solve(3) = solve(1)
     solve_sol => solve(1:)
     call soldef(4_ip)

  else
     !
     ! Master: allocate minimum memory
     !
     allocate(vdrop(1,1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'VELOC','memunk',vdrop)
     allocate(cdrop(1,1),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'PRESS','memunk',cdrop)

  end if

end subroutine got_memall
