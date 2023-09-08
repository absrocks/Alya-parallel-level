subroutine par_arrgro(inter)
  !-------------------------------------------------------------------------------
  !****f* parall/par_arrgro
  ! NAME
  !    par_arrgro
  ! DESCRIPTION
  ! INPUT
  ! OUTPUT
  ! USED BY
  !   par_grogro
  !***
  !-------------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_parall
  use def_solver
  use mod_memchk
  implicit none
  integer(ip), intent(in)  :: inter
  integer(ip)              :: offsetI,offsetB,ipart,vv,ii,b_ind,jj
  integer(ip)              :: nbGroupInter,nbGroupBound,igrou,jgrou,kk,jpart
  integer(ip)              :: dsize,ginde(4,npart_par+1)
  integer(4)               :: istat
  integer(ip), pointer     :: lnin2_par(:),invpI(:),invpB(:)
  !
  ! Allocate local memory
  !
  dsize = 2*(npoin/ npart_par)
  allocate(lnin2_par(comle(icoml)%ngrou),stat=istat)
  call memchk(zero,istat,mem_servi(1:2,servi),'LNIN2_PAR','par_memory',lnin2_par)
  allocate(invpI(dsize),stat=istat)
  call memchk(zero,istat,mem_servi(1:2,servi),'invpI','par_arrays',invpI)
  allocate(invpB(dsize),stat=istat)
  call memchk(zero,istat,mem_servi(1:2,servi),'invpB','par_arrays',invpB)
  !
  ! For each subdomain: count total interior and boundary nodes
  !     
  comle(icoml)%gni = 0
  comle(icoml)%gnb = 0
  offsetI          = 0
  offsetB          = inter

  do ipart = 1,npart_par

     nbGroupInter = 1
     nbGroupBound = 1

     do vv = 1,comle(icoml)%ngrou
        if( comle(icoml)%lnpar_par(vv) == ipart ) then
           invpI(nbGroupInter) = vv
           nbGroupInter        = nbGroupInter + 1 

        else if( comle(icoml)%lnpar_par(vv) == -ipart ) then
           invpB(nbGroupBound) = vv
           nbGroupBound        = nbGroupBound + 1

        end if
     end do

     nbGroupInter     = nbGroupInter - 1
     nbGroupBound     = nbGroupBound - 1
     ginde(1,ipart)   = comle(icoml)%gni + 1
     ginde(2,ipart)   = comle(icoml)%gnb + 1
     ginde(3,ipart)   = nbGroupInter
     ginde(4,ipart)   = nbGroupBound
     comle(icoml)%gni = comle(icoml)%gni + nbGroupInter
     comle(icoml)%gnb = comle(icoml)%gnb + nbGroupBound 
     !
     ! Number interior groups
     !
     do ii = 1, nbGroupInter
        kk            = invpI(ii)
        jj            = ii + offsetI
        lnin2_par(jj) = kk
     end do
     offsetI = offsetI + nbGroupInter 
     !
     ! Number boundary groups
     !
     do ii = 1, nbGroupBound
        kk            = invpB(ii)
        jj            = ii + offsetB
        lnin2_par(jj) = kk
     end do
     offsetB = offsetB + nbGroupBound

  end do
  !
  ! BADJ, BDOM, BPOIN
  !
  call par_memgro(3_ip)

  comle(icoml)%badj(1) = 1
  b_ind                = 1
  igrou                = comle(icoml)%gni
  kk                   = 1

  do ipart= 1, npart_par

     do ii= 1, ginde(4,ipart)

        igrou = igrou + 1        ! Interior/boundary numbering (strarting at end of interior)
        jgrou = lnin2_par(igrou) ! Original numbering

        do jj= 1, comle(icoml)%ndomi(jgrou)
           jpart                     = comle(icoml)%domli(jj,jgrou)
           comle(icoml)%bdom(b_ind)  = jpart
           comle(icoml)%bpoin(b_ind) = jgrou
           b_ind                     = b_ind + 1
        end do
        kk                    = kk + 1
        comle(icoml)%badj(kk) = b_ind

     end do

  end do
  !
  ! Deallocate local memory
  !
  call memchk(two,istat,mem_servi(1:2,servi),'INVPB','par_arrgro',invpb)
  deallocate(invpb,stat=istat)
  if(istat/=0) call memerr(two,'invpb','par_arrgro',0_ip)

  call memchk(two,istat,mem_servi(1:2,servi),'INVPI','par_arrgro',invpI)
  deallocate(invpI,stat=istat)
  if(istat/=0) call memerr(two,'invpI','par_arrgro',0_ip)

  call memchk(two,istat,mem_servi(1:2,servi),'LNIN2_PAR','par_arrgro',lnin2_par)
  deallocate(lnin2_par,stat=istat)
  if(istat/=0) call memerr(two,'lnin2_par','par_arrgro',0_ip)

end subroutine par_arrgro
