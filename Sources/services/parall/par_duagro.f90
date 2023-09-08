subroutine par_duagro(nbDual)
!-------------------------------------------------------------------------------
!****f* parall/duagra
! NAME
!    par_duagro
! DESCRIPTION
!
! INPUT
!    nelem
!    ia
!    ja
! OUTPUT
!    nelemDual
!    transl
!    iaDual
!    comle(icoml)%jaDual
! USED BY
!    metisPermute
!***
!-------------------------------------------------------------------------------
  use def_parame
  use def_parall
  use mod_memchk
  use def_master
  implicit none
  integer(ip), intent(out) :: nbDual
  integer(ip)              :: ii, jj, vv, ww, t1, t2
  integer(4)               :: istat
  logical(lg)              :: fin
  integer(ip), allocatable :: iwa(:)
  !
  ! Construir vector para asignar un id a cada arista y asi traducir las 
  ! aristas como nodos del nuevo grafo 
  !
  allocate(comle(icoml)%translDual(comle(icoml)%xadjDom(npart_par+1)-1),stat=istat)
  call memchk(zero,istat,mem_servi(1:2,servi),'comle(icoml)%translDual','par_metis',comle(icoml)%translDual)
  !
  ! Edge numbering for each adjacency: TRANSDUAL
  !
  nbDual = 0
  do vv = 1, npart_par
     do ii = comle(icoml)%xadjDom(vv), comle(icoml)%xadjDom(vv+1)-1
        ww = comle(icoml)%adjDom(ii)
        if (vv < ww) then
           !
           ! Add unique edge vv-ww ( domain vv < domain ww )
           !
           nbDual = nbDual + 1
           comle(icoml)%translDual(ii) = nbDual
        else
           !
           ! Point to already existing edge vv-ww
           !
           fin = .false.
           jj  = comle(icoml)%xadjDom(ww)
           do while ( (jj < comle(icoml)%xadjDom(ww+1) ) .and. ( .not. fin) )
              if ( comle(icoml)%adjDom(jj) == vv) then
                 comle(icoml)%translDual(ii) = comle(icoml)%translDual(jj)
                 fin = .true.             
              endif
              jj = jj + 1
           enddo
        endif
     enddo
  enddo

  !----------------------------------------------------------------------
  !
  ! Dual graph JADUAL and IADUAL
  !
  !----------------------------------------------------------------------

  !
  ! Calcular tamaño comle(icoml)%jaDual (nodos)
  !
  allocate(iwa(nbDual),stat=istat)
  call memchk(zero,istat,mem_servi(1:2,servi),'iwa','par_metis',iwa)

  do ii = 1, nbDual
     iwa(ii) = 0
  enddo

  do ii = 1, comle(icoml)%xadjDom(npart_par+1)-1
     vv = comle(icoml)%adjDom(ii)
     t1 = comle(icoml)%translDual(ii)
     do jj = comle(icoml)%xadjDom(vv), comle(icoml)%xadjDom(vv+1)-1
        if ( comle(icoml)%translDual(jj) /= t1 ) then
           iwa(t1) = iwa(t1) + 1
        endif
     enddo
  enddo
  !
  ! Rellenar comle(icoml)%iaDual (nodos)
  !
  allocate(comle(icoml)%iaDual(nbDual+1),stat=istat)
  call memchk(zero,istat,mem_servi(1:2,servi),'comle(icoml)%iaDual','par_metis',comle(icoml)%iaDual)

  comle(icoml)%iaDual(1) = 1
  do ii = 1, nbDual
     comle(icoml)%iaDual(ii+1)   = comle(icoml)%iaDual(ii) + iwa(ii)
     iwa(ii) = comle(icoml)%iaDual(ii)
  enddo
  !
  ! Rellenar comle(icoml)%jaDual (no ordenado)
  !
  allocate(comle(icoml)%jaDual(comle(icoml)%iaDual(nbDual+1)-1),stat=istat)
  call memchk(zero,istat,mem_servi(1:2,servi),'comle(icoml)%jaDual','par_metis',comle(icoml)%jaDual)

  do ii = 1, comle(icoml)%xadjDom(npart_par+1)-1
     vv = comle(icoml)%adjDom(ii)
     t1 = comle(icoml)%translDual(ii)
     do jj = comle(icoml)%xadjDom(vv), comle(icoml)%xadjDom(vv+1)-1
        t2 = comle(icoml)%translDual(jj)
        if (t1 /= t2) then
           comle(icoml)%jaDual(iwa(t1)) = t2
           iwa(t1)         = iwa(t1) + 1
        endif
     enddo
  enddo

  call memchk(two,istat,mem_servi(1:2,servi),'iwa','par_duagro',iwa)
  deallocate(iwa,stat=istat)
  if(istat/=0) call memerr(two,'iwa','par_duagro',0_ip)

end subroutine par_duagro
