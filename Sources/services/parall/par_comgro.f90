subroutine par_comgro(nbdual)
  !-------------------------------------------------------------------------------
  !****f* parall/par_comgro
  ! NAME
  !    par_comgro
  ! DESCRIPTION
  !
  ! INPUT
  !    npart_par
  !    nbcol
  !    nbdual
  !    ia
  !    ja
  !    transl
  !    colourDual
  ! OUTPUT
  !    communSort
  ! USED BY
  !    metisPermute
  !***
  !-------------------------------------------------------------------------------
  use def_parame
  use def_parall
  use def_master
  use mod_memchk
  use mod_par_memchk
  use def_master
  implicit none
  integer(ip),    intent(in)  :: nbdual
  integer(ip)                 :: adjncy,ii,jj,vv,ww
  integer(4)                  :: istat
  type(tAdj_par), allocatable :: invTransl(:)

  call par_memgro(4_ip)
  allocate(invTransl(nbdual),stat=istat)
  call par_memchk(zero,istat,mem_servi(1:2,servi),'INVTRANSL','par_comgro',invTransl)
  !
  ! Construir la inversa de transl. For each edge number JJ, keep VV-WW
  !
  do vv = 1, npart_par
     do ii = comle(icoml)%xadjDom(vv), comle(icoml)%xadjDom(vv+1)-1
        ww = comle(icoml)%adjDom(ii)
        if ( vv < ww ) then
           jj                  = comle(icoml)%translDual(ii)
           invTransl(jj)%node1 = vv
           invTransl(jj)%node2 = ww
        endif
     enddo
  enddo
  !
  ! inicializar tabla de comunicaciones
  !
  do ii = 1, comle(icoml)%nbcol
     do jj = 1, npart_par
        comle(icoml)%lcomm_par(ii,jj) = 0
     enddo
  enddo
  !
  ! construir la tabla de comunicaciones:
  ! VV communication with WW at iteration II
  ! WW communication with VV at iteration II
  !
  do ii = 1 , comle(icoml)%nbcol
     do adjncy = 1 , nbdual
        if ( comle(icoml)%colours(adjncy) == ii ) then
           vv                            = invTransl(adjncy)%node1  
           ww                            = invTransl(adjncy)%node2
           comle(icoml)%lcomm_par(ii,vv) = ww
           comle(icoml)%lcomm_par(ii,ww) = vv
        endif
     enddo
  enddo

  call par_memchk(two,istat,mem_servi(1:2,servi),'INVTRANSL','par_comgro',invTransl)
  deallocate(invTransl,stat=istat)
  if(istat/=0) call memerr(two,'invTransl','par_comgro',0_ip)

end subroutine par_comgro

