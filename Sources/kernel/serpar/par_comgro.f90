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
  use mod_memory
  use def_master
  use mod_parall, only : par_memor
  implicit none
  integer(ip),    intent(in)  :: nbdual
  integer(ip)                 :: adjncy,ii,jj,vv,ww
  integer(4)                  :: istat
  type(tAdj_par), pointer     :: invTransl(:)

  nullify(invTransl)
  call par_memgro(4_ip)
  call memory_alloca(par_memor,'INVTRANSL','par_comgro',invTransl,nbdual,'DO_NOT_INITIALIZE')
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
  
  call memory_deallo(par_memor,'INVTRANSL','par_comgro',invTransl)

end subroutine par_comgro

