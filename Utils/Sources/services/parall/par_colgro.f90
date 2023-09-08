subroutine par_colgro(ndual_par)
  !-------------------------------------------------------------------------------
  !****f* parall/par_colgro
  ! NAME
  !    par_colgro
  ! DESCRIPTION
  !
  ! INPUT
  !    nelem
  !    ia
  !    ja
  ! OUTPUT
  !    nbColours
  !    colour
  ! USED BY
  !    metisPermute
  !***
  !-------------------------------------------------------------------------------
  use def_parame
  use def_parall
  use mod_memchk
  use def_master
  implicit none
  !--------------------------------------------------------------- Input Variables
  integer(ip) :: ndual_par
  !--------------------------------------------------------------- Local Variables
  integer(ip),parameter :: maxColors = 50
  integer(ip) :: colNode, ii, vv, ww, node
  integer(4)  :: istat
  logical     :: colourFound

  type :: tSortVect 
     integer(ip) :: node
     integer(ip) :: nbAdj
  end type tSortVect

  type(tSortVect) :: aux
  type(tSortVect),allocatable :: sortVec(:)
  integer(ip),    allocatable :: mask(:)
  !------------------------------------------------------------------------- BEGIN

  allocate(comle(icoml)%colours(ndual_par),stat=istat)
  call memchk(zero,istat,mem_servi(1:2,servi),'comle(icoml)%colours','par_metis',comle(icoml)%colours)

  allocate(sortVec(ndual_par),stat=istat)
  !call memchk(zero,istat,mem_servi(1:2,servi),'sortVec','par_metis',sortVec)

  allocate(mask(0:maxColors),stat=istat)
  call memchk(zero,istat,mem_servi(1:2,servi),'mask','par_metis',mask)

  do vv = 1, ndual_par
     sortVec(vv)%node   = vv
     sortVec(vv)%nbAdj  = comle(icoml)%iaDual(vv+1) - comle(icoml)%iaDual(vv)
  enddo

  !     Ordenar sortVec, de momento burbuja
  do vv = ndual_par, 1, -1
     do ww = 1, vv-1
        if (sortVec(ww)%nbAdj > sortVec(ww+1)%nbAdj) then
           aux = sortVec(ww) 
           sortVec(ww) = sortVec(ww+1)
           sortVec(ww+1) = aux
        endif
     enddo
  enddo

  !     Asignación de colores a cada nodo 
  do vv = 1, ndual_par
     comle(icoml)%colours(vv) = 0
  enddo
  comle(icoml)%nbcol = 0

  !     Usar minimo numero de colores, de tal manera que cada nodo no tenga el mismo color que los nodos vecinos 
  do ii= 0, maxColors
     mask(ii) = 0
  enddo

  comle(icoml)%nbcol = 1
  comle(icoml)%colours(sortVec(1)%node) = 1

  DO vv= 2, ndual_par
     node = sortVec(vv)%node

     do ii= comle(icoml)%iaDual(node), comle(icoml)%iaDual(node+1)-1
        mask(comle(icoml)%colours(comle(icoml)%jaDual(ii))) = vv
     enddo

     colourFound = .false.
     ii= 1
     do while (ii<=comle(icoml)%nbcol .and. .not.colourFound)
        if (mask(ii) /= vv) then
           colNode        = ii
           colourFound = .true.
        endif

        ii = ii + 1
     enddo

     if (.not.colourFound) then
        comle(icoml)%nbcol = comle(icoml)%nbcol + 1
        colNode   = comle(icoml)%nbcol
     endif

     comle(icoml)%colours(sortVec(vv)%node) = colNode
  ENDDO

  call memchk(two,istat,mem_servi(1:2,servi),'mask','par_colgro',mask)
  deallocate(mask,stat=istat)
  if(istat/=0) call memerr(two,'mask','par_colgro',0_ip)
  !call memchk(two,istat,mem_servi(1:2,servi),'sortVec','par_colgro',sortVec)
  deallocate(sortVec,stat=istat)
  !if(istat/=0) call memerr(two,'sortVec','par_colgro',0_ip)

  !--------------------------------------------------------------------------- END
end subroutine par_colgro
