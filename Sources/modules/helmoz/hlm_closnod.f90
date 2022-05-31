subroutine hlm_closnod(tpoin,nn,clnod1,clnod2)

  !-----------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_closnod.f90
  ! NAME
  !    hlm_closnod
  ! DESCRIPTION
  !    This routine finds N closest nodes of a mesh to a test point.
  ! INPUT ARGUMENTS
  !    TPOIN ... Cartesian coordinates of a test point (xt, yt, zt)
  !    NN ...... Number of closest mesh nodes to a test point, N
  ! OUTPUT ARGUMENTS
  !    CLNOD1 ... Global numbers of N closest mesh nodes to a test point
  !    CLNOD2 ... Distances between each of the N closest mesh nodes 
  !               and a test point
  ! USES
  ! USED BY
  !    hlm_inivar
  !-----------------------------------------------------------------------

  use def_master
  use def_domain
  use def_helmoz
  use def_parame
  use def_kermod
  use mod_elsest,         only : elsest_host_element

  implicit none

  real(rp),    intent(in)  :: tpoin(3)
  integer(ip), intent(in)  :: nn
  integer(ip), intent(out) :: clnod1(nn)
  real(rp),    intent(out) :: clnod2(nn)

  logical     :: FINISH	
  integer(ip) :: ii,jj,ielem,pelty,pnode,inode,ipoin,nn2,start,fin,value1
  real(rp)    :: shapf(mnode),deriv(ndime,mnode),coloc(3),x,y,z,value2
  integer(ip) :: aux1(4_ip*nn)
  real(rp)    :: aux2(4_ip*nn),dista

  !Find the host element for a test point
  !call elsest(2_ip,1_ip,ielse,mnode,ndime,npoin,nelem,nnode(1:),pelpo,&
  !            lelpo,lnods,ltype,ltopo,coord,tpoin,relse,ielem,shapf,deriv,coloc,dummi)

  !print *,'size(coord)',size(coord)
  !print *,'coord?',coord(1,1)
 
  call elsest_host_element(&
       ielse,relse,1_ip,meshe(ndivi),tpoin,ielem,&
       shapf,deriv,coloc,dista)

  !call elsest(2_ip,1_ip,ielse,mnode,ndime,npoin,nelem,nnode(1:),&
  !            lnods,ltype,ltopo,coord,tpoin,relse,ielem,&
  !            shapf,deriv,coloc,dummi)
  !write(*,*) 'sale' 


  if (ielem == 0_ip) call runend('COULD NOT FIND THE HOST ELEMENT FOR A TEST POINT')
  pelty = ltype(ielem)       !Type of the host element
  pnode = nnode(pelty)       !Number of nodes of the host element
  do inode = 1,pnode
    aux1(inode) = lnods(inode,ielem)       !Put the global number of a node of the host element in array
  enddo
  nn2 = 4_ip * nn
  start = 1_ip      
  fin = pnode
  call hlm_recu(nn2,aux1,start,fin)     !Create array of NN2 closest mesh nodes 
  !Create array of distances between each of the NN2 closest mesh nodes and a test point 
  do ii = 1,nn2
    x = coord(1,aux1(ii)) - tpoin(1)
    y = coord(2,aux1(ii)) - tpoin(2)
    z = coord(3,aux1(ii)) - tpoin(3)
    aux2(ii) = sqrt(x*x+y*y+z*z)
!    if (airpl_hlm /= 1.0e30_rp .and. ((coord(3,aux1(ii))<airpl_hlm .and. tpoin(3)>airpl_hlm) .or. (coord(3,aux1(ii))>airpl_hlm .and. tpoin(3)<airpl_hlm))) then
!        aux2(ii) = 1.0e30_rp            !Ban the point since it's on the other side of the air-plane
!		write (*,*) 'Point banned for', tpoin(1), tpoin(2), tpoin(3)
!    endif
  enddo
  !Insertion sort of array of distances between each of the NN2 closest mesh nodes and a test point (aux2) and
  !of array of NN2 closest mesh nodes (aux1)
  do ii = 2,nn2
	value1 = aux1(ii)
	value2 = aux2(ii)

	jj = ii - 1_ip
	FINISH = .false.

	do while (.not. FINISH)
		if (aux2(jj) > value2) then
			aux1(jj+1) = aux1(jj)
			aux2(jj+1) = aux2(jj)
			jj = jj - 1_ip
			if (jj < 1_ip) FINISH = .true.
		else
			FINISH = .true.
		endif             
	enddo 
	aux1(jj+1) = value1
	aux2(jj+1) = value2
  enddo

  !Create array of N closest mesh nodes
  do ii = 1,nn
    clnod1(ii) = aux1(ii)
    clnod2(ii) = aux2(ii)
  enddo

end subroutine hlm_closnod

recursive subroutine hlm_recu(nn2,aux,start,fin)

  !-----------------------------------------------------------------------
  ! NAME
  !    hlm_recu
  ! DESCRIPTION
  !    This routine creates array of NN2 closest mesh nodes to a test point
  ! INPUT/OUTPUT ARGUMENTS
  !    NN2 ..... Number of mesh nodes which will be checked for being closest to a test point 
  !    AUX ..... Array of NN2 closest mesh nodes
  ! USES
  ! USED BY
  !    hlm_closnod
  !-----------------------------------------------------------------------

  use def_master
  use def_domain
  use def_helmoz
  use def_parame

  implicit none

  integer(ip), intent(in)    :: nn2
  integer(ip), intent(inout) :: aux(nn2)
  integer(ip), intent(inout) :: start,fin

  integer(ip) :: ii,jj,kk,ielem,pelty,pnode,inode,ipoin,fin1,fin2,already

  fin1 = fin
  fin2 = fin
  points: do ii = start,fin
	elements: do jj = pelpo(aux(ii)),pelpo(aux(ii)+1_ip)-1_ip
		ielem = lelpo(jj)
		pelty = ltype(ielem)
		pnode = nnode(pelty)
		nodes: do inode = 1,pnode
			ipoin = lnods(inode,ielem)
			already = 0_ip
			search: do kk = 1,fin2
				if (aux(kk) == ipoin) then
					already = 1_ip
					exit search
				endif
			enddo search
			if (already == 0_ip .and. fin1 < nn2) then
				fin1 = fin1 + 1_ip
				aux(fin1) = ipoin
			endif
			if (fin1 >= nn2) return
		enddo nodes
		fin2 = fin1
	enddo elements
  enddo points
  start = fin + 1_ip
  fin = fin1

  call hlm_recu(nn2,aux,start,fin)

end subroutine hlm_recu
