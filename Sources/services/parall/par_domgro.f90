subroutine par_domgro(inter)
  !-------------------------------------------------------------------------------
  !****f* parall/par_domgro
  ! NAME
  !    par_domgro
  ! DESCRIPTION
  !     Half-matrix with the neighbourhood domains
  !     NEIGHDOM(I,J)... I =J : total # boundary nodes (including 
  !                             repeated communications)
  !                      I/=J : # boundary node shared between I and J
  !     LNPAR_PAR(I) ... Subdomain # if I=internal node
  !                  ... 0           if I=boundary node
  !     LNEIG_PAR(I) ... Number of neighbors of subdomain I
  !     ADJDOM(I) ...... Adjacency pointer for subdomain I
  !     XADJDOM ........ Adjacancy array
  ! INPUT
  !    LEPAR_PAR
  !    PELPO
  !    LELPO
  ! OUTPUT
  !    LNPAR_PAR
  !    NEIGHDOM
  ! USED BY
  !    metisPermute
  !***
  !-------------------------------------------------------------------------------
  use def_parame
  use def_kintyp
  use def_domain
  use def_parall
  use def_solver
  use mod_memchk
  use def_master
  implicit none
  integer(8),  intent(out) :: inter
  integer(ip)              :: igrou,ndomi
  integer(ip)              :: domin, dom1, dom2, ii, jj, kk, idomi, ipoin, ipart
  integer(4)               :: istat
  integer(ip), pointer     :: domli(:)
  logical(lg)              :: ifoun,ispli
  !
  ! Allocate memory
  !
  ispli = .false.
  call par_memgro(1_ip)
  !
  ! Create a half-matrix with the neighbourhood domains (including main diagonal)
  !
  do dom1 = 1, npart_par
     do dom2 = 1, dom1
        comle(icoml)%neighDom( (dom1*(dom1-1))/2 + dom2) = 0
     enddo
  enddo

  allocate(domli(npart_par),stat=istat)
  call memchk(zero,istat,mem_servi(1:2,servi),'domli','par_domgro',domli)

  comle(icoml)%ngrou_par(1:npart_par) = 0

  !front = 0
  inter = 0
  !
  ! NDOMI: number of subdomain to which igrou belongs
  ! DOMLI: list of subdomains to which igrou belongs
  !
  do ipoin = 1,npoin
     igrou   = comle(icoml)%lgrou(ipoin)
     if( igrou > 0 ) then
        ndomi = 0
        call par_domlis( pelpo, lelpo, ipoin, lepar_par, ndomi, domli )
        do idomi = 1,ndomi
           ifoun = .false. 
           popo: do ipart = 1,comle(icoml)%ndomi(igrou)
              if( comle(icoml)%domli(ipart,igrou) == domli(idomi) ) then
                 ifoun= .true.
                 exit popo
              end if
           end do popo
           if( .not. ifoun ) then
              comle(icoml)%ndomi(igrou)                           = comle(icoml)%ndomi(igrou) + 1 
              comle(icoml)%domli(comle(icoml)%ndomi(igrou),igrou) = domli(idomi) 
           end if
        end do
     end if
  end do
  !
  ! NEIGHDOM: number of boundary groups in common for DOM1 and DOM2 (index KK)
  !
  do igrou = 1, comle(icoml)%ngrou

     if ( comle(icoml)%ndomi(igrou) == 1 ) then
        ! This is a internal group
        inter = inter + 1

        domin                         = comle(icoml)%domli(1,igrou)
        comle(icoml)%lnpar_par(igrou) = domin
        comle(icoml)%ngrou_par(domin) = comle(icoml)%ngrou_par(domin) + 1  
    
     else
        ! This is a boundary group
        !      iff i.ne.j   boundary nodes shared by domains i and j
        !      iff i.eq.j   boundary nodes of domain i
        !front = front + 1
        do ii= 1, comle(icoml)%ndomi(igrou)
           dom1 = comle(icoml)%domli(ii,igrou)
           if( ispli ) then
              comle(icoml)%lnpar_par(igrou) = 0
           else
              comle(icoml)%lnpar_par(igrou) = -dom1
           end if
           comle(icoml)%ngrou_par(dom1)  = comle(icoml)%ngrou_par(dom1) + 1
           do jj= 1, ii
              dom2 = comle(icoml)%domli(jj,igrou)
              if (dom1>dom2) then
                 kk = (dom1*(dom1-1))/2 + dom2
              else
                 kk = (dom2*(dom2-1))/2 + dom1
              endif
              comle(icoml)%neighDom(kk) = comle(icoml)%neighDom(kk) + 1
           enddo
        enddo

     endif

  enddo

  call memchk(two,istat,mem_servi(1:2,servi),'domli','par_domgro',domli)
  deallocate(domli,stat=istat)
  if(istat/=0) call memerr(two,'domli','par_domgro',0_ip)
  !
  ! The sum of the nodes that will have every domain
  !
  comle(icoml)%ngrou_total = 0
  do domin= 1, npart_par
     comle(icoml)%ngrou_total = comle(icoml)%ngrou_total + comle(icoml)%ngrou_par(domin)
  enddo
  !
  ! Split boundary
  !
  !if( ispli ) call par_bougro( )
  !
  ! Count the number of neighbours per domain 
  !
  do dom1 = 1, npart_par
     comle(icoml)%lneig_par(dom1) = 0
  enddo

  do dom1 = 1, npart_par
     do dom2 = 1, dom1-1
        if (comle(icoml)%neighDom((dom1*(dom1-1))/2 + dom2) /= 0) then
           comle(icoml)%lneig_par(dom1) = comle(icoml)%lneig_par(dom1) + 1
           comle(icoml)%lneig_par(dom2) = comle(icoml)%lneig_par(dom2) + 1
        endif
     enddo
  enddo
  !
  ! Build the domain interconection graph     
  !
  comle(icoml)%xadjDom(1) = 1
  do dom1 = 1, npart_par
     comle(icoml)%xadjDom(dom1+1) = comle(icoml)%xadjDom(dom1) + comle(icoml)%lneig_par(dom1)
  enddo

  call par_memgro(2_ip)
  !
  ! iwa is a vector with indexes to insert adjacencies
  !
  do dom1= 1, npart_par
     do dom2= 1, dom1-1
        if (comle(icoml)%neighDom((dom1*(dom1-1))/2 + dom2) /= 0) then
           comle(icoml)%adjDom(comle(icoml)%xadjDom(dom1)) = dom2
           comle(icoml)%xadjDom(dom1)                      = comle(icoml)%xadjDom(dom1) + 1
           comle(icoml)%adjDom(comle(icoml)%xadjDom(dom2)) = dom1
           comle(icoml)%xadjDom(dom2)                      = comle(icoml)%xadjDom(dom2) + 1
        endif
     enddo
  enddo

  comle(icoml)%xadjDom(1) = 1
  do dom1 = 1, npart_par
     comle(icoml)%xadjDom(dom1+1) = comle(icoml)%xadjDom(dom1) + comle(icoml)%lneig_par(dom1)
  enddo

end subroutine par_domgro
