subroutine renum3(nadj,nact,noda,newnn,nodes,i)
!-----------------------------------------------------------------------
!
! Resequence nodes for minimum profile
!
!-----------------------------------------------------------------------
  use def_kintyp
  implicit none
  integer(ip) :: nodes,i
  integer(ip) :: nadj(1:*),newnn(1:*),nact(1:*),noda(1:*)
  integer(ip) :: large,j,nac,maxfrt,npj,nad,k,minnew,lmin,iact
  integer(ip) :: mini,newi,n,next,nif,npos,ilast

  nif   = 0
  large = 5**5
!
! King's scheme
!
  do j = 1,nodes
     newnn(j) = 0
     noda(j) = 0
  end do
  newnn(i) = 1
  nac = 0
!
! Negate all ndeg entries for nodes which are
! adjacent to starting node i
!
  maxfrt = 0
  nad = nadj(i)
  do while (nad>0)
     maxfrt = maxfrt + 1
     npj = nadj(nad)
     if(noda(npj)==0) then
        nac = nac + 1
        noda(npj) = nac
        nact(nac) = npj
     end if
     nad = nadj(nad+1)
  end do
  noda(i) = large
!
! Loop over nodes to be renumbered
!
  do k = 2,nodes
     minnew = large
     lmin = large
!
! Loop over active nodes
! Skip to next node if old node is already renumbered
!
     do iact = 1,nac
        j = nact(iact)
        if(newnn(j)<=0) then
           newi = -1
           mini = large
!
! Compute the increment in active nodes for each node j
! Compute when this node was first active by checking for renumbered
! neighbours with lowest numbers
!
           nad = nadj(j)
           do while (nad>0)
              n = nadj(nad)
              if(noda(n)==0) newi = newi + 1
              if(newnn(n)/=0) then
                 if(newnn(n)<mini)mini = newnn(n)
              end if
              nad = nadj(nad+1)
           end do
!
! Select node with smallest increment in active nodes
! in the case of a tie, select node which has been longest active
!
           if(newi<=minnew) then
              if((newi/=minnew).or.(mini<lmin)) then
                 minnew = newi
                 lmin = mini
                 next = j
              end if
           end if
        end if
     end do
!
! Renumber node and compute number of active nodes
!
     newnn(next) = k
     nif = nif+minnew
     if(nif>maxfrt) maxfrt = nif
!
! Set nodes which are adjacent to the node just renumbered
! as actives nodes, deactivate next
!
     npos = noda(next)
     ilast = nact(nac)
     noda(ilast) = npos
     nact(npos) = ilast
     nac = nac - 1
     
     if(minnew /= -1) then
        nad = abs(nadj(next))
        do while (nad>0)
           n = nadj(nad)
           if(noda(n)==0) then
              nac = nac + 1
              noda(n) = nac
              nact(nac) = n
           end if
           nad = nadj(nad+1)
        end do
     end if
  end do
  
end subroutine renum3
