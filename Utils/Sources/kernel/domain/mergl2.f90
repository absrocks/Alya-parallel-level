subroutine mergl2( kk, lista, lsize, nodes, nnode, me )
  !-----------------------------------------------------------------------
  !****f* Domain/mergl2
  ! NAME
  !    mergl2
  ! DESCRIPTION
  !    This routine merges to list of nodes
  ! OUTPUT
  !    LISTA
  !    LSIZE
  ! USED BY
  !    domgra
  !***
  !-----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip), intent(inout) :: lsize,lista(*)
  integer(ip), intent(in)    :: kk,nnode,me
  integer(ip), intent(in)    :: nodes(nnode)
  integer(ip)                :: ii, jj, n1, n2
  logical(lg)                :: noEncontrado

  if (me<0) then
     do ii= 1, nnode
        n1 = nodes(ii)
        if( n1 /= kk ) then
           jj = 1
           noEncontrado = .true.
           do while( jj<=lsize .and. noEncontrado)
              n2 = lista(jj)
              if (n1==n2) then
                 noEncontrado = .false.
              end if
              jj = jj + 1
           end do
           if (noEncontrado) then
              lsize = lsize + 1
              lista(lsize) = n1
           end if
        end if
     end do
  else
     do ii= 1, nnode
        n1 = nodes(ii)
        if( n1 /= kk ) then
           if (n1/=me) then
              jj = 1
              noEncontrado = .true.
              do while( jj<=lsize .and. noEncontrado)
                 n2 = lista(jj)
                 if (n1==n2) then
                    noEncontrado = .false.
                 endif
                 jj = jj + 1
              enddo
              if (noEncontrado) then
                 lsize = lsize + 1
                 lista(lsize) = n1
              endif
           endif
        end if
     enddo
  end if

end subroutine mergl2
