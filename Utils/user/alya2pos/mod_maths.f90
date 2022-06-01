module mod_maths

  use def_kintyp
  implicit none
  
contains

  subroutine maths_heapsort_real(itask,ndofn,nrows,ivin,ivou,SAVING)

    integer(ip),  intent(in)             :: itask
    integer(ip),  intent(in)             :: ndofn
    integer(ip),  intent(in)             :: nrows
    integer(ip),  intent(inout)          :: ivin(nrows)
    real(rp),     intent(inout)          :: ivou(ndofn,*)
    logical(lg),  intent(in),   optional :: SAVING
    integer(ip),  pointer                :: ivin_sav(:)
    integer(ip)                          :: len,ir,ii,jj,iaux,krows
    real(rp)                             :: iau1(ndofn)
    logical(lg)                          :: if_saving

    if_saving = .false.
    if( present(SAVING) ) if_saving = SAVING
    if( if_saving ) then
       allocate(ivin_sav(nrows))
       ivin_sav(1:nrows)=ivin(1:nrows)
    end if
       
    select case ( itask )

    case ( 1_ip )
       !
       ! Decreasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       len = (nrows/2) + 1
       ir  = nrows

100    continue

       if( len > 1 ) then
          len = len - 1
          iaux = ivin(len)
          iau1(1:ndofn) = ivou(1:ndofn,len)
       else
          iaux = ivin(ir)
          iau1(1:ndofn) = ivou(1:ndofn,ir)
          ivin(ir) = ivin(1)
          ivou(1:ndofn,ir) = ivou(1:ndofn,1)

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             ivou(1:ndofn,1) = iau1(1:ndofn)
             goto 500
          end if
       end if

       ii = len
       jj = len + len

200    if( jj <= ir ) then
          if( jj < ir ) then
             if( ivin(jj) > ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux > ivin(jj) ) then
             ivin(ii) = ivin(jj)
             ivou(1:ndofn,ii) = ivou(1:ndofn,jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          endif

          goto 200
       end if

       ivin(ii) = iaux
       ivou(1:ndofn,ii) = iau1(1:ndofn)

       goto 100

    case ( 2_ip )
       !
       ! Increasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       len = (nrows/2) + 1
       ir  = nrows

300    continue

       if( len > 1 ) then
          len              = len - 1
          iaux             = ivin(len)
          iau1(1:ndofn)    = ivou(1:ndofn,len)
       else
          iaux             = ivin(ir)          
          ivin(ir)         = ivin(1)
          iau1(1:ndofn)    = ivou(1:ndofn,ir)
          ivou(1:ndofn,ir) = ivou(1:ndofn,1)

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1)         = iaux
             ivou(1:ndofn,1) = iau1(1:ndofn)
             goto 500
          end if
       end if

       ii = len
       jj = len + len

400    if( jj <= ir ) then
          if( jj < ir ) then
             if ( ivin(jj) < ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if
          
          if( iaux < ivin(jj) ) then
             ivin(ii)         = ivin(jj)
             ivou(1:ndofn,ii) = ivou(1:ndofn,jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          end if
          
          goto 400
       end if
       
       ivin(ii)         = iaux
       ivou(1:ndofn,ii) = iau1(1:ndofn)
       
       goto 300
       
    end select

500 continue

    if( if_saving ) then
       ivin(1:nrows)=ivin_sav(1:nrows)
       deallocate(ivin_sav)
    end if

  end subroutine maths_heapsort_real

  subroutine maths_heapsort_int(itask,ndofn,nrows,ivin,ivou,SAVING)

    integer(ip),  intent(in)              :: itask
    integer(ip),  intent(in)              :: ndofn
    integer(ip),  intent(in)              :: nrows
    integer(ip),  intent(inout)           :: ivin(*)
    integer(ip),  intent(inout), optional :: ivou(ndofn,*)
    logical(lg),  intent(in),    optional :: SAVING
    integer(ip),  pointer                 :: ivin_sav(:)
    integer(ip)                           :: len,ir,ii,jj,iaux,krows
    integer(ip)                           :: iau1(ndofn)
    logical(lg)                           :: if_saving

    if_saving = .false.
    if( present(SAVING) ) if_saving = SAVING
    if( if_saving ) then
       allocate(ivin_sav(nrows))
       ivin_sav(1:nrows)=ivin(1:nrows)
    end if

    select case ( itask )

    case ( 1_ip )
       !
       ! Decreasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       len = (nrows/2) + 1
       ir  = nrows

100    continue

       if( len > 1 ) then
          len = len - 1
          iaux = ivin(len)
          if(present(ivou) ) iau1(:) = ivou(:,len)
       else
          iaux = ivin(ir)
          if(present(ivou) ) iau1(:) = ivou(:,ir)
          ivin(ir) = ivin(1)
          if(present(ivou) ) ivou(:,ir) = ivou(:,1)

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             if(present(ivou) ) ivou(:,1) = iau1(:)
             goto 500
          end if
       end if

       ii = len
       jj = len + len

200    if( jj <= ir ) then
          if( jj < ir ) then
             if( ivin(jj) > ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux > ivin(jj) ) then
             ivin(ii) = ivin(jj)
             if(present(ivou) ) ivou(:,ii) = ivou(:,jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          endif

          goto 200
       end if

       ivin(ii) = iaux
       if(present(ivou) ) ivou(:,ii) = iau1(:)

       goto 100

    case ( 2_ip )
       !
       ! Increasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       len = (nrows/2) + 1
       ir  = nrows

300    continue

       if( len > 1 ) then
          len = len - 1
          iaux = ivin(len)
          if(present(ivou) ) iau1(:) = ivou(:,len)
       else
          iaux     = ivin(ir)
          ivin(ir) = ivin(1)
          iau1(:)     = ivou(:,ir)
          if(present(ivou) ) ivou(:,ir) = ivou(:,1)

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             if(present(ivou) ) ivou(:,1) = iau1(:)
             goto 500
          end if
       end if

       ii = len
       jj = len + len


400    if( jj <= ir ) then
          if( jj < ir ) then
             if ( ivin(jj) < ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if
          
          if( iaux < ivin(jj) ) then
             ivin(ii) = ivin(jj)
             if(present(ivou) ) ivou(:,ii) = ivou(:,jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          end if
          
          goto 400
       end if
       
       ivin(ii) = iaux
       if(present(ivou) ) ivou(:,ii) = iau1(:)
       
       goto 300
       
    end select

500 continue
    
    if( if_saving ) then
       ivin(1:nrows)=ivin_sav(1:nrows)
       deallocate(ivin_sav)
    end if

  end subroutine maths_heapsort_int

  subroutine maths_heap_sort_I1(itask,nrows,ivin,message,ivo1,ivo2,PERMUTATION)

    integer(ip),  intent(in)            :: itask
    integer(ip),  intent(inout)         :: nrows
    integer(ip),  intent(inout)         :: ivin(*)
    character(*), intent(in),  optional :: message
    integer(ip),  intent(inout), optional :: ivo1(*)
    integer(ip),  intent(inout), optional :: ivo2(*)
    integer(ip),  intent(inout), optional :: PERMUTATION(*)
    integer(ip)                         :: len,ir,ii,jj,iaux,krows
    integer(ip)                         :: iau1,iau2

    select case ( itask )

    case ( 1_ip )
       !
       ! Decreasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       len = (nrows/2) + 1
       ir  = nrows

100    continue

       if( len > 1 ) then
          len = len - 1
          iaux = ivin(len)
          if( present(ivo1) ) iau1 = ivo1(len)
          if( present(ivo2) ) iau2 = ivo2(len)
       else
          iaux = ivin(ir)
          if( present(ivo1) ) iau1 = ivo1(ir)
          if( present(ivo2) ) iau2 = ivo2(ir)
          ivin(ir) = ivin(1)
          if( present(ivo1) ) ivo1(ir) = ivo1(1)
          if( present(ivo2) ) ivo2(ir) = ivo2(1)

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             if( present(ivo1) ) ivo1(1) = iau1
             if( present(ivo2) ) ivo2(1) = iau2
             goto 500
          end if
       end if

       ii = len
       jj = len + len

200    if( jj <= ir ) then
          if( jj < ir ) then
             if( ivin(jj) > ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux > ivin(jj) ) then
             ivin(ii) = ivin(jj)
             if( present(ivo1) ) ivo1(ii) = ivo1(jj)
             if( present(ivo2) ) ivo2(ii) = ivo2(jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          endif

          goto 200
       end if

       ivin(ii) = iaux
       if( present(ivo1) ) ivo1(ii) = iau1
       if( present(ivo2) ) ivo2(ii) = iau2

       goto 100

    case ( 2_ip )
       !
       ! Increasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       len = (nrows/2) + 1
       ir  = nrows

300    continue

       if( len > 1 ) then
          len = len - 1
          iaux = ivin(len)
          if( present(ivo1) ) iau1 = ivo1(len)
          if( present(ivo2) ) iau2 = ivo2(len)
       else
          iaux     = ivin(ir)
          ivin(ir) = ivin(1)
          if( present(ivo1) ) then
             iau1     = ivo1(ir)
             ivo1(ir) = ivo1(1)
          end if
          if( present(ivo2) ) then
             iau2     = ivo2(ir)
             ivo2(ir) = ivo2(1)
          end if

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             if( present(ivo1) ) ivo1(1) = iau1
             if( present(ivo2) ) ivo2(1) = iau2
             goto 500
          end if
       end if

       ii = len
       jj = len + len

       if( present(PERMUTATION) ) then
401       if( jj <= ir ) then
             if( jj < ir ) then
                if ( PERMUTATION(ivin(jj)) < PERMUTATION(ivin(jj+1)) ) then
                   jj = jj + 1
                end if
             end if

             if( PERMUTATION(iaux) < PERMUTATION(ivin(jj)) ) then
                ivin(ii) = ivin(jj)
                if( present(ivo1) ) ivo1(ii) = ivo1(jj)
                if( present(ivo2) ) ivo2(ii) = ivo2(jj)
                ii = jj
                jj = jj + jj
             else
                jj = ir + 1
             end if

             goto 401
          end if          
       else
400       if( jj <= ir ) then
             if( jj < ir ) then
                if ( ivin(jj) < ivin(jj+1) ) then
                   jj = jj + 1
                end if
             end if

             if( iaux < ivin(jj) ) then
                ivin(ii) = ivin(jj)
                if( present(ivo1) ) ivo1(ii) = ivo1(jj)
                if( present(ivo2) ) ivo2(ii) = ivo2(jj)
                ii = jj
                jj = jj + jj
             else
                jj = ir + 1
             end if

             goto 400
          end if
       end if

       ivin(ii) = iaux
       if( present(ivo1) ) ivo1(ii) = iau1
       if( present(ivo2) ) ivo2(ii) = iau2

       goto 300

    end select

    return
    !
    ! Eliminate duplicates
    !
500 continue
    if( present(message) ) then
       if( trim(message) == 'ELIMINATE DUPLICATES') then
          krows = nrows
          do ii = 1,krows-1
             if( ivin(ii+1) == ivin(ii) ) then
                do jj = ii+1,krows-1
                   ivin(jj) = ivin(jj+1)
                   if( present(ivo1) ) ivo1(jj) = ivo1(jj+1)
                   if( present(ivo2) ) ivo2(jj) = ivo2(jj+1)
                end do
                krows = krows - 1
             end if
          end do
          nrows = krows
       end if
    end if


  end subroutine maths_heap_sort_I1

  subroutine maths_heapsort(itask,nrows,ndofn,ivin,ivo1)

    integer(ip),  intent(in)    :: itask
    integer(ip),  intent(inout) :: nrows
    integer(ip),  intent(in)    :: ndofn
    integer(ip),  intent(inout) :: ivin(*)
    real(rp),     intent(out)   :: ivo1(ndofn,*)
    integer(ip)                 :: len,ir,ii,jj,iaux,krows
    real(rp)                    :: iau1(ndofn)
    !integer(ip), allocatable :: ivin_cpy(:)

    !allocate(ivin_cpy(nrows))
    !ivin_cpy(1:nrows) = ivin(1:nrows)
    
    select case ( itask )

    case ( 1_ip )
       !
       ! Decreasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       len = (nrows/2) + 1
       ir  = nrows

100    continue

       if( len > 1 ) then
          len = len - 1
          iaux = ivin(len)
          iau1(:) = ivo1(:,len)
       else
          iaux = ivin(ir)
          iau1(:) = ivo1(:,ir)
          ivin(ir) = ivin(1)
          ivo1(:,ir) = ivo1(:,1)
 
          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             ivo1(:,1) = iau1(:)
             goto 500
          end if
       end if

       ii = len
       jj = len + len

200    if( jj <= ir ) then
          if( jj < ir ) then
             if( ivin(jj) > ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux > ivin(jj) ) then
             ivin(ii) = ivin(jj)
             ivo1(:,ii) = ivo1(:,jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          endif

          goto 200
       end if

       ivin(ii) = iaux
       ivo1(:,ii) = iau1(:)

       goto 100

    case ( 2_ip )
       !
       ! Increasing order
       !
       if( nrows < 2 ) then
          goto 500
       end if

       len = (nrows/2) + 1
       ir  = nrows

300    continue

       if( len > 1 ) then
          len = len - 1
          iaux = ivin(len)
          iau1(:) = ivo1(:,len)
       else
          iaux     = ivin(ir)
          ivin(ir) = ivin(1)
          iau1(:)    = ivo1(:,ir)
          ivo1(:,ir) = ivo1(:,1)

          ir = ir - 1

          if( ir == 1 ) then
             ivin(1) = iaux
             ivo1(:,1) = iau1(:)
             goto 500
          end if
       end if

       ii = len
       jj = len + len

400    if( jj <= ir ) then
          if( jj < ir ) then
             if ( ivin(jj) < ivin(jj+1) ) then
                jj = jj + 1
             end if
          end if

          if( iaux < ivin(jj) ) then
             ivin(ii) = ivin(jj)
             ivo1(:,ii) = ivo1(:,jj)
             ii = jj
             jj = jj + jj
          else
             jj = ir + 1
          end if

          goto 400
       end if

       ivin(ii) = iaux
       ivo1(:,ii) = iau1(:)

       goto 300

    end select

500 continue
    !ivin(1:nrows) = ivin_cpy(1:nrows)
    
  end subroutine maths_heapsort
  
end module mod_maths

