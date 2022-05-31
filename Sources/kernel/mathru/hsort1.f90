subroutine hsort1(itask,nrows,ivi1,ivou)
  !------------------------------------------------------------------------
  !****f* mathru/hsort1
  ! NAME
  !    hsort1
  ! DESCRIPTION
  !    Quick sorting of IVOU using IVI1. The element in ivi1 are sorting in:
  !    ITASK = 1 ... Decreasing value, i.e., ivi1(1) > ivi1(2) > ...
  !    ITASK = 2 ... Increasing value, i.e., ivi1(1) < ivi1(2) < ...
  ! INPUT
  !    ITASK ... 1,2 for decreasing, increasing order
  !    NROWS ... Size of IVI1
  !    IVI1 .... Array to be ordered
  !    IVOU .... Array to be ordered
  ! OUTPUT
  !    IVI1 .... Ordered array
  !    IVOU .... Ordered array
  ! USED BY
  !    
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg  
  implicit none
  integer(ip), intent(in)    :: itask,nrows
  integer(ip), intent(inout) :: ivi1(*) 
  integer(ip), intent(inout) :: ivou(*) 
  integer(ip)                :: leni, ir, ii, jj, iau1, jaux

  select case(itask)

  case(1_ip)

     call runend('NOT CODED')

  case (2_ip)

     !-------------------------------------------------------------------
     !
     ! Increasing order
     !
     !-------------------------------------------------------------------

     if(nrows<2) then
        return
     end if

     leni = (nrows/2) + 1
     ir  = nrows

300  continue

     if ( leni > 1 ) then
        leni = leni - 1
        iau1 = ivi1(leni)
        jaux = ivou(leni)
     else
        iau1     = ivi1(ir)
        ivi1(ir) = ivi1(1)

        jaux     = ivou(ir)
        ivou(ir) = ivou(1)

        ir = ir - 1

        if (ir==1) then
           ivi1(1) = iau1
           ivou(1) = jaux
           goto 301
        endif
     end if

     ii = leni
     jj = leni + leni

400  if (jj<=ir) then
        if (jj<ir) then
           if ( ivi1(jj)<ivi1(jj+1) ) then
              jj = jj + 1
           endif
        endif

        if (iau1<ivi1(jj) ) then
           ivi1(ii) = ivi1(jj)
           ivou(ii) = ivou(jj)

           ii = jj
           jj = jj + jj
        else
           jj = ir + 1
        endif

        goto 400
     end if

     ivi1(ii) = iau1
     ivou(ii) = jaux

     goto 300

301 return

  end select

end subroutine hsort1
