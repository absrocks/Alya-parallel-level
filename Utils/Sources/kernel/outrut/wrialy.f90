subroutine wrialy(itask,wopos)
  !-----------------------------------------------------------------------
  !****f* outrut/wrialy
  ! NAME 
  !    wrialy
  ! DESCRIPTION
  !    This routine initializes Ensight:
  !    ITASK=1 ... Case for geometry only
  !    ITASK=2 ... Header of the postprocess case
  !    ITASK=3 ... Variable
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use mod_iofile
  use def_postpr
  implicit none
  integer(ip),  intent(in) :: itask
  character(5), intent(in) :: wopos(2)
  integer(ip)              :: istpp

  select case(itask)

  case(2)
     !
     ! Header of the case
     !
     write(nunam_pos,'(i7)') ittim      !   <<-- to write an integer to a character              
     if( ittim < 10 ) then
        write(nunam_pos,'(a,i1)') '000000',ittim
     else if( ittim < 100 ) then
        write(nunam_pos,'(a,i2)') '00000',ittim
     else if( ittim < 1000 ) then
        write(nunam_pos,'(a,i3)') '0000',ittim
     else if( ittim < 10000 ) then
        write(nunam_pos,'(a,i4)') '000',ittim
     else if( ittim < 100000 ) then
        write(nunam_pos,'(a,i5)') '00',ittim
     else if( ittim < 1000000 ) then
        write(nunam_pos,'(a,i6)') '0',ittim
     end if

  case(3)
     !
     ! Check if variable already exists
     !
     istpp = 1
     do while( istpp <= 100 )        
        if( wopos(1) == varna_pos(1,istpp) ) then
           istpp = 100
        else if( varna_pos(1,istpp) == 'NULL' ) then
           varna_pos(1,istpp) = wopos(1)
           varna_pos(2,istpp) = wopos(2)
           varnu_pos          = varnu_pos+1
           istpp = 100
        end if
        istpp = istpp + 1
     end do

  case(4)
     !
     ! Write case end of file
     !
     if( ncoun_pos /= 0 ) then
        nppva_ens     = 0
        kfl_statu_ens = 2
        ncoun_pos     = 0
     end if

  end select

50 format(2a)
60 format(a,3x,i4,3x,a)
70 format(a,4x,i4,4x,a,4x,a)
80 format(a,3x,i4)

end subroutine wrialy
