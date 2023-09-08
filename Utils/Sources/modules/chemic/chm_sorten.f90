subroutine chm_sorten(ndime,nrows,xnerg_chm,lnerg_chm)
  !------------------------------------------------------------------------
  !****f* mathru/chm_sorten
  ! NAME
  !    chm_sorten
  ! DESCRIPTION
  !    Quik sorting. The element in ivin are sorting in:
  !    Decreasing value, i.e., ivin(1) > ivin(2) > ...
  ! INPUT
  ! OUTPUT
  !    IVIN .... Ordered array
  ! USED BY
  !    
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg  
  use def_chemic
  implicit none
  integer(ip), intent(in)    :: ndime,nrows
  real(rp),    intent(inout) :: xnerg_chm(nspec_chm+1,ndime) 
  integer(ip), intent(inout) :: lnerg_chm(nspec_chm) 
  integer(ip)                :: len, ir, ii, jj, idime, icrit
  real(rp)                   :: raux(ndime)
  !
  ! Decreasing order
  !
  do ii = 1,nspec_chm
     lnerg_chm(ii) = ii
  end do
  if( nrows < 2 ) then
     return
  end if

  len = (nrows/2) + 1
  ir  = nrows

100 continue

  if ( len > 1 ) then
     len   = len - 1
     icrit = lnerg_chm(len)
     do idime = 1,ndime
        raux(idime) = xnerg_chm(len,idime)
     end do
  else
     icrit         = lnerg_chm(ir)
     lnerg_chm(ir) = lnerg_chm(1)
     do idime = 1,ndime
        raux(idime)         = xnerg_chm(ir,idime)
        xnerg_chm(ir,idime) = xnerg_chm( 1,idime)
     end do

     ir = ir - 1

     if ( ir == 1 ) then      
        lnerg_chm(1) = icrit
        do idime = 1,ndime
           xnerg_chm(1,idime) = raux(idime)
        end do
        return
     end if
  end if

  ii = len
  jj = len + len

200 if ( jj <= ir ) then
     if ( jj < ir ) then
        if ( xnerg_chm(jj,1) > xnerg_chm(jj+1,1) ) then
           jj = jj + 1
        end if
     end if

     if ( raux(1) > xnerg_chm(jj,1) ) then
        lnerg_chm(ii) = lnerg_chm(jj) 
        do idime = 1,ndime
           xnerg_chm(ii,idime) = xnerg_chm(jj,idime)
        end do
        ii = jj
        jj = jj + jj
     else
        jj = ir + 1
     end if

     goto 200
  end if

  lnerg_chm(ii) = icrit
  do idime = 1,ndime
     xnerg_chm(ii,idime) = raux(idime)
  end do
  goto 100

end subroutine chm_sorten
