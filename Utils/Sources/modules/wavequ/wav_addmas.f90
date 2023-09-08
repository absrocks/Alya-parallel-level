subroutine wav_addmas(npoin,a,m,f,ja,ia) 
  !----------------------------------------------------------------------
  !****f* mathru/bcsrab
  ! NAME 
  !     bcsrab
  ! DESCRIPTION
  !     Multiply a non symmetric matrix stored in BCSR by a vector
  !     y = A b 
  ! INPUT
  !    NPOIN ... Number of equations
  !    NDOFN ..... Number of variables
  !    A ......... Matrix
  !    JA ........ List of elements
  !    IA ........ Pointer to list of elements
  !    B ......... Vector
  ! OUTPUT
  !    Y ......... result vector
  ! USES
  ! USED BY
  !***
  !----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip), intent(in)    :: npoin
  real(rp),    intent(in)    :: m(*),f
  real(rp),    intent(inout) :: a(*)
  integer(ip), intent(in)    :: ja(*),ia(*)
  integer(ip)                :: ii,jj,col

  do ii=1,npoin
     jj=ia(ii)-1
     do while(jj<ia(ii+1)-1)
        jj=jj+1
        col=ja(jj)
        if(col==ii) then
           a(jj)=a(jj)+f*m(ii)
           jj=ia(ii+1)
        end if
     end do
  end do

end subroutine wav_addmas
