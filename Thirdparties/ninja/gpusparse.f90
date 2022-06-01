module gpusparse

  interface gpubsrmv
     module procedure gpudoublebsrmv,gpusinglebsrmv
  end interface gpubsrmv

  interface gpuwarpcsrmv
     module procedure gpudoublewarpcsrmv,gpusinglewarpcsrmv
  end interface gpuwarpcsrmv

  interface gpuvectorcsrmv
     module procedure gpudoublevectorcsrmv,gpusinglevectorcsrmv
  end interface gpuvectorcsrmv
  
contains

  subroutine gpudoublebsrmv(hand,ND,MD,NNZB,alp,desc,d_val,d_row,d_col,V,d_x,bet,d_y,str)
    implicit none
    integer*8    :: hand,str,desc
    integer*4    :: ND,MD,NNZB,V
    real*8       :: alp,bet
    integer*8    :: d_val,d_row,d_col,d_x,d_y
    
  end subroutine gpudoublebsrmv

  subroutine gpusinglebsrmv(hand,ND,MD,NNZB,alp,desc,d_val,d_row,d_col,V,d_x,bet,d_y,str)
    implicit none
    integer*8    :: hand,str,desc
    integer*4    :: ND,MD,NNZB,V
    real*4       :: alp,bet
    integer*8    :: d_val,d_row,d_col,d_x,d_y
    
  end subroutine gpusinglebsrmv
  
  subroutine gpudoublewarpcsrmv(hand,ND,MD,NNZB,alp,desc,d_val,d_row,d_col,V,d_x,bet,d_y,str)
    implicit none
    integer*8    :: hand,str,desc
    integer*4    :: ND,MD,NNZB,V
    real*8       :: alp,bet
    integer*8    :: d_val,d_row,d_col,d_x,d_y
    
  end subroutine gpudoublewarpcsrmv

  subroutine gpusinglewarpcsrmv(hand,ND,MD,NNZB,alp,desc,d_val,d_row,d_col,V,d_x,bet,d_y,str)
    implicit none
    integer*8    :: hand,str,desc
    integer*4    :: ND,MD,NNZB,V
    real*4       :: alp,bet
    integer*8    :: d_val,d_row,d_col,d_x,d_y

  end subroutine gpusinglewarpcsrmv

    subroutine gpudoublevectorcsrmv(hand,ND,MD,NNZB,alp,desc,d_val,d_row,d_col,V,d_x,bet,d_y,str)
    implicit none
    integer*8    :: hand,str,desc
    integer*4    :: ND,MD,NNZB,V
    real*8       :: alp,bet
    integer*8    :: d_val,d_row,d_col,d_x,d_y
    
  end subroutine gpudoublevectorcsrmv

  subroutine gpusinglevectorcsrmv(hand,ND,MD,NNZB,alp,desc,d_val,d_row,d_col,V,d_x,bet,d_y,str)
    implicit none
    integer*8    :: hand,str,desc
    integer*4    :: ND,MD,NNZB,V
    real*4       :: alp,bet
    integer*8    :: d_val,d_row,d_col,d_x,d_y

  end subroutine gpusinglevectorcsrmv

end module gpusparse
