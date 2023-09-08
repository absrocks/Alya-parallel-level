subroutine GPUDMV(hand,ND,MD,NNZB,V,alp,d_val,d_row,d_col,d_x,bet,d_y,str)
  use gpumem, only     : lightflag,desc,spmvformatflag,desccsr,bdrhandle,bdrstream
  use gpumem, only     : d_colcsr,d_rowcsr,d_valcsr,d_rc,d_hybmat,d_rcbdr,no_sms
  use gpumanager
  use def_master, only : npoi1
  implicit none
  integer*4            :: ND,NNZB,V,MD,exlen,off
  integer*8            :: hand,d_val,d_col,d_row,d_x,d_y,d_off,d_rowoff,d_yoff
  integer*8            :: tsize,str
  real*8               :: alp,bet

  if( lightflag == 0 ) then     
     !
     ! NVIDIA library
     !
     d_rowoff = d_row + npoi1*4
     d_yoff = d_y + V*npoi1*8
     call CUDADBSRMV(bdrhandle,0,0,(ND-npoi1),MD,NNZB,alp,desc,&
          d_val,d_rowoff,d_col,V,d_x,bet,d_yoff)
     call CUDADBSRMV(hand,0,0,npoi1,MD,NNZB,alp,desc,&
          d_val,d_row,d_col,V,d_x,bet,d_y)
     exlen = (ND-npoi1)*V
     off = npoi1 * V
     d_off = d_y + off*8
     call PAR_INTERFACE_NODE_EXCHANGE_GPU(bdrstream,V,exlen,d_off,'SUMI','IN MY CODE','SYNCHRONOUS')  
     call GPUSYNC()     

  else
     !
     ! Vishal's library
     !
     tsize = 1
     tsize = tsize*4
     call GPUMEMSET(d_rc,0,tsize)
     call GPUMEMSET(d_rcbdr,0,tsize)
     
     if( lightflag == 2 ) then
        !
        ! Good one
        !
        no_sms   = 64
        d_rowoff = d_rowcsr + V*npoi1*4
        d_yoff   = d_y + V*npoi1*8        
        call LIGHTSPMVWARP(d_rcbdr,(ND-npoi1)*V,d_rowoff,d_colcsr,d_valcsr,d_x,d_yoff,bdrstream,2)   ! Boundary
        call LIGHTSPMVWARP(d_rc,npoi1*V,d_rowcsr,d_colcsr,d_valcsr,d_x,d_y,str,no_sms)               ! Interior
        exlen = (ND-npoi1)*V
        off   = npoi1 * V
        d_off = d_y + off*8
        call PAR_INTERFACE_NODE_EXCHANGE_GPU(bdrstream,V,exlen,d_off,'SUMI','IN MY CODE','SYNCHRONOUS')
        call GPUSYNC()
        
     else if( lightflag == 1 ) then
        !
        ! Not better!
        !
        d_rowoff = d_rowcsr + V*npoi1*4
        d_yoff   = d_y + V*npoi1*8        
        call LIGHTSPMVVECTOR(d_rcbdr,(ND-npoi1)*V,d_rowoff,d_colcsr,d_valcsr,d_x,d_yoff,bdrstream,2) ! Boundary
        call LIGHTSPMVVECTOR(d_rc,npoi1*V,d_rowcsr,d_colcsr,d_valcsr,d_x,d_y,str,no_sms)             ! Interior
        exlen = (ND-npoi1)*V
        off   = npoi1 * V
        d_off = d_y + off*8
        call PAR_INTERFACE_NODE_EXCHANGE_GPU(bdrstream,V,exlen,d_off,'SUMI','IN MY CODE','SYNCHRONOUS')
        call GPUSYNC()
        
     else
        
        call runend("NO Sparse Matrix Vector routine defined. Use 0 Cusparse, 1 lightspmv vector,2 ligthspmv warp")
        
     end if   
  end if
  
end subroutine GPUDMV



subroutine preparecsr(hand,ND,MD,NNZB,V,d_val,d_row,d_col,str)
  use gpumem, only     : lightflag,desc,spmvformatflag,desccsr,bdrstream,bdrhandle
  use gpumem, only     : d_colcsr,d_rowcsr,d_valcsr,d_rc,d_hybmat,d_rcbdr,csrband,bandratio
  use gpumanager
  implicit none
  integer*4                :: ND,NNZB,V,MD,NDV
  integer*8                :: hand,d_val,d_col,d_row
  integer*8                :: d_max
  integer*8                :: tsize,str
  integer*4,allocatable    :: dummr(:)
  if(lightflag == 0) then
     if(spmvformatflag == 0) then
        call CUSPARSE_MAT_DESCR_CREATE(desc)
        call CUSPARSE_SET_MAT_TYPE(desc,0)
        call CUSPARSE_SET_MAT_BASE(desc,1)
        call CUSPARSEHANDLE(bdrhandle)
        call CUSPARSEATTACHSTREAM(bdrhandle,bdrstream)
        spmvformatflag = 1
     end if
     
  else
     
     if(spmvformatflag == 0 ) then
        
        tsize = 1
        tsize = tsize*4
        call GPUMEMASSIGN(d_rc,tsize)
        call GPUMEMASSIGN(d_rcbdr,tsize)
        tsize = NNZB
        tsize = tsize*V*V*4
        call GPUMEMASSIGN(d_colcsr,tsize)
        tsize = ND*V + 1
        tsize = tsize*4
        call GPUMEMASSIGN(d_rowcsr,tsize)
        
        call CUSPARSE_MAT_DESCR_CREATE(desccsr)
        call CUSPARSE_SET_MAT_TYPE(desccsr,0)
        call CUSPARSE_SET_MAT_BASE(desccsr,0)
        
        if ( V > 1 ) then
           tsize = NNZB
           tsize = tsize*V*V*8
           call GPUMEMASSIGN(d_valcsr,tsize)  
           
           call CUSPARSE_MAT_DESCR_CREATE(desc)
           call CUSPARSE_SET_MAT_TYPE(desc,0)
           call CUSPARSE_SET_MAT_BASE(desc,1)
           
           call CUDADBSR2CSR(hand,0,ND,MD,desc,d_val,d_row,&
                d_col,V,desccsr,d_valcsr,d_rowcsr,d_colcsr)
        else
           d_valcsr = d_val
           call tobasezero(d_row,d_col,d_rowcsr,d_colcsr,ND,MD,NNZB,str)
        end if
        
!        if(lightflag == 3) then
!           call CUSPARSE_MAT_HYB_CREATE(d_hybmat)
!           call CUDADCSR2HYB(hand,ND,MD,desccsr,d_valcsr,d_rowcsr,d_colcsr,d_hybmat,0,0)
!        endif
        
        spmvformatflag =1
     end if
     NDV = ND * V
     tsize = 4
     call GPUMEMASSIGN(d_max,tsize)
     allocate(dummr(1))
     call GPUMEMSET(d_max,0,tsize)     
     call find_bandwidth(d_rowcsr,d_colcsr,d_max,NDV)
     call MEMCPYFROMGPU(dummr,d_max,4,"bandwidth")
     csrband = dummr(1)     
     deallocate(dummr)
     call GPUMEMUNASSIGN(d_max,tsize)
     bandratio = (ND*V)/csrband
  end if

end subroutine preparecsr

subroutine generatediag(ND,V,d_diag)
  use gpumem, only     : d_colcsr,d_rowcsr,d_valcsr
  implicit none
  integer*4            :: ND,V
  integer*8            :: d_diag

  call DDIAGONALGEN(ND*V,d_valcsr,d_rowcsr,d_colcsr,d_diag,0)

end subroutine generatediag
