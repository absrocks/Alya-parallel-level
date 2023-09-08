subroutine gpuqrfactorization(ND,V,amatr,ia,ja,perm)
  use gpuvarqr
  use gpumanager
  implicit none
  integer*4          :: ia(*),ja(*),perm(*)
  real*8             :: amatr(*)
  integer*4          :: ND,V,ii,N
  integer*8          :: sz,str
 
  !
  ! RAS: Restricted Additive Schwarz
  !   
  
  N = NNZB * V * V
  nrows = ND
  nvar = V
  str = 0
  call gpumeminitiate(ND,V)         
  if (INOTMASTER) then
     NNZB = ia(ND+1) -1
     sz = N
     sz = sz*8
     call GPUMEMASSIGN(d_rasmat,sz)
     call MEMCPYTOGPU(d_rasmat,amatr,sz,"direct solver matrix")
     sz = NNZB
     sz = sz*4
     call GPUMEMASSIGN(d_rasja,sz)
     call MEMCPYTOGPU(d_rasja,ja,sz,"direct solver col")
     sz = ND+1
     sz = sz*4
     call GPUMEMASSIGN(d_rasia,sz)
     call MEMCPYTOGPU(d_rasia,ia,sz,"direct solver row")
     
     
     if (V == 1) then
        
        if (rasprepflag == 0 ) then
           
           call cusolversphandle(handsol)
           call cusparsehandle(handsp)
           call cusolverspattachstream(handsol,str)
           call CUSPARSE_MAT_DESCR_CREATE(matdesc)
           call CUSPARSE_SET_MAT_TYPE(matdesc,0)
           call CUSPARSE_SET_MAT_BASE(matdesc,1)
           
        end if
        
        call cusolverspcreatecsrqrinfo(solinfo)
        call cusolverspxcsrqranalysis(handsol,nrows,nrows,NNZB,matdesc,d_rasia,d_rasja,solinfo)
        call cusolverspdcsrqrbufferinfo(handsol,nrows,nrows,NNZB,matdesc,d_rasmat,d_rasia,d_rasja,solinfo,datain,workin)
        call GPUMEMASSIGN(qrbuff,workin)        
        call cusolverspdcsrqrsetup(handsol,nrows,nrows,NNZB,matdesc,d_rasmat,d_rasia,d_rasja,0,solinfo)
        call cusolverspdcsrqrfactor(handsol,nrows,nrows,NNZB,solinfo,qrbuff)
     else
        
        call runend("more than one degree of freedom in GPU QR Direct solver")
        
     end if

  end if
end subroutine gpuqrfactorization


subroutine gpuqrsolve(rhs,unk)
  use gpuvarqr
  use gpumanager
  implicit none
  real*8             :: rhs(*),unk(*)
  integer*8          :: tsize,in,out
  
  tsize = nrows*nvar
  tsize = tsize*8
  if (INOTMASTER) then
     call GPUMEMASSIGN(in,tsize)
     call GPUMEMASSIGN(out,tsize)
     
     call MEMCPYTOGPU(in,rhs,tsize,"Direct solver rhs")
     
     call cusolverspdcsrqrsolve(handsol,nrows,nrows,in,out,solinfo,qrbuff)
     
     call MEMCPYFROMGPU(unk,out,tsize,"Direct solver rhs")
  end if
     call gpumemreset()
end subroutine gpuqrsolve

subroutine gpuqrdestroy()
  use gpuvarqr
  implicit none
  call cusolverspdestroycsrqrinfo(solinfo)
  
end subroutine gpuqrdestroy
