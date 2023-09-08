subroutine nsi_element_operations_cuda(VECTOR_DIM,pnode,pgaus,list_elements,time1)
  use gpumanager
  use nastinvar
  implicit none
  integer(ip), intent(in)          :: VECTOR_DIM                                       !< Number of nodes
  integer(ip), intent(in)          :: pnode                                            !< Number of nodes
  integer(ip), intent(in)          :: pgaus                                            !< Number of Gauss points
  integer(ip), intent(in)          :: list_elements(VECTOR_DIM)                        !< List of elements
  real(rp),    intent(inout)       :: time1(10)                                        ! Timings

end subroutine nsi_element_operations_cuda

subroutine nastincudainit(kfl_savda,coord,ltype,lnods,gravi_nsi,elmda_gpvol,elmda_gpcar)
  use gpumanager
  use nastinvar
  
  implicit none
  integer(ip)  :: kfl_savda
  integer*8    :: sz
  real(rp)     :: coord(*),gravi_nsi(*) !gravi_nsi(3)
  real(rp)     :: elmda_gpvol(*),elmda_gpcar(*)
  integer(ip)  :: ltype(*),lnods(*)
  

  call gpumeminitiate()
  
  if (INOTMASTER) then
     sz = NDIME*NPOIN
     sz = sz*8
     call GPUMEMASSIGN(dcoord,sz,1)
     call MEMCPYTOGPU(dcoord,coord,sz,"coord cudanastin")
     
     sz = NELEM
     sz = sz*4
     call GPUMEMASSIGN(dltype,sz,1)
     call MEMCPYTOGPU(dltype,ltype,sz,"ltype cudanastin")
     
     sz = NELEM*MNODE
     sz = sz*4
     call GPUMEMASSIGN(dlnods,sz,1)
     call MEMCPYTOGPU(dlnods,lnods,sz,"lnods cudanastin")
     
     sz = 3
     sz = sz*8
     call GPUMEMASSIGN(dgravi_nsi,sz,1)
     call MEMCPYTOGPU(dgravi_nsi,gravi_nsi,sz,"gravi_nsi cudanastin")
     
     if( kfl_savda == 0 ) then
        
     else
        
        sz = mgaus * nelem
        sz = sz * 8
        call GPUMEMASSIGN(delmda_gpvol,sz,1)
        call MEMCPYTOGPU(delmda_gpvol,elmda_gpvol,sz,"elmda_gpvol cudanastin")
        
        sz = mgaus * nelem
        sz = sz * ndime * mnode
        sz = sz * 8
        call GPUMEMASSIGN(delmda_gpcar,sz,1)
        call MEMCPYTOGPU(delmda_gpcar,elmda_gpvol,sz,"elmda_gpvol cudanastin")
        
     end if
  endif

end subroutine nastincudainit

subroutine cuda_nsi_elmope(dt_rho_nsi,mu_rho_nsi,rhsid,veloc,flag)
  use gpumanager
  use nastinvar
  implicit none
  integer*8    :: sz
  real(rp)     :: veloc(*), rhsid(*),dt_rho_nsi(*),mu_rho_nsi(*)
  integer*4    :: flag
  
  if (INOTMASTER) then

     if(flag == 0) then
        sz = npoin * ndime
        sz = sz * 8
        call GPUMEMASSIGN(drhsid,sz)
        call MEMCPYTOGPU(drhsid,rhsid,sz,"rhsid cudanastin")
        
        sz = npoin * ndime
        sz = sz * ncomp_nsi * 8
        call GPUMEMASSIGN(dveloc,sz)
        call MEMCPYTOGPU(dveloc,veloc,sz,"veloc cudanastin")
        
        sz = npoin
        sz = sz * 8
        call GPUMEMASSIGN(ddt_rho_nsi,sz)
        call GPUMEMASSIGN(dmu_rho_nsi,sz)
        call MEMCPYTOGPU(ddt_rho_nsi,dt_rho_nsi,sz,"dt_rho_nsi cudanastin")
        call MEMCPYTOGPU(dmu_rho_nsi,mu_rho_nsi,sz,"mu_rho_nsi cudanastin")
        
     else
        
        sz = npoin * ndime
        sz = sz * 8
        call MEMCPYFROMGPU(rhsid,drhsid,sz,"rhsid from cudanastin")
        
        sz = npoin
        sz = sz * 8
        call MEMCPYFROMGPU(dt_rho_nsi,ddt_rho_nsi,sz,"dt_rho_nsi from cudanastin")
        call MEMCPYFROMGPU(mu_rho_nsi,dmu_rho_nsi,sz,"mu_rho_nsi from cudanastin")
        call GPUMEMRESET()
     
     endif
     
  endif

end subroutine cuda_nsi_elmope
