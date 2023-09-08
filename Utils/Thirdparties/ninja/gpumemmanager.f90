MODULE GPUMEM
  use def_master, only : IMASTER, INOTSLAVE,npoi1
  use mod_communications, only : PAR_MIN,PAR_MAX,PAR_INTERFACE_NODE_EXCHANGE,PAR_INTERFACE_NODE_EXCHANGE_VALUE,COMMUNICATOR_VALUE
  use def_master, only : kfl_paral
  use def_parall, only : kfl_cores_per_gpu
  use mod_parall, only : PAR_WORLD_SIZE
  use def_elmtyp
  use def_kintyp
  use iso_c_binding
 IMPLICIT NONE
  integer*8 :: d_start,d_stop,hansp,desc,desccsr
  integer*8 :: d_colcsr,d_valcsr,d_rowcsr,d_row,d_rc,d_hybmat,d_rcbdr
  integer*8 :: d_memst, d_memallocst,d_memcurr, d_memend,bdrstream,bdrhandle,d_tempbuff=0
  integer*4 :: alignm,id,offset,gpuend,cpustart,no_sms,csrband
  integer*8 :: d_totbytes,d_totfreebytes,d_bytesleft,d_bytesmin,tsize
  real*8 :: temp,time,allocper,bandratio
  INTEGER :: gpumemflag = 0  
  INTEGER :: gpurdmaflag = 0
  INTEGER :: gpurdmaassflag = 0
  INTEGER :: lightflag = 0
  INTEGER :: spmvformatflag = 0
  INTEGER :: hybridflag = 0
  INTEGER :: hostmemflag = 0
  
  real*8,pointer :: ps1(:),ps2(:),tempbuff(:)
  type(c_ptr) :: cptr

  integer*8 :: d_send,d_recv,d_perm
  integer                                   :: bound_dim
  integer,                       pointer    :: bound_perm(:)
  integer,                       pointer    :: bound_size(:)
END MODULE GPUMEM

module gpumanager

  interface gpumemassign
     module procedure gpumemassign_ip, gpumemassign_2ip, gpumemassign_perm
  end interface gpumemassign
  interface gpumemunassign
     module procedure  gpumemunassign_ip, gpumemunassign_2ip
  end interface gpumemunassign
  interface memcpytogpu
     module procedure memcpytogpu_dp, memcpytogpu_2dp, memcpytogpu_ip,memcpytogpu_2ip
     module procedure memcpytogpu_dp_2d, memcpytogpu_2dp_2d, memcpytogpu_ip_2d,memcpytogpu_2ip_2d
  end interface memcpytogpu
  
  interface memcpyfromgpu
     module procedure memcpyfromgpu_dp,memcpyfromgpu_2dp,memcpyfromgpu_ip,memcpyfromgpu_2ip
     module procedure memcpyfromgpu_dp_2d,memcpyfromgpu_2dp_2d,memcpyfromgpu_ip_2d,memcpyfromgpu_2ip_2d
  end interface memcpyfromgpu
  
  interface memcpytogpuasync
     module procedure memcpytogpuasync_dp, memcpytogpuasync_2dp, memcpytogpuasync_ip,&
          memcpytogpuasync_2ip,memcpytogpuasync_pinned
  end interface memcpytogpuasync
  
  interface memcpyfromgpuasync
     module procedure memcpyfromgpuasync_dp,memcpyfromgpuasync_2dp,memcpyfromgpuasync_ip,&
          memcpyfromgpuasync_2ip,memcpyfromgpuasync_pinned
  end interface memcpyfromgpuasync
  
  interface par_interface_node_exchange_gpu
     module procedure par_interface_node_exchange_gpu_nhyb,par_interface_node_exchange_gpu_hyb
  end interface par_interface_node_exchange_gpu

  interface gpumeminitiate
!    module procedure gpumeminitiate_param,gpumeminitiate_noparam, &
    module procedure gpumeminitiate_param_sfc,gpumeminitiate_noparam_sfc
  end interface gpumeminitiate

  interface cpumemregister
     module procedure cpumemregister_rp,cpumemregister_ip
  end interface cpumemregister

  interface cpumemunregister
     module procedure cpumemunregister_rp,cpumemunregister_ip
  end interface cpumemunregister
  
contains
  
  subroutine par_interface_node_exchange_gpu_nhyb(stream,V,exlen,d_xx,what,where,wsynch,dom_k)
    use GPUMEM
    implicit none
    integer*8 :: d_xx,stream
    integer*4 :: V,exlen
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: where
    character(*),         optional, intent(in)    :: wsynch
    integer,              optional, intent(in)    :: dom_k
    real*8                                        :: alp = 1.0
    
    if (gpurdmaflag == 1) then
       !
       ! GPU direct RDMA
       !
       if( INOTSLAVE ) return
              
       if (gpurdmaassflag == 0) then
          call communicator_value(bound_dim,bound_perm,bound_size,where)
          tsize = bound_dim
          tsize = tsize*V*8
          
          call gpumemassign(d_send,tsize)
          call gpumemassign(d_recv,tsize)
          tsize = bound_dim
          tsize = tsize*4 
          call gpumemassign(d_perm,tsize)

          call memcpytogpu(d_perm,bound_perm,bound_dim*4,'perm array')
          gpurdmaassflag = 1
       end if
       
       call initialize_boundaries(npoi1,V,bound_dim,d_perm,d_send,d_recv,d_xx)
       
       call int2cptr(d_send,cptr)
       call c_f_pointer(cptr,ps1,[bound_dim*V])

       call int2cptr(d_recv,cptr)
       call c_f_pointer(cptr,ps2,[bound_dim*V])
       
       call PAR_INTERFACE_NODE_EXCHANGE_VALUE(V,exlen,ps1,ps2,what,where,wsynch,dom_k)     
       
       if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' .or. trim(what) == 'SUMI' ) then
                    
          call permute_sum(npoi1,V,bound_dim,d_perm,d_send,d_recv,d_xx)          
          
       else
          call runend('PAR_INTERFACE_NODE_EXCHANGE_GPU: Undefined operation only SUM supported')
       endif
       
    else
       !
       ! Boundary is sent back through PCI transfer
       !
       if( INOTSLAVE ) return

       tsize = exlen
       tsize = tsize*8
       if ( hostmemflag == 0 .and. exlen /= 0) THEN
          call gpumallochost(d_tempbuff,tsize)
          call int2cptr(d_tempbuff,cptr)
          call c_f_pointer(cptr,tempbuff,[exlen])
          hostmemflag = 1
       end if
       
       if ( stream == 0 ) then
          !
          ! Synchronous
          !
          call MEMCPYFROMGPU(tempbuff,d_xx,tsize,"bdr ex")
          call PAR_INTERFACE_NODE_EXCHANGE(V,tempbuff,what,where,wsynch)     
          call MEMCPYTOGPU(d_xx,tempbuff,tsize,"bdr ex")

       else
          !
          ! Asynchronous
          !
          call MEMCPYFROMGPUASYNC(tempbuff,d_xx,tsize,stream,"bdr ex")
          call GPUSTREAMSYNC(stream)
          call PAR_INTERFACE_NODE_EXCHANGE(V,tempbuff,what,where,wsynch)     
          call MEMCPYTOGPUASYNC(d_xx,tempbuff,tsize,stream,"bdr ex")       
          
       end if
    end if

  end subroutine par_interface_node_exchange_gpu_nhyb

  subroutine par_interface_node_exchange_gpu_hyb(V,xx,what,where,wsynch,dom_k)
    use GPUMEM
    use def_master,only                           : npoi1
    implicit none
    real*8                                        ::  xx(*)
    integer*4 :: V,exlen
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: where
    character(*),         optional, intent(in)    :: wsynch
    integer,              optional, intent(in)    :: dom_k

    call PAR_INTERFACE_NODE_EXCHANGE(V,xx,what,where,wsynch)
    
  end subroutine par_interface_node_exchange_gpu_hyb
   
  subroutine gpumeminitiate_param(ND,V)
    use GPUMEM
    implicit none
    integer*4    ::  ND,V
    
    if (gpumemflag == 0) then
       call gpuset(kfl_paral)
       !call gpureset()
       call gpumemgetinfo(d_totfreebytes,d_totbytes)
       open(unit = 4,file = "GPUconfig.dat")
       read(4,*)allocper,gpurdmaflag,lightflag
       close(4)
       
       if ( gpurdmaflag == hybridflag .and. hybridflag /= 0) then
          call runend('Cannot enable GPU Direct RDMA and Hybrid computation at same time. !!!CHECK CONFIG!!!')
       end if
       if(kfl_paral == 0) allocper = 0.1
       d_totfreebytes = (allocper * d_totbytes)

       temp = d_totfreebytes
       call getgpualignment(alignm)
       call gpugetsms(no_sms)
       
       call gpumalloc(d_memallocst,d_totfreebytes,"manager alloc")
       
       d_bytesmin = d_totfreebytes
       d_bytesleft = d_totfreebytes
       d_memst = d_memallocst
       d_memcurr = d_memst
       d_memend = d_memst + d_totfreebytes
       gpumemflag = 1

       if( INOTSLAVE )  then
          open(unit = 5,file = "GPUstats.txt")
          write (5,*)"GPUMEM MANAGER Allocated gigabytes ---->>>> ",(temp/(1024*1024*1024))
          close(5)
       end if
       call GPUSTREAMCREATE(bdrstream)
       call GPUEVENTCREATE(d_start)
       call GPUEVENTCREATE(d_stop)

    end if
    call calculateoffset(ND,V)
    call GPUEVENTRECORD(d_start)
  end subroutine gpumeminitiate_param

  subroutine gpumeminitiate_noparam()
    use GPUMEM
    implicit none
    if (gpumemflag == 0) then
       call gpuset(kfl_paral)
       !call gpureset()
       call gpumemgetinfo(d_totfreebytes,d_totbytes)
       open(unit = 4,file = "GPUconfig.dat")
       read(4,*)allocper,gpurdmaflag,lightflag
       close(4)
       
       if ( gpurdmaflag == hybridflag .and. hybridflag /= 0) then
          call runend('Cannot enable GPU Direct RDMA and Hybrid computation at same time. !!!CHECK CONFIG!!!')
       end if  
     
       d_totfreebytes = (allocper * d_totbytes)


       temp = d_totfreebytes
       call getgpualignment(alignm)
       call gpugetsms(no_sms)
       
       call gpumalloc(d_memallocst,d_totfreebytes,"manager alloc")
       
       d_bytesmin = d_totfreebytes
       d_bytesleft = d_totfreebytes
       d_memst = d_memallocst
       d_memcurr = d_memst
       d_memend = d_memst + d_totfreebytes
       gpumemflag = 1

       if( INOTSLAVE )  then
          open(unit = 5,file = "GPUstats.txt")
          write (5,*)"GPUMEM MANAGER Allocated gigabytes ---->>>> ",(temp/(1024*1024*1024))
          close(5)
       end if
       call GPUSTREAMCREATE(bdrstream)
       call GPUEVENTCREATE(d_start)
       call GPUEVENTCREATE(d_stop)

    end if
    call GPUEVENTRECORD(d_start)
  end subroutine gpumeminitiate_noparam

  subroutine gpumeminitiate_param_sfc(ND,V)
    use GPUMEM

   implicit none
    integer*4    ::  ND,V
    real*8, allocatable   :: lcorr(:)      
    integer(ip)             :: ii, ipart
    logical(lg)             :: auxl 
    integer(ip)             :: npart
    integer                 :: ncpucore,mpierr

    npart=PAR_WORLD_SIZE-1


    if (gpumemflag == 0) then
       if(kfl_cores_per_gpu==0) then
        call gpuset(kfl_paral)
       else
        call gpusetmany(kfl_paral,kfl_cores_per_gpu)
       end if



       !call gpureset()
       call gpumemgetinfo(d_totfreebytes,d_totbytes)

       open(unit = 4,file = "GPUconfig.dat")
       read(4,*)allocper,gpurdmaflag,lightflag
       close(4)


      !To calculate the variable allocation based on rank-element    
      allocate(lcorr(npart))
      inquire( file="rank-elements.dat", exist=auxl )

         
         if(auxl) then
                   open(unit=1405875, action='READ' ,file='rank-elements.dat')
            do ii=1, npart
               read(1405875,*) ipart, lcorr(ii)
            end do
            lcorr(npart) = npart - sum(lcorr(1:npart-1))
            close(1405875)
         else
            lcorr(1:npart) = 1_rp
         endif
   
       if ( gpurdmaflag == hybridflag .and. hybridflag /= 0) then
          call runend('Cannot enable GPU Direct RDMA and Hybrid computation at same time. !!!CHECK CONFIG!!!')
       end if

       if(kfl_cores_per_gpu==0) then
           d_totfreebytes = (allocper * d_totbytes)
       else
            
           if(kfl_paral == 0) then
              d_totfreebytes = (0.08 * d_totbytes)
           else
             !combines GPUconfig.dat with rank-element 
              d_totfreebytes = lcorr(kfl_paral)*((allocper * d_totbytes)/kfl_cores_per_gpu)
           end if               
 
      end if

       temp = d_totfreebytes
       call getgpualignment(alignm)
       call gpugetsms(no_sms)
       
       call gpumalloc(d_memallocst,d_totfreebytes,"manager alloc")
       
       d_bytesmin = d_totfreebytes
       d_bytesleft = d_totfreebytes
       d_memst = d_memallocst
       d_memcurr = d_memst
       d_memend = d_memst + d_totfreebytes
       gpumemflag = 1

       if( INOTSLAVE )  then
          open(unit = 5,file = "GPUstats.txt")
          write (5,*)"GPUMEM MANAGER Allocated gigabytes ---->>>> ",(temp/(1024*1024*1024))
          close(5)
       end if
       call GPUSTREAMCREATE(bdrstream)
       call GPUEVENTCREATE(d_start)
       call GPUEVENTCREATE(d_stop)

    end if
    call calculateoffset(ND,V)
    call GPUEVENTRECORD(d_start)
  end subroutine gpumeminitiate_param_sfc

  subroutine gpumeminitiate_noparam_sfc()
    use GPUMEM
    implicit none
    real*8, allocatable       :: lcorr(:)      
    integer             :: ii, ipart
    logical             :: auxl 
    integer             :: npart
    
   npart=PAR_WORLD_SIZE-1
   if (gpumemflag == 0) then
       if(kfl_cores_per_gpu==0) then
        call gpuset(kfl_paral)
       else
        call gpusetmany(kfl_paral,kfl_cores_per_gpu)
       end if

       call gpumemgetinfo(d_totfreebytes,d_totbytes)
       open(unit = 4,file = "GPUconfig.dat")
       read(4,*)allocper,gpurdmaflag,lightflag
       close(4)
       
       if ( gpurdmaflag == hybridflag .and. hybridflag /= 0) then
          call runend('Cannot enable GPU Direct RDMA and Hybrid computation at same time. !!!CHECK CONFIG!!!')
       end if  
 

      !To calculate the variable allocation based on rank-element    
      allocate(lcorr(npart))

      inquire( file="rank-elements.dat", exist=auxl )

      if(auxl) then
         open(unit=1405873, action='READ', file='rank-elements.dat')
         do ii=1, npart
            read(1405873,*) ipart, lcorr(ii)
         end do
         lcorr(npart) = npart - sum(lcorr(1:npart-1))
      else
         lcorr(1:npart) = 1_rp
      endif

       if(kfl_cores_per_gpu==0) then
                d_totfreebytes = (allocper * d_totbytes)
       else
            
               if(kfl_paral == 0) then
                 d_totfreebytes = (0.1 * d_totbytes)
                else
                 !combines GPUconfig.dat with rank-element 
                d_totfreebytes = lcorr(kfl_paral)*((allocper * d_totbytes)/kfl_cores_per_gpu)
                end if               
 
      end if



       temp = d_totfreebytes
       call getgpualignment(alignm)
       call gpugetsms(no_sms)
       
       call gpumalloc(d_memallocst,d_totfreebytes,"manager alloc")
       
       d_bytesmin = d_totfreebytes
       d_bytesleft = d_totfreebytes
       d_memst = d_memallocst
       d_memcurr = d_memst
       d_memend = d_memst + d_totfreebytes
       gpumemflag = 1

       if( INOTSLAVE )  then
          open(unit = 5,file = "GPUstats.txt")
          write (5,*)"GPUMEM MANAGER Allocated gigabytes ---->>>> ",(temp/(1024*1024*1024))
          close(5)
       end if
       call GPUSTREAMCREATE(bdrstream)
       call GPUEVENTCREATE(d_start)
       call GPUEVENTCREATE(d_stop)

    end if
    call GPUEVENTRECORD(d_start)
  end subroutine gpumeminitiate_noparam_sfc

  subroutine calculateoffset(ND,V)
    use GPUMEM
    implicit none
    integer*4    :: ND,V
   
    gpuend = (ND*(100-offset))/100
    cpustart = gpuend + 1

    if( gpuend > npoi1 ) gpuend = npoi1
    
  end subroutine calculateoffset
  
  subroutine gpumemassign_ip(poin,sss)
    use GPUMEM
    implicit none
    integer*8 :: poin
    integer*4 :: sss,md
    poin = d_memcurr
    md = mod(sss,alignm)
    !!----GPU address alignment----------
    if (md /= 0) then
       sss = sss + (alignm-md)
    end if
    d_bytesleft = d_bytesleft - sss
    temp = d_bytesleft
    
    d_memcurr = d_memcurr + sss 
    if (d_memcurr >= d_memend) then
       print *,'GPU ran out of memory'
       stop
    end if
    
  end subroutine gpumemassign_ip

  subroutine gpumemassign_2ip(poin,sss)
    use GPUMEM
    implicit none
    integer*8 :: poin,sss
    integer*4 :: md,somesss
    poin = d_memcurr
    somesss = sss
    md = mod(somesss,alignm)
    !!----GPU address alignment----------
    if (md /= 0) then
       sss = sss + (alignm-md)
    end if
    d_bytesleft = d_bytesleft - sss
    temp = d_bytesleft
    
    d_memcurr = d_memcurr + sss 
    if (d_memcurr >= d_memend) then
       print *,'GPU ran out of memory'
       stop
    end if
    
  end subroutine gpumemassign_2ip
  
  subroutine gpumemassign_perm(poin,sss,flag)
    use GPUMEM
    implicit none
    integer*8 :: poin,sss
    integer*4 :: md,flag,somesss

    if(d_memst .ne. d_memcurr) then
       call runend('permanent gpu memory routine cannot be called in middle of pool assignments. Call immediately after gpuinit')
    end if
    somesss = sss
    poin = d_memcurr
    md = mod(somesss,alignm)
    !!----GPU address alignment----------
    if (md /= 0) then
       sss = sss + (alignm-md)
    end if
    d_bytesleft = d_bytesleft - sss
    temp = d_bytesleft
    
    d_memcurr = d_memcurr + sss

    if (flag == 1) d_memst = d_memcurr
    
    if (d_memcurr >= d_memend) then
       print *,'GPU ran out of memory'
       stop
    end if

  end subroutine gpumemassign_perm
  
  subroutine gpumemreset()
    use GPUMEM

    implicit none
    call GPUEVENTRECORD(d_stop)
    call PAR_MIN(d_bytesleft,'IN MY CODE')
    call GPUSYNC()

    if(d_tempbuff /= 0 )    call gpufreehost(d_tempbuff)
    call GPUEVENTELAPSEDTIME(time,d_start,d_stop)
    call PAR_MAX(time,'IN MY CODE')

    if ( INOTSLAVE ) then
       temp = d_bytesleft
       open(5, file="GPUstats.txt", status="old", position="append", action="write")  
       write (5,*)"GPUMEM MANAGER Reseting, Unused Memory (GBs)---->>>> ",(temp/(1024*1024*1024))
       write (5,*)"Solver Timing in milliseconds,              ---->>>> ",time
       close(5)
    end if
    d_memcurr = d_memst
    d_bytesleft = d_totfreebytes
    gpurdmaassflag = 0
    spmvformatflag = 0
    hostmemflag = 0
          
  end subroutine gpumemreset



  subroutine memcpytogpu_dp(dst,src,sz,str)
    implicit none
    integer*8 :: dst
    real*8,intent(in) :: src(*)
    integer*4 :: sz
    character(len=*),intent(in) :: str
    call memcpytogpuapi(dst,src,sz,str)
    
  end subroutine memcpytogpu_dp
  
  subroutine memcpytogpu_2dp(dst,src,sz,str)
    implicit none
    integer*8 :: dst
    real*8,intent(in) :: src(*)
    integer*8 :: sz
    character(len=*),intent(in) :: str
    call memcpytogpuapi2(dst,src,sz,str)
    
  end subroutine memcpytogpu_2dp

  subroutine memcpytogpu_ip(dst,src,sz,str)
    implicit none
    integer*8 :: dst
    integer*4,intent(in) :: src(*)
    integer*4 :: sz
    character(len=*),intent(in) :: str
    call memcpytogpuapi(dst,src,sz,str)
    
  end subroutine memcpytogpu_ip
  
  subroutine memcpytogpu_2ip(dst,src,sz,str)
    implicit none
    integer*8 :: dst
    integer*4,intent(in) :: src(*)
    integer*8 :: sz
    character(len=*),intent(in) :: str
    call memcpytogpuapi2(dst,src,sz,str)
    
  end subroutine memcpytogpu_2ip


    subroutine memcpytogpu_dp_2d(dst,src,sz,str)
    implicit none
    integer*8 :: dst
    integer*4 :: sz
    real*8,intent(in) :: src((sz/8),1)
    character(len=*),intent(in) :: str
    call memcpytogpuapi(dst,src,sz,str)
    
  end subroutine memcpytogpu_dp_2d
  
  subroutine memcpytogpu_2dp_2d(dst,src,sz,str)
    implicit none
    integer*8 :: dst
    integer*8 :: sz
    real*8,intent(in) :: src((sz/8),1)
    character(len=*),intent(in) :: str
    call memcpytogpuapi2(dst,src,sz,str)
    
  end subroutine memcpytogpu_2dp_2d

  subroutine memcpytogpu_ip_2d(dst,src,sz,str)
    implicit none
    integer*8 :: dst
    integer*4 :: sz
    integer*4,intent(in) :: src((sz/4),1)
    character(len=*),intent(in) :: str
    call memcpytogpuapi(dst,src,sz,str)
    
  end subroutine memcpytogpu_ip_2d
  
  subroutine memcpytogpu_2ip_2d(dst,src,sz,str)
    implicit none
    integer*8 :: dst
    integer*8 :: sz
    integer*4,intent(in) :: src((sz/4),1)
    character(len=*),intent(in) :: str
    call memcpytogpuapi2(dst,src,sz,str)
    
  end subroutine memcpytogpu_2ip_2d

  
  subroutine memcpyfromgpu_dp(dst,src,sz,str)
    implicit none
    integer*8 :: src
    real*8:: dst(*)
    integer*4 :: sz
    character(len=*),intent(in) :: str
    call memcpyfromgpuapi(dst,src,sz,str)
    
  end subroutine memcpyfromgpu_dp
  
  subroutine memcpyfromgpu_2dp(dst,src,sz,str)
    implicit none
    integer*8 :: src
    real*8 :: dst(*)
    integer*8 :: sz
    character(len=*),intent(in) :: str
    call memcpyfromgpuapi2(dst,src,sz,str)
    
  end subroutine memcpyfromgpu_2dp

  subroutine memcpyfromgpu_ip(dst,src,sz,str)
    implicit none
    integer*8 :: src
    integer*4:: dst(*)
    integer*4 :: sz
    character(len=*),intent(in) :: str
    call memcpyfromgpuapi(dst,src,sz,str)
    
  end subroutine memcpyfromgpu_ip
  
  subroutine memcpyfromgpu_2ip(dst,src,sz,str)
    implicit none
    integer*8 :: src
    integer*4 :: dst(*)
    integer*8 :: sz
    character(len=*),intent(in) :: str
    call memcpyfromgpuapi2(dst,src,sz,str)

  end subroutine memcpyfromgpu_2ip

  subroutine memcpyfromgpu_dp_2d(dst,src,sz,str)
    implicit none
    integer*8 :: src
    integer*4 :: sz
    real*8:: dst((sz/8),1)
    character(len=*),intent(in) :: str
    call memcpyfromgpuapi(dst,src,sz,str)
    
  end subroutine memcpyfromgpu_dp_2d
  
  subroutine memcpyfromgpu_2dp_2d(dst,src,sz,str)
    implicit none
    integer*8 :: src
    integer*8 :: sz
    real*8 :: dst((sz/8),1)
    character(len=*),intent(in) :: str
    call memcpyfromgpuapi2(dst,src,sz,str)
    
  end subroutine memcpyfromgpu_2dp_2d

  subroutine memcpyfromgpu_ip_2d(dst,src,sz,str)
    implicit none
    integer*8 :: src
    integer*4 :: sz
    integer*4:: dst((sz/4),1)
    character(len=*),intent(in) :: str
    call memcpyfromgpuapi(dst,src,sz,str)
    
  end subroutine memcpyfromgpu_ip_2d
  
  subroutine memcpyfromgpu_2ip_2d(dst,src,sz,str)
    implicit none
    integer*8 :: src
    integer*8 :: sz
    integer*4 :: dst((sz/4),1)
    character(len=*),intent(in) :: str
    call memcpyfromgpuapi2(dst,src,sz,str)

  end subroutine memcpyfromgpu_2ip_2d

  subroutine memcpytogpuasync_dp(dst,src,sz,stream,str)
    implicit none
    integer*8 :: dst,stream
    real*8 :: src(*)
    integer*4 :: sz
    character(len=*),intent(in) :: str
    call memcpytogpuasyncapi(dst,src,sz,stream,str)
    
  end subroutine memcpytogpuasync_dp
  
  subroutine memcpytogpuasync_2dp(dst,src,sz,stream,str)
    implicit none
    integer*8 :: dst,stream
    real*8 :: src(*)
    integer*8 :: sz
    character(len=*),intent(in) :: str
    call memcpytogpuasyncapi2(dst,src,sz,stream,str)
    
  end subroutine memcpytogpuasync_2dp

  subroutine memcpytogpuasync_ip(dst,src,sz,stream,str)
    implicit none
    integer*8 :: dst,stream
    integer*4 :: src(*)
    integer*4 :: sz
    character(len=*),intent(in) :: str
    call memcpytogpuasyncapi(dst,src,sz,stream,str)
    
  end subroutine memcpytogpuasync_ip
  
  subroutine memcpytogpuasync_2ip(dst,src,sz,stream,str)
    implicit none
    integer*8 :: dst,stream
    integer*4 :: src(*)
    integer*8 :: sz
    character(len=*),intent(in) :: str
    call memcpytogpuasyncapi2(dst,src,sz,stream,str)
    
  end subroutine memcpytogpuasync_2ip

  subroutine memcpytogpuasync_pinned(dst,src,sz,stream,str)
    implicit none
    integer*8 :: dst,stream,src
    integer*8 :: sz
    character(len=*),intent(in) :: str
    call memcpytogpuasyncapi3(dst,src,sz,stream,str)
    
  end subroutine memcpytogpuasync_pinned

  subroutine memcpyfromgpuasync_pinned(dst,src,sz,stream,str)
    implicit none
    integer*8 :: src,stream,dst
    integer*4 :: sz
    character(len=*),intent(in) :: str
    call memcpyfromgpuasyncapi3(dst,src,sz,stream,str)
    
  end subroutine memcpyfromgpuasync_pinned
  
  subroutine memcpyfromgpuasync_dp(dst,src,sz,stream,str)
    implicit none
    integer*8 :: src,stream
    real*8:: dst(*)
    integer*4 :: sz
    character(len=*),intent(in) :: str
    call memcpyfromgpuasyncapi(dst,src,sz,stream,str)
    
  end subroutine memcpyfromgpuasync_dp
  
  subroutine memcpyfromgpuasync_2dp(dst,src,sz,stream,str)
    implicit none
    integer*8 :: src,stream
    real*8 :: dst(*)
    integer*8 :: sz
    character(len=*),intent(in) :: str
    call memcpyfromgpuasyncapi2(dst,src,sz,stream,str)
    
  end subroutine memcpyfromgpuasync_2dp

  subroutine memcpyfromgpuasync_ip(dst,src,sz,stream,str)
    implicit none
    integer*8 :: src,stream
    integer*4:: dst(*)
    integer*4 :: sz
    character(len=*),intent(in) :: str
    call memcpyfromgpuasyncapi(dst,src,sz,stream,str)
    
  end subroutine memcpyfromgpuasync_ip
  
  subroutine memcpyfromgpuasync_2ip(dst,src,sz,stream,str)
    implicit none
    integer*8 :: src,stream
    integer*4 :: dst(*)
    integer*8 :: sz
    character(len=*),intent(in) :: str
    call memcpyfromgpuasyncapi2(dst,src,sz,stream,str)
    
  end subroutine memcpyfromgpuasync_2ip

  subroutine gpumemunassign_ip(poin,sss)
    use GPUMEM
    implicit none
    integer*8 :: poin
    integer*4 :: sss,md
    poin = 0
    md = mod(sss,alignm)
    !!----GPU address alignment----------
    if (md /= 0) then
       sss = sss + (alignm-md)
    end if
    d_bytesleft = d_bytesleft + sss
    temp = d_bytesleft
    
    d_memcurr = d_memcurr - sss 
      
  end subroutine gpumemunassign_ip

  subroutine gpumemunassign_2ip(poin,sss)
    use GPUMEM
    implicit none
    integer*8 :: poin,sss
    integer*4 :: md,somesss
    poin = 0
    somesss = sss
    md = mod(somesss,alignm)
    !!----GPU address alignment----------
    if (md /= 0) then
       sss = sss + (alignm-md)
    end if
    d_bytesleft = d_bytesleft + sss
    temp = d_bytesleft
    
    d_memcurr = d_memcurr - sss
    
  end subroutine gpumemunassign_2ip

  subroutine cpumemregister_ip(ptr,sz)
    use GPUMEM
    implicit none
    integer*4      :: ptr(*)
    integer*8      :: sz
    call gpuhostregister(ptr,sz)

  end subroutine cpumemregister_ip

  subroutine cpumemregister_rp(ptr,sz)
    use GPUMEM
    implicit none
    real*8         :: ptr(*)
    integer*8      :: sz
    call gpuhostregister(ptr,sz)

  end subroutine cpumemregister_rp

  subroutine cpumemunregister_ip(ptr)
    use GPUMEM
    implicit none
    integer*4      :: ptr(*)
    call gpuhostunregister(ptr)
    
  end subroutine cpumemunregister_ip

  subroutine cpumemunregister_rp(ptr)
    use GPUMEM
    implicit none
    real*8         :: ptr(*)
    call gpuhostunregister(ptr)

  end subroutine cpumemunregister_rp

end module gpumanager
