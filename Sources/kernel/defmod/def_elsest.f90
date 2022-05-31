module def_elsest 

  use def_kintyp
  !type i1p
  !   integer(ip), pointer :: l(:)           ! pointer to allocate new memory dynam.
  !end type i1p
  type poiarr
     type(octbox), pointer     :: o
  end type poiarr

  type octbox
     integer(ip)               :: id          ! My global ID
     integer(ip)               :: level       ! Generation
     integer(ip)               :: npoinbox    ! Number of nodes
     integer(ip)               :: nelembox    ! Number of elements
     integer(ip)               :: childid     ! Child ID (1->4 or 1->8)
     integer(ip)               :: whoiam      ! Father or have nodes
     integer(ip),  pointer     :: nodes(:)    
     integer(ip),  pointer     :: elems(:)    ! List of elements
     real(rp)                  :: minc(3)     ! Min coordinates
     real(rp)                  :: maxc(3)     ! Max coordinates
     type(octbox), pointer     :: parent      ! Pointer to parent
     type(octbox), pointer     :: children(:) ! Pointer to children
  end type octbox
  type(octbox)   :: octbox_init = octbox(&
       0_ip,&                      ! id         
       0_ip,&                      ! level      
       0_ip,&                      ! npoinbox   
       0_ip,&                      ! nelembox   
       0_ip,&                      ! childid    
       0_ip,&                      ! whoiam     
       null(),&                    ! nodes(:)   
       null(),&                    ! elems(:)   
       (/0.0_rp,0.0_rp,0.0_rp/),&  ! minc(3)    
       (/0.0_rp,0.0_rp,0.0_rp/),&  ! maxc(3)    
       null(),&                    ! parent     
       null())                     ! children(:)
  !
  ! General data
  !
  integer(ip)                  :: iunit(3)
  integer(ip)                  :: nthre,ndime_els,nelem_els,npoin_els
  integer(ip)                  :: kfl_memor(3)
  integer(ip),     pointer     :: pelpo_els(:)
  integer(ip),     pointer     :: lelpo_els(:)
  real(rp),        pointer     :: coord_els(:,:)
  !
  ! Bin structures ---------------------------------------------------------
  !
  type bintype
     integer(ip)           :: nboxe
     integer(ip)           :: mboel
     integer(ip)           :: nboxx(3)
     integer(ip)           :: dataf
     integer(ip)           :: iallo       
     integer(ip)           :: iwhat
     integer(ip), pointer  :: lnobo(:)
     integer(ip), pointer  :: lboel(:)
     integer(ip), pointer  :: pboel(:)
     type(i1p),   pointer  :: tboel(:)
     integer(ip), pointer  :: memor(:,:) 
     integer(ip)           :: memax
     integer(ip), pointer  :: kstat(:,:) 
     integer(ip), pointer  :: ksear(:)    
     integer(ip), pointer  :: kfirs(:)    
     integer(ip), pointer  :: kseco(:)    
     real(rp)              :: delta(3)
     real(rp)              :: comin(3)
     real(rp)              :: comax(3)
     real(rp)              :: lmini
     real(rp)              :: lmaxi
     real(rp),    pointer  :: elcod(:,:,:)
     real(rp),    pointer  :: cputi(:,:)    
  end type bintype
  type octtype
     integer(ip)           :: mboel
     integer(ip)           :: iallo
     integer(ip)           :: iwhat
     integer(ip), pointer  :: lnobo(:)
     integer(ip), pointer  :: nbono(:)
     integer(ip), pointer  :: pbono(:)
     type(poiarr),pointer  :: current(:)
     type(octbox),pointer  :: tree_root
     integer(ip), pointer  :: memor(:,:) 
     integer(ip)           :: memax
     integer(ip), pointer  :: kstat(:,:) 
     integer(ip), pointer  :: ksear(:)    
     integer(ip), pointer  :: kfirs(:)    
     integer(ip), pointer  :: kseco(:)    
     integer(ip)           :: limit
     integer(ip)           :: divmax
     real(rp)              :: comin(3)
     real(rp)              :: comax(3)
     real(rp)              :: lmini
     real(rp)              :: lmaxi
     real(rp),    pointer  :: elcod(:,:,:)
     integer(ip), pointer  :: lboel(:)
     real(rp),    pointer  :: cputi(:,:)    
  end type octtype

  type(bintype),  pointer  :: bin_struc(:)          
  type(octtype),  pointer  :: oct_struc(:)  
        
  type(bintype),  parameter :: bin_struc_init=bintype(&
       0_ip,&                      ! nboxe         
       0_ip,&                      ! mboel         
       (/0_ip,0_ip,0_ip/),&        ! nboxx(3)      
       0_ip,&                      ! dataf         
       0_ip,&                      ! iallo         
       0_ip,&                      ! iwhat         
       null(),&                    ! lnobo(:)      
       null(),&                    ! lboel(:)      
       null(),&                    ! pboel(:)      
       null(),&                    ! tboel(:)      
       null(),&                    ! memor(:,:)    
       0_ip,&                      ! memax         
       null(),&                    ! kstat(:,:)    
       null(),&                    ! ksear(:)      
       null(),&                    ! kfirs(:)      
       null(),&                    ! kseco(:)      
       (/0.0_rp,0.0_rp,0.0_rp/),&  ! delta(3)      
       (/0.0_rp,0.0_rp,0.0_rp/),&  ! comin(3)      
       (/0.0_rp,0.0_rp,0.0_rp/),&  ! comax(3)      
       0.0_rp,&                    ! lmini         
       0.0_rp,&                    ! lmaxi         
       null(),&                    ! elcod(:,:,:)  
       null())                     ! cputi(:,:)    
  type(octtype), parameter :: oct_struc_init=octtype(&
       0_ip,&                      ! mboel         
       0_ip,&                      ! iallo               
       0_ip,&                      ! iwhat              
       null(),&                    ! lnobo(:)        
       null(),&                    ! nbono(:)        
       null(),&                    ! pbono(:)        
       null(),&                    ! current(:)      
       null(),&                    ! tree_root       
       null(),&                    ! memor(:,:)      
       0_ip,&                      ! memax                  
       null(),&                    ! kstat(:,:)      
       null(),&                    ! ksear(:)        
       null(),&                    ! kfirs(:)        
       null(),&                    ! kseco(:)        
       0_ip,&                      ! limit               
       0_ip,&                      ! divmax             
       0.0_rp,&                    ! comin(3)        
       0.0_rp,&                    ! comax(3)      
       0.0_rp,&                    ! lmini         
       0.0_rp,&                    ! lmaxi     
       null(),&                    ! elcod(:,:,:)    
       null(),&                    ! lboel(:)        
       null())                     ! cputi(:,:)      
      
        
  integer(ip),    pointer  :: nboxe
  integer(ip),    pointer  :: mboel
  integer(ip),    pointer  :: nboxx(:)
  integer(ip),    pointer  :: dataf
  integer(ip),    pointer  :: iallo
  integer(ip),    pointer  :: iwhat
  integer(ip),    pointer  :: lboel(:)
  integer(ip),    pointer  :: nbono(:)
  integer(ip),    pointer  :: pboel(:)
  type(i1p),      pointer  :: tboel(:)
  type(poiarr),   pointer  :: current(:)
  type(octbox),   pointer  :: tree_root
  integer(ip),    pointer  :: chbox(:,:)  
  integer(ip),    pointer  :: limit
  integer(ip),    pointer  :: divmax
  integer(ip),    pointer  :: chdel(:,:)  
  integer(ip),    pointer  :: nexbo(:,:)  
  integer(ip),    pointer  :: memor(:,:)  
  integer(ip),    pointer  :: memax
  integer(ip),    pointer  :: kstat(:,:)  
  integer(ip),    pointer  :: ksear(:)    
  integer(ip),    pointer  :: kfirs(:)    
  integer(ip),    pointer  :: kseco(:)    
  real(rp),       pointer  :: delta(:)
  real(rp),       pointer  :: comin(:)
  real(rp),       pointer  :: comax(:)
  real(rp)                 :: lmini
  real(rp)                 :: lmaxi
  real(rp),       pointer  :: elcod(:,:,:)
  real(rp),       pointer  :: cputi(:,:)   

  interface elsest_memchk
     module procedure memrp1,memrp2,memrp3,&   ! reals
          &           memrp4,&
          &           memip1,memip2,memip3,&   ! 4-byte integers 
          &           memi81,memi82,memi83,&   ! 8-byte integers 
          &           mei1p1,mei1p2,&          ! types
          &           memlg1                   ! logical
  end interface 

contains
 

  function elsest_intost(integ)
    !-------------------------------------
    !
    !  Convert an integer(ip) to a string
    !
    !-------------------------------------
    implicit none
    integer(ip)   :: integ
    character(20) :: elsest_intost
    character(20) :: intaux

    write(intaux,*) integ
    elsest_intost=adjustl(intaux)

  end function elsest_intost

  subroutine memrp1(itask,ithre,istat,cumem,vanam,vacal,varia)
    !
    ! Real(rp)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    real(rp)                    :: varia(:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0.0_rp
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if 
    call memsum(cumem,ithre,lbyts)

  end subroutine memrp1

  subroutine memrp2(itask,ithre,istat,cumem,vanam,vacal,varia)
    !
    ! Real(rp)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    real(rp)                    :: varia(:,:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0.0_rp
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memrp2

  subroutine memrp3(itask,ithre,istat,cumem,vanam,vacal,varia)
    !
    ! Real(rp)(:,:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    real(rp)                    :: varia(:,:,:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0.0_rp
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memrp3

  subroutine memrp4(itask,ithre,istat,cumem,vanam,vacal,varia)
    !
    ! Real(rp)(:,:,:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    real(rp)                    :: varia(:,:,:,:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0.0_rp
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memrp4

  subroutine memip1(itask,ithre,istat,cumem,vanam,vacal,varia)
    !
    ! Integer(4)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    integer(4)                  :: varia(:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memip1

  subroutine memip2(itask,ithre,istat,cumem,vanam,vacal,varia)
    !
    ! Integer(4)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    integer(4)                  :: varia(:,:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memip2

  subroutine memip3(itask,ithre,istat,cumem,vanam,vacal,varia)
    !
    ! Integer(4)(:,:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    integer(4)                  :: varia(:,:,:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memip3

  subroutine mei1p1(itask,ithre,istat,cumem,vanam,vacal,varia)
    !
    ! Type(i1p)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    type(i1p)                   :: varia(:)
    integer(ip)                 :: lbyts,isize,nsize
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          lbyts=nsize*ip
          do isize=1,nsize
             nullify(varia(isize)%l)
          end do
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*ip
    end if
    call memsum(cumem,ithre,lbyts)

  end subroutine mei1p1

  subroutine mei1p2(itask,ithre,istat,cumem,vanam,vacal,varia)
    !
    ! Type(i1p)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    type(i1p)                   :: varia(:,:)
    integer(ip)                 :: lbyts,nsize,isiz1,nsiz1,isiz2,nsiz2
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          nsiz1=size(varia,1)
          nsiz2=size(varia,2)
          lbyts=nsize*ip
          do isiz2=1,nsiz2
             do isiz1=1,nsiz1
                nullify(varia(isiz1,isiz2)%l)
             end do
          end do
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*ip
    end if
    call memsum(cumem,ithre,lbyts)

  end subroutine mei1p2

  subroutine memi81(itask,ithre,istat,cumem,vanam,vacal,varia)
    !
    ! Integer(8)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    integer(8)                  :: varia(:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memi81

  subroutine memi82(itask,ithre,istat,cumem,vanam,vacal,varia)
    !
    ! Integer(8)(:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    integer(8)                  :: varia(:,:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memi82

  subroutine memi83(itask,ithre,istat,cumem,vanam,vacal,varia)
    !
    ! Integer(8)(:,:,:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    integer(8)                  :: varia(:,:,:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=0
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memi83

  subroutine memlg1(itask,ithre,istat,cumem,vanam,vacal,varia)
    !
    ! Logical(lg)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(ip),  intent(in)    :: istat,itask,ithre
    integer(ip),  intent(inout) :: cumem
    logical(lg)                 :: varia(:)
    integer(ip)                 :: lbyts
    if(itask==0) then
       if(istat==0) then
          lbyts=size(varia)*kind(varia)
          varia=.false.
       else
          call elsest_memerr(itask,vanam,vacal,istat)
       end if
    else
       lbyts=-size(varia)*kind(varia)
    end if    
    call memsum(cumem,ithre,lbyts)

  end subroutine memlg1
  
  subroutine memsum(cumem,ithre,lbyts)
    implicit none
    integer(ip), intent(in)  :: ithre,lbyts
    integer(ip), intent(out) :: cumem
    integer(ip)              :: tomem,i
    
    cumem = cumem+lbyts
    tomem = 0
    do i=1,nthre
       tomem = tomem + memor(1,i)+memor(2,i)+memor(3,i)+memor(4,i)&
            &  +memor(5,i)+memor(6,i)+memor(7,i)+memor(8,i)&
            &  +memor(9,i)+memor(10,i)  
    end do
  !$OMP CRITICAL(chkmemax)
    memax = max(memax,tomem)
  !$OMP END CRITICAL(chkmemax)

  end subroutine memsum

end module def_elsest
