subroutine chm_reasou()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_reasou
  ! NAME 
  !    chm_reasou
  ! DESCRIPTION
  !    This routine reads source file.
  !    The source units must be [Kg/s.m^3].
  ! USES
  !    chm_updunk
  ! USED BY
  !    chemic
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_memchk
  use mod_elmgeo, only  :  elmgeo_natural_coordinates 
  implicit none
!
  logical(lg)           :: go_on
  character(20)         :: wline  
  integer(ip)           :: iclas,ipoin,iline,isour,ielem,inode,idime,igaus
  integer(ip)           :: pelty,pnode,ptopo,ifoun
  integer(ip)           :: nsour,nclas,ndims_src(2),my_nsour
  real(rp)              :: time1,time2,dummr
  real(rp)              :: elcod(3,64),deriv(3,64),shapf(64),coglo(3),coloc(3)
  real(rp)              :: volum,xjacm(3,3),xjaci(3,3),gpcar(3,64),detjm
!
   integer(ip), allocatable :: ihave(:) 
   real(rp),    allocatable :: srtmp(:,:) 
!
!  Debug source option 
!
  if (kfl_sourc_chm == -1 ) then      
      !
      if( INOTMASTER ) then
      !
      !  Sequential or slave
      !
         do iclas=1,nclas_chm
            do ipoin=1,npoin
              tmrat_chm(iclas,ipoin) = 0.0_rp       
              if( (coord(1,ipoin).eq.300.).and. &
                  (coord(2,ipoin).eq.300.).and. &
                  (coord(3,ipoin).eq.1200.) ) then     ! Point 5165 global mesh
                  tmrat_chm(iclas,ipoin) = 10.0_rp  
              end if 
            end do 
        end do
      end if
     
      !
      tsour_chm = 1e10
      return
  end if
!
!  Other source options
!
!  1. Sequential or master. Read the source file (in the required format) for the current time slab.
!     Source variables are stored on 2 arrays to minimize communications master-slave
!     INTEGER
!        ndims_src(1)     nsour
!        ndims_src(2)     INT(tsour_chm)
!     REAL
!        srtmp(nsour,1)              x-coordinate
!        srtmp(nsour,2)              y-coordinate
!        srtmp(nsour,3)              z-coordinate
!        srtmp(nsour,3+nclas_chm)    tmrat (in kg/s)
!
  if( INOTSLAVE ) then
     !
     if (kfl_sourc_chm == -2 ) then 
     !
     !  ASCII. File format is
     !        time1 time2
     !        nsour nclass
     !        xsrc ysr zsrc tmrat(nclas)
     !        ....   
     !  
       go_on = .true.
       rewind(lun_resou_chm)
       iline=0
       do while(go_on)
          iline=iline+1
          read(lun_resou_chm,*,err=101) time1,time2
          iline=iline+1
          read(lun_resou_chm,*,err=101) nsour,nclas
          if(nclas/=nclas_chm) call runend(' chm_reasou: wrong number of classes while reading the source file')
          if( oltim>=time1 .and. oltim<=time2 ) then
            !
            !  Time interval found
            !     
            go_on = .false.
            allocate(srtmp(nsour,3+nclas_chm))
            do isour=1,nsour
               iline=iline+1
               read(lun_resou_chm,*,err=101) srtmp(isour,1),srtmp(isour,2),srtmp(isour,3),(srtmp(isour,3+iclas),iclas=1,nclas_chm)
            end do
         else
            !
            !  Time interval not found
            !   
            do isour=1,nsour
               iline=iline+1
               read(lun_resou_chm,*,err=101) dummr
            end do       
         end if
       end do
       ndims_src(1) = nsour
       ndims_src(2) = INT(time2)
       tsour_chm    = time2
!    
    else if (kfl_sourc_chm == -3 ) then 
     !
     !  netCDF    
     !
      call runend('chm_reasrc: netcdf not implemented yet')
     !
    end if
!
!
!
 end if
!  
!
!  2. In parallel flow, master sends information to slaves
!
!
  if( IPARALL ) then
    !
    call parari('BCT',0_ip,SIZE(ndims_src),ndims_src) 
    nsour     = ndims_src(1)
    tsour_chm = 1.0_rp*ndims_src(2)
    !
    if(ISLAVE) then
       allocate(srtmp(ndims_src(1),3+nclas_chm))
    end if
    call pararr('BCT',0_ip,ndims_src(1)*(3+nclas_chm),srtmp) 
    !
  end if
!
!
!  3. Sequential or slave. Interpolation of source points to nodes
!
!
  if( INOTMASTER ) then
     !
     !  Initialization
     !
     do iclas=1,nclas_chm
        do ipoin=1,npoin
           tmrat_chm(iclas,ipoin) = 0.0_rp       
        end do   
     end do
     !
     !  Allocate a vector used to check if 2 slaves share a source point
     !  (source may lay at the frontier)
     !
     allocate(ihave(nsour))
     ihave(1:nsour) = 0
     !
     !  Loop over source points
     !
     do isour = 1,nsour
        
        do idime = 1,ndime
           coglo(idime) = srtmp(isour,idime)
        end do
        !
        !  Loop over elements
        !
        ielem = 0
        go_on = .true.
        do while(go_on)
           ielem = ielem + 1
           pelty = ltype(ielem)
           ptopo = ltopo(pelty) 
           pnode = nnode(pelty)
!
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              do idime = 1,ndime
                 elcod(idime,inode) = coord(idime,ipoin)
              end do  
           end do
!
           call elmgeo_natural_coordinates(      &
                ndime,pelty,pnode,elcod,shapf,   &
                deriv,coglo,coloc,ifoun)
!
           if( ifoun == 1 ) then 
               go_on = .false.
               ihave(isour) = 1
               !
               !  Computes the element volume
               !
               volum = 0.0_rp
               do igaus = 1,ngaus(pelty)
                  call jacobi(&
                       ndime,pnode,elcod,elmar(pelty)%deriv(1,1,igaus),&
                       xjacm,xjaci,gpcar,detjm)
                  volum = volum + elmar(pelty)%weigp(igaus)*detjm
               end do
               !
               !  Extrapolates source to the nodes
               !  Source term is in kg/s so it is necessary to divide by element volume
               !
               do inode = 1,pnode
                  ipoin = lnods(inode,ielem)
                  do iclas = 1,nclas_chm
                     tmrat_chm(iclas,ipoin) = tmrat_chm(iclas,ipoin) + shapf(inode)*srtmp(isour,ndime+iclas)/volum   
                  end do       
               end do
              
           end if
           if( ielem == nelem ) go_on = .false.  
        end do   ! ielem = 1,nelem
      end do     ! isour
      !  
      !  Exchange between slaves for potential sources at the border elements
      !
      call pararr('SLX',NPOIN_TYPE,nclas_chm*npoin,tmrat_chm)
!
  end if   ! end of the interpolation
!
!   4. Check that no source points lay in the border by doing a parallel reduce
!
   if( INOTMASTER ) then
      my_nsour = sum(ihave)
   else
      my_nsour = 0
   end if
   if( IPARALL )  call parari('SUM',0_ip,1_ip,my_nsour)
   if( my_nsour /= nsour ) then
      call runend('chm_reasou: source at the border. Not yet implemented')
   end if   
!
!  5. Everyone release memeory
!
   deallocate(srtmp)
   if( INOTMASTER )  deallocate(ihave)
!
  return

101 continue
    wline = intost(iline)
    call runend('ERROR WHILE READING SOURCE FILE AT LINE '//trim(wline))
!
end subroutine chm_reasou
