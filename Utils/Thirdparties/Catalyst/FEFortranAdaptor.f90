module tcp
  use iso_c_binding
  implicit none
  public
  interface tcp_adaptor
    module procedure testcoprocessor
    end interface
    contains

      subroutine testcoprocessor(itste,ttime,npoin,coord,nelem,lnods,ltype,rank,veloc,press)
      !subroutine testcoprocessor(itste,ttime,npoin,coord,nelem,lnods,ltype,rank)  
        use iso_c_binding
        use def_kintyp
        !use def_domain

        implicit none
        integer(ip), intent(in)                    :: itste,nelem,npoin,rank
        real(rp),    intent(in)                    :: ttime
        real(rp),    dimension(:,:), intent (in)   :: coord
        integer(ip), dimension(:,:), intent (in)        :: lnods
        integer(ip), dimension(:),   intent (in)        :: ltype
        integer(ip)                                     :: flag
        !---------------------------------------------------------
        real(rp),    dimension(:,:,:), intent (in)   :: veloc
        real(rp),    dimension(:,:), intent (in)     :: press
        !--------------------------------------------------------- 
        integer(ip)             :: inode,ielem,idime,ipoin
        integer(ip)             :: offset,istat,icount,pnode
        !-----------------------------------------------
        !ojo 
        !lnods_tmp and offset_tmp had to be integer double precision 8 bytes
        !ltype_tmp had to be integer simple precison 2 bytes
        !  
        !----------------------------------------------- 
        integer(8), pointer       :: lnods_tmp(:),offset_tmp(:)
        integer(4), pointer       :: ltype_tmp(:)
        !---------------------------------------------------------
        real(rp), pointer         :: pressure(:),velocity(:)
        !---------------------------------------------------------

        !if (rank==0)then
        !   return
        !endif
        
        nullify(lnods_tmp)
        nullify(ltype_tmp)
        nullify(offset_tmp)
        !
        nullify(pressure)
        nullify(velocity)
        
        
        offset=1
        icount=0
        do ielem=1,nelem
           if       ( ltype(ielem) == 30 ) then ! TETRA
              icount=icount+4+1
           else if  ( ltype(ielem) == 37 ) then ! HEXA
              icount=icount+8+1
           else if  ( ltype(ielem) == 32 ) then ! PYRA
              icount=icount+5+1
           else if  ( ltype(ielem) == 34 ) then ! PENTA   
              icount=icount+6+1
           endif
        enddo
        !
        ! allocate memory temporary
        !
        allocate (lnods_tmp(icount),stat=istat)
        allocate (offset_tmp(nelem),stat=istat)
        allocate (ltype_tmp(nelem),stat=istat)
        !
        allocate (pressure(npoin),stat=istat)
        allocate (velocity(3*npoin),stat=istat)
        !
        ! CONNEC
        !
        icount=0
        do ielem=1,nelem
           if       ( ltype(ielem) == 30 ) then ! TETRA
              pnode = 4
              ltype_tmp(ielem) = 10
           else if  ( ltype(ielem) == 37 ) then ! HEXA
              pnode = 8
              ltype_tmp(ielem) = 12
           else if  ( ltype(ielem) == 32 ) then ! PYRA
              pnode = 5
              ltype_tmp(ielem) = 14
           else if  ( ltype(ielem) == 34 ) then ! PENTA
              pnode = 6
              ltype_tmp(ielem) = 13
           end if
           icount=icount+1
           lnods_tmp(icount)=pnode
           do inode=1,pnode
              icount=icount+1
              lnods_tmp(icount)=lnods(inode,ielem)-1  
           end do
           offset_tmp(ielem) = offset
           offset = offset + pnode +1
        enddo
        !
        ! attribut test
        !
        do ipoin=1,npoin
           pressure(ipoin)=press(ipoin,1)
        enddo
        do ipoin=1,npoin
           velocity(ipoin)=veloc(1,ipoin,1)
           velocity(ipoin+npoin)=veloc(2,ipoin,1)
           velocity(ipoin+npoin*2)=veloc(3,ipoin,1)
        enddo
        !
        !   cheating way
        !
        !do ipoin=1,npoin
        !   pressure(ipoin)=real(rank)
        !enddo
        !do ipoin=1,npoin
        !   velocity(ipoin)=real(rank)
        !   velocity(ipoin+npoin)=real(rank+1)
        !   velocity(ipoin+npoin*2)=real(rank+2)
        !enddo
        !-----check----------------------------------------------
        !write(*,*)'nelem=',nelem
        !write(*,*)'npoin=',npoin
        !write(*,*)'icount=',icount
        !write(*,*)'size(lype_tmp)=',size(ltype_tmp)
        !write(*,*)'lype_tmp=',ltype_tmp
        !write(*,*)'size(lnods_tmp)=',size(lnods_tmp)
        !write(*,*)'lnods_tmp=',lnods_tmp
        !write(*,*)'size(offset_tmp)=',size(offset_tmp)
        !write(*,*)'offset_tmp=',offset_tmp
        !----------------------------------------
        !-----writing ouptup / partition--------------------------
        !write(rank+700,*)'nelem=',nelem
        !write(rank+700,*)'npoin=',npoin
        !write(rank+700,*)'icount=',icount
        !write(rank+700,*)'lype_tmp=',ltype_tmp
        !write(rank+700,*)'lnods_tmp=',lnods_tmp
        !write(rank+700,*)'coord=',coord
        !if (itste==10)write(rank+700,*)'velocity=',veloc(1:3,1:1,1)
        !if (itste==10)write(rank+700,*)'pressure=',press(1:1,1)
        !----------------------------------------
        call requestdatadescription(itste,ttime,flag)
        if (flag .ne. 0) then
           call catalystcoprocess(npoin, coord, nelem, ltype_tmp, icount, lnods_tmp, offset_tmp, ttime, itste, pressure, velocity)
           call coprocess()
        end if                         
        
        !--------------------------------------------------------
        deallocate(lnods_tmp)
        deallocate(offset_tmp)
        deallocate(ltype_tmp)
        !
        deallocate(pressure)
        deallocate(velocity)

        return

      end subroutine testcoprocessor
end module tcp

