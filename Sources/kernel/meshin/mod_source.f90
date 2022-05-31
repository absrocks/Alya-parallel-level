module mod_source

contains

  subroutine readsour(ndim,nsour,rsuni,rscal,rsour,rsgeo,isizcrit,ismoo,&
             iinter,icartout,irefsurf,rtolsca,isizuni)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip), intent(in)   :: ndim
    integer(ip),intent(inout) :: nsour,isizcrit,ismoo,iinter,icartout
    integer(ip),intent(inout) :: irefsurf,isizuni
    real(rp),intent(inout)    :: rsuni,rscal,rtolsca 
    real(rp),pointer          :: rsour(:,:),rsgeo(:,:,:) 
    integer(ip)               :: isour,iostatus
    character*80              :: ntext
    integer(4)                :: istat
    !
    !     This subroutine reads the source parameters 
    !
    !
    !     rsuni: the size of the first level of cells
    !     rscal: scaling factor for the size
    !     isizcrit: size criterion for cartesian mesh:   
    !             - 1 keep the size as is
    !             - 2 apply minimum between cell size and uniform size 
    !             - 3 apply uniform size
    !             - 4 apply source size
    !     ismoo: do we want to smooth the size? 
    !     iinter: do we want to check intersection between faces and
    !             cells
    !     icartout: do we want to output the cartesian mesh
    !     isizuni: do we want to apply a uniform size to the cartesian mesh 
    !     irefsurf: do we want to refine the cartesian mesh according 
    !     to the surface mesh
    !     rtolsca: tolerance for scalar product  
    ! 
    !
    open(unit=77,file='source.dat',status='old',iostat=iostatus)
    if(iostatus>0)then
       write(*,*)'Source file not found in readsour'
       stop
    endif
    rewind 77

    write(*,*)'Reading nsour.dat'
    !
    !     Read points
    !
    read(77,*)ntext,ntext,rsuni
    read(77,*)ntext,ntext,rscal
    read(77,*)ntext,ntext,isizcrit
    read(77,*)ntext,ismoo
    read(77,*)ntext,iinter
    read(77,*)ntext,icartout
    read(77,*)ntext,isizuni
    read(77,*)ntext,irefsurf
    read(77,*)ntext,rtolsca
    read(77,*)ntext,nsour

    write(*,*)'rsuni=',rsuni 
    write(*,*)'rscal=',rscal 
    write(*,*)'isizcrit=',isizcrit 
    write(*,*)'ismoo=',ismoo 
    write(*,*)'iinter=',iinter 
    write(*,*)'icartout=',icartout 
    write(*,*)'irefsurf=',irefsurf 
    write(*,*)'isizuni=',isizuni 
    write(*,*)'rtolsca=',rtolsca 
    write(*,*)'nsour=',nsour 

    if(nsour<0)then
       write(*,*)'Error in readsour, nsour negative'
       stop
    endif

    if(nsour>0)then

       allocate(rsgeo(ndim,ndim,nsour),stat=istat)
       call memchk(zero,istat,memor_msh,'RSGEO','readsour',rsgeo)

       allocate(rsour(3,nsour),stat=istat)
       call memchk(zero,istat,memor_msh,'RSOUR','readsour',rsour)

       do isour=1,nsour

          !
          !     Read coordinates of the source triangles
          !
          read(77,*)rsgeo(1,1,isour),rsgeo(2,1,isour),rsgeo(3,1,isour)
          read(77,*)rsgeo(1,2,isour),rsgeo(2,2,isour),rsgeo(3,2,isour)
          read(77,*)rsgeo(1,3,isour),rsgeo(2,3,isour),rsgeo(3,3,isour)
          !
          !     Read parameters of the source:
          !     - minimum element size
          !     - radius of constant size
          !     - linear parameter
          !
          read(77,*)rsour(1,isour),rsour(2,isour),rsour(3,isour)

       enddo

    else

       allocate(rsgeo(1,1,1),stat=istat)
       call memchk(zero,istat,memor_msh,'RSGEO','readsour',rsgeo)
       allocate(rsour(1,1),stat=istat)
       call memchk(zero,istat,memor_msh,'RSOUR','readsour',rsour)

    endif

    write(*,*)'End reading nsour.dat'

    close(77)

  end subroutine readsour

  subroutine chksour(ndim,nsour,rsgeo,rsour)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)     :: nsour,ndim
    real(rp),intent(in)        :: rsour(3,nsour)
    real(rp), intent(inout)    :: rsgeo(ndim,ndim,nsour)
    integer(ip)                  :: isour,lcont(3),ncont
    real(rp)   :: rpax,rpay,rpaz,rpbx,rpby,rpbz,rpcx,rpcy,rpcz,rdabx,rdaby,rdabz
    real(rp)   :: rdacx,rdacy,rdacz,rdcbx,rdcby,rdcbz
    real(rp)   :: rdistab,rdistac,rdistcb,epsil,rtol
    !
    !     This sub modifies the sources for point or line degeneracy
    !

    rtol=1.0d-08

    do isour=1,nsour
       rpax=rsgeo(1,1,isour)
       rpay=rsgeo(2,1,isour)
       rpaz=rsgeo(3,1,isour)
       rpbx=rsgeo(1,2,isour)
       rpby=rsgeo(2,2,isour)
       rpbz=rsgeo(3,2,isour)
       rpcx=rsgeo(1,3,isour)
       rpcy=rsgeo(2,3,isour)
       rpcz=rsgeo(3,3,isour)

       rdabx=rpax-rpbx
       rdaby=rpay-rpby
       rdabz=rpaz-rpbz
       rdistab=sqrt(rdabx*rdabx+rdaby*rdaby+rdabz*rdabz)

       rdacx=rpax-rpcx
       rdacy=rpay-rpcy
       rdacz=rpaz-rpcz
       rdistac=sqrt(rdacx*rdacx+rdacy*rdacy+rdacz*rdacz)

       rdcbx=rpcx-rpbx
       rdcby=rpcy-rpby
       rdcbz=rpcz-rpbz
       rdistcb=sqrt(rdcbx*rdcbx+rdcby*rdcby+rdcbz*rdcbz)
       !
       !     Compute the local epsilon
       !
       epsil=rsour(1,isour)*rtol
       !
       !     Compares to local distance
       !
       ncont=0_ip

       if(rdistab<epsil)then
          ncont=ncont+1
          lcont(ncont)=1_ip
       endif

       if(rdistac<epsil)then
          ncont=ncont+1
          lcont(ncont)=2_ip
       endif

       if(rdistcb<epsil)then
          ncont=ncont+1
          lcont(ncont)=3_ip
       endif


       !
       !     Do we have some degeneracy?
       !

       if(ncont==0)cycle

       !
       !     Do we have a point degeneracy?
       !
       if(ncont==3)then

          !
          !     Consider the first point and moves the others from epsil
          !





          !
          !     Do we have a line degeneracy?
          !
       else if(ncont==1)then 



       endif








    enddo


  end subroutine chksour


end module mod_source


