module def_postpr
  !-----------------------------------------------------------------------
  !    
  ! Heading for the postprocess subroutines
  !
  !-----------------------------------------------------------------------
  use def_master
  use def_kermod
  use def_domain, only : ndime
  implicit none
  !
  ! General
  ! 
  integer(ip)              :: &
       npart_pos,             &      ! # of geometry part
       kfl_ivari(5)
  !
  ! Element name for post
  !
  character(15)            :: &
       cepos(60)                     ! List of element names (upper case)
  character(5)             :: &
       varna_pos(2,100)              ! Variables names to postprocess
  character(5)             :: &
       wopow(2)                      ! Variables names to postprocess
  integer(ip)              :: &      
       varnu_pos                     ! Number of variables to postprocess
  character(20)            :: &
       nunam_pos                     ! Postprocess name
  !
  ! Types
  ! 
  type partt_pos 
     character(50) :: name
     integer(ip)   :: numepart
     integer(ip)   :: npoin
  end type partt_pos
  !
  ! Parts data base
  !
  type(partt_pos)  :: parts_pos(10)
  !
  ! ENSIGHT special variables
  !
  integer(ip)              :: &
       nppva_ens,             &      ! Ensight postprocess time counter
       kfl_statu_ens                 ! Ensight postprocess status           
  character(30)            :: &      ! Ensight variable type
       ensty_ens(50)
  character(15)            :: &      ! Ensight variable name
       ensva_ens(50)
  character(150)           :: &      ! Ensight variable file name
       ensfi_ens(2,50)
  !
  ! VU special variables
  !
  integer(ip)              :: &
       ncoun_pos                     ! VU postprocess position counter
  !
  ! VTK special variables
  !
  integer(ip)              :: &
       vtk_time                     ! VTK postprocess flag time
  !
  ! Variables for voxel output (Blender format)
  !
  integer(ip)              :: &
       totvx,                 &      ! Total number of voxes detected inside mesh
       nvoxl                         ! Number of voxels inside the domain
  integer(ip),pointer      :: &    
       nslvx(:),              &      ! Number of voxels that each slave has
       dplvx(:)                      ! Displacement list to locate data coming from slaves

  integer(ip),pointer      :: &
       idxvx(:),              &     ! Master index of voxels
       elevx(:)                     ! Element that contains each voxel
       
  real(rp),pointer         :: &
       vxdat(:)                     ! Voxel data

  interface voxpos
     module procedure tr2pos,id2pos
  end interface

contains

  !
  ! Functions for voxel indexing and positioning
  !

  function idx2tr(indexi) result(triple)
    integer(ip),intent(in) :: indexi
    integer(ip)            :: triple(3)

    triple(3) = 1 + int((indexi-1)/(resvx(1)*resvx(2)))
    triple(2) = 1 + int((indexi-1-(triple(3)-1)*resvx(1)*resvx(2))/resvx(1))
    triple(1) = indexi-( (triple(2)-1)*resvx(1)+(triple(3)-1)*resvx(1)*resvx(2) )

  end function idx2tr

  function tr2idx(triad) result(indexi)
    integer(ip)            :: indexi
    integer(ip),intent(in) :: triad(3)

    indexi = triad(1) + (triad(2)-1)*resvx(1) + (triad(3)-1)*resvx(1)*resvx(2)

  end function tr2idx

  function tr2pos(triple) result(vxpos)
    real(rp)               :: vxpos(3)
    integer(ip),intent(in) :: triple(3)

    vxpos(1) = bobvx(1,1) + real(triple(1)-1_ip,rp) * ( bobvx(2,1)-bobvx(1,1) ) / real(resvx(1)-1,rp)
    vxpos(2) = bobvx(1,2) + real(triple(2)-1_ip,rp) * ( bobvx(2,2)-bobvx(1,2) ) / real(resvx(2)-1,rp)
    vxpos(3) = bobvx(1,3) + real(triple(3)-1_ip,rp) * ( bobvx(2,3)-bobvx(1,3) ) / real(resvx(3)-1,rp)

  end function tr2pos

  function id2pos(indexi) result(vxpos)
    integer(ip),intent(in) :: indexi
    real(rp)               :: vxpos(3)
    integer(ip)            :: triple(3)
  
    triple(3) = 1+int((indexi-1)/(resvx(1)*resvx(2)))
    triple(2) = 1+int ( (indexi-1-(triple(3)-1)*resvx(1)*resvx(2))/resvx(1))
    triple(1) = indexi-( (triple(2)-1)*resvx(1)+(triple(3)-1)*resvx(1)*resvx(2) )

    vxpos(1)  = bobvx(1,1) + (triple(1)-1) * ( bobvx(2,1)-bobvx(1,1) ) / (max(resvx(1),2_ip)-1) + travx(1)
    vxpos(2)  = bobvx(1,2) + (triple(2)-1) * ( bobvx(2,2)-bobvx(1,2) ) / (max(resvx(2),2_ip)-1) + travx(2)
    vxpos(3)  = bobvx(1,3) + (triple(3)-1) * ( bobvx(2,3)-bobvx(1,3) ) / (max(resvx(3),2_ip)-1) + travx(3)

  end function id2pos

  function nvoxsl(slave) result(numvx) !Computes the number of voxels in a given slave
    integer(ip),intent(in) :: slave
    integer(ip)            :: numvx
    numvx = nslvx(slave+1)-nslvx(slave)
  end function nvoxsl


end module def_postpr
