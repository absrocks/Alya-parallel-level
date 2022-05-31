!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_par_affinity.f90
!> @author  houzeaux
!> @date    2018-06-23
!> @brief   Affinity
!> @details Displays affinity of Alya
!-----------------------------------------------------------------------

module mod_par_affinity
  use def_master
  use mod_parall
  use mod_lcpsort
  use mod_communications, only : PAR_GATHER
  use mod_outfor, only : outfor
  implicit none
  private

  public :: par_affinity
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-23
  !> @brief   Process affinity
  !> @details This is currently under construction
  !> 
  !-----------------------------------------------------------------------

  subroutine par_affinity()

    character(20)             :: hostname
    character(20), pointer    :: hostname_gat(:)
    integer(ip)               :: nsize,ii,jj,kk,ihost,nhosts,ipart
    type(lptype), allocatable :: ITEM(:),SCRATCH(:)

    return
    
    nsize = npart+1
    nullify(hostname_gat)
    if( IMASTER ) then
       allocate(hostname_gat(0:npart))
    end if
    !
    ! Get hostname and avoid null name
    !
    call GET_ENVIRONMENT_VARIABLE('HOSTNAME',hostname)
    if(trim(hostname)=='') hostname = 'login'
    !
    ! Gather hostnames
    !
    call PAR_GATHER(hostname,hostname_gat)

    if( IMASTER ) then
       call lcpsort_initialization(ITEM,nsize,SCRATCH)
       do ii = 1,nsize
          ITEM(ii) % lcpl =  0
          ITEM(ii) % str  => hostname_gat(ii-1)
       end do
       call lcpsort_sort(ITEM,nsize,SCRATCH)
       call lcpsort_concatenate(ITEM,nsize,nhosts)

       ioutp(1) = nhosts
       call outfor(92_ip,lun_outpu,' ')
       do ihost = 1,nhosts
          kk = 0
          ipart = 1
          do jj = 0,npart
             if( trim(hostname_gat(jj)) == trim(item(ihost) % str) ) then
                if( kk == 0 ) then
                   kk = 1
                   coutp = trim(hostname_gat(jj))
                   call outfor(93_ip,lun_outpu,' ')
                end if
                ipart = ipart + 1
                ioutp(1) = jj
                call outfor(94_ip,lun_outpu,' ')
                if( ipart == 9 ) then
                   write(lun_outpu,*)
                   ipart=1
                end if
             end if
          end do
       end do
       deallocate(hostname_gat)
       call lcpsort_deallocate(ITEM,SCRATCH)
    end if

  end subroutine par_affinity
  
end module mod_par_affinity
  !> @}
  
