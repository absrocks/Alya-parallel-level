module mod_ker_tendencies
  use def_kintyp, only : rp, ip

  implicit none
  save
  logical    , public               :: kfl_tendencies_ker

  ! pointers
  real(rp),    public, pointer      :: ten_ugeos(:,:,:)  ! 2,nz,nsteps
  real(rp),    public, pointer      :: ten_uadve(:,:,:)  ! 2,nz,nsteps
  real(rp),    public, pointer      :: ten_thadv(:,:)    ! nz, nsteps
  real(rp),    public, pointer      :: ten_theta(:,:)    ! nz, nsteps
  real(rp),    public, pointer      :: ten_heigh(:,:)    ! nz, nsteps
  real(rp),    public, pointer      :: ten_time(:)       ! time in secs, nsteps
  integer(ip), public, pointer      :: ten_nzlev(:)      ! number of zlevels nsteps

  ! module variables
  integer(ip), public               :: ten_nstep, &      ! number of steps
       ten_nzlev_max      ! number of z levels

  integer(ip), public               :: up_istep_last = 1 ! initialization

  ! functions
  public  :: read_tendencies
  public  :: get_tendencies_u
  private :: get_tendencies_u_istep

contains

  subroutine  read_tendencies()
    ! reads tendencies and allocates structures
    use def_inpout,              only : words, param, exists
    use mod_ecoute,              only : ecoute
    use mod_memory
    use def_master,              only : mem_modul, modul
    implicit none

    integer(ip)    :: istep, izlev

    ten_nzlev_max = 500 ! initialization
    
    do while( words(1) /= 'STEPS' )
       call ecoute('ker_readat')  ! STEPS
    end do
    
    ten_nstep = param(1)         ! number of steps
    if (words(2) == 'NZLEV')  ten_nzlev_max = param(2)
    ! 
    ! allocate structures
    !
    call memory_alloca(mem_modul(1:2,modul),'TEN_UGEOS','read_tendencies',ten_ugeos, 2_ip, ten_nzlev_max,ten_nstep)
    call memory_alloca(mem_modul(1:2,modul),'TEN_UADVE','read_tendencies',ten_uadve, 2_ip, ten_nzlev_max,ten_nstep)
    call memory_alloca(mem_modul(1:2,modul),'TEN_THADV','read_tendencies',ten_thadv, ten_nzlev_max,ten_nstep)
    call memory_alloca(mem_modul(1:2,modul),'TEN_THETA','read_tendencies',ten_theta, ten_nzlev_max,ten_nstep)
    call memory_alloca(mem_modul(1:2,modul),'TEN_HEIGH','read_tendencies',ten_heigh, ten_nzlev_max,ten_nstep)
    call memory_alloca(mem_modul(1:2,modul),'TEN_TIME' ,'read_tendencies',ten_time , ten_nstep)
    call memory_alloca(mem_modul(1:2,modul),'TEN_NZLEV','read_tendencies',ten_nzlev, ten_nstep)
! read from:
!    /gpfs/scratch/bsc21/bsc21811/Alya-Runs/Hornamossen> more Hornamossen_tendencies/Hornamossen_tendencies
    istep = 0
    call ecoute('ker_readat')   ! TIME
    do while(words(1) /= 'ENDST' )
       do while(words(2) /= 'TIME' )
          call ecoute('ker_readat')   ! TIME
       end do
       istep = istep +1 
       ten_time(istep) = param(2)  ! time  ! NEW FILE
       call ecoute('ker_readat')   ! blank  
       izlev = 0
       do while(words(2) /= 'ENDTI' )
          izlev = izlev +1 
          call ecoute('ker_readat')
          ten_heigh(izlev, istep) = param(1) ! should be given in increasing order
          ten_ugeos(1:2,izlev, istep) = param(2:3)
          ten_uadve(1:2,izlev, istep) = param(4:5)
          ten_thadv(izlev,istep)      = param(6)
          ten_theta(izlev,istep)      = param(7)          
       end do
       ten_nzlev(istep) = izlev
!       do izlev =1, ten_nzlev(istep)
          
!       end do
       call ecoute('ker_readat')
    end do      

    
  end subroutine read_tendencies
  subroutine  get_tendencies_u_istep(zcoor,istep,gpten_geo, gpten_uadv)
    ! interpolate tendency value to zcoor height, for step istep
    !
    !  for interpolation, height values are supposed to be given in increasing order
    !
    implicit none
    real(rp), intent(in)    :: zcoor   ! z coord
    real(rp), intent(out)   :: gpten_geo(2), gpten_uadv(2)
    integer(ip), intent(in) :: istep
    ! local variables
    integer(ip)            ::  iz, jz, kz
    real(rp)               ::  facto

    

    ! for a given time step
    if (zcoor.lt.ten_heigh(1, istep)) then ! coorditate less than minimum 
       gpten_geo =0.0_rp
       gpten_uadv = 0.0_rp
!    else if(zcoor.gt.z_top )  then 
    else
       iz =1
       jz =ten_nzlev(istep)
       kz =jz/2
       do while ((jz-iz).gt.1)  ! bisection method
          if (zcoor.lt.ten_heigh(kz,istep)) then
             jz = kz                  
          else
             iz = kz
          end if
          kz = (iz+jz)/2
       end do
       !interpolation factor
       facto = (zcoor-ten_heigh(iz, istep))/(ten_heigh(jz, istep) -  ten_heigh(iz, istep))
       gpten_geo(1:2) = ten_ugeos(1:2,iz, istep) + facto*( ten_ugeos(1:2,jz, istep) - ten_ugeos(1:2,iz, istep) )
    end if
    
  end subroutine get_tendencies_u_istep

  subroutine  get_tendencies_u(zcoor,ctime,gpten_geo, gpten_uadv)
    ! interpolate tendency value to zcoor height and ctime
    !
    !
    implicit none
    real(rp), intent(in)    :: zcoor, ctime   ! z coord
    real(rp), intent(out)   :: gpten_geo(2), gpten_uadv(2)

    real(rp)                :: geo_old(2), uadv_old(2), geo_new(2), uadv_new(2)
    real(rp)                :: facto
    !
    ! time interpolation
    !
    do while (ctime.gt.ten_time(up_istep_last))  ! check if new time limits
         up_istep_last= up_istep_last +1
    end do
   
    facto = (ctime - ten_time(up_istep_last-1))/  (ten_time(up_istep_last) - ten_time(up_istep_last-1))

    call get_tendencies_u_istep(zcoor, up_istep_last -1, geo_old, uadv_old)
    call get_tendencies_u_istep(zcoor, up_istep_last   , geo_new, uadv_new)

    gpten_geo = facto*geo_new  + (1.0_rp-facto)*geo_old
    gpten_uadv = facto*uadv_new  + (1.0_rp-facto)*uadv_old
    
  end subroutine get_tendencies_u
  
  
end module mod_ker_tendencies
