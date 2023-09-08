!----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_color_communicators.f90
!> @author  Guillaume Houzeaux
!> @date    01/02/2014
!> @brief   Max number of zones and subdomains
!> @details Compute the max number of zones and subdomains in order
!>          to define the mapping (code,zone,subd) => color
!>
!>          \verbatim
!>
!>          I_AM_IN_COLOR(ICOLO) .................................... if I have color ICOLO (TRUE/FALSE)
!>          PAR_COMM_COLOR(0:MCOLO,0:MCOLO) ......................... Intercolor communicator
!>          PAR_COMM_COLOR_ARRAY(0:MCOLO) ........................... Color communication arrays
!>          PAR_CPU_TO_COLOR(0:PAR_WORLD_SIZE) % L(:) ............... List of colors for each world partition
!>          PAR_COLOR_TO_CPU(0:MCOLO) % L(:) ........................ List of world partitions for each color
!>          PAR_COMM_COLOR_PERM(0:MCOLO,0:MCOLO,0:PAR_WORLD_SIZE) ... Ranks for each communicator
!>          PAR_COMM_WORLD_TO_CODE_PERM(2,0:PAR_WORLD_SIZE) ......... Rank permutation from world to code
!>
!>          \endverbatim 
!> @} 
!----------------------------------------------------------------------

subroutine par_color_communicators()
  use def_kintyp,         only    :  ip,rp,lg
  use mod_parall,         only    :  PAR_WORLD_SIZE
  use mod_parall,         only    :  PAR_MY_CODE_RANK, PAR_COMM_CURRENT
  use mod_parall,         only    :  PAR_COMM_WORLD_TO_CODE_PERM
  use mod_parall,         only    :  I_AM_IN_COLOR
  use mod_parall,         only    :  PAR_COMM_COLOR
  use mod_parall,         only    :  PAR_CPU_TO_COLOR
  use mod_parall,         only    :  PAR_COLOR_TO_CPU
  use mod_parall,         only    :  PAR_COMM_COLOR_PERM
  use mod_parall,         only    :  PAR_COMM_COLOR_ARRAY
  use mod_parall,         only    :  PAR_COMM_MY_CODE
  use mod_parall,         only    :  PAR_COMM_WORLD
  use mod_parall,         only    :  PAR_MY_WORLD_RANK
  use mod_parall,         only    :  PAR_COMM_MY_CODE_ARRAY
  use mod_parall,         only    :  PAR_INITIALIZE_COMMUNICATION_ARRAY
  use mod_parall,         only    :  PAR_COPY_COMMUNICATION_ARRAY
  use mod_parall,         only    :  mcode
  use mod_parall,         only    :  msubd
  use mod_parall,         only    :  mzone
  use mod_parall,         only    :  mcolo
  use mod_parall,         only    :  ncolo
  use mod_parall,         only    :  par_code_zone_subd_to_color
  use mod_parall,         only    :  par_color_to_subd
  use mod_parall,         only    :  par_color_to_zone
  use mod_parall,         only    :  par_color_to_code
  use mod_parall,         only    :  par_memor
  use mod_communications, only    :  PAR_MAX
  use mod_communications, only    :  PAR_SUM
  use mod_communications, only    :  PAR_ALLGATHER
  use mod_communications, only    :  PAR_ALLGATHERV
  use mod_communications, only    :  PAR_COMM_SPLIT
  use mod_communications, only    :  PAR_SUM_ALL
  use mod_communications, only    :  PAR_MAX_ALL
  use def_master,         only    :  current_code, kfl_paral
  use def_master,         only    :  ioutp,ISEQUEN,lun_outpu
  use def_master,         only    :  I_AM_IN_ZONE
  use def_master,         only    :  I_AM_IN_SUBD
  use def_domain,         only    :  nzone,nsubd
  use mod_memory,         only    :  memory_alloca
  use mod_outfor,         only    :  outfor
  use mod_messages,       only    :  livinf
  use mod_iofile,         only    :  iofile_flush_unit

  use def_master
  use mod_communications

  implicit none
  integer(ip)                     :: icode,izone,isubd,icolo,jcolo,ierro
  integer(ip)                     :: my_new_rank,ipart,kcolo,icolor
  integer(ip)                     :: jsubd,jzone,m1,cz,ii,jcode,jpart
  logical(lg)                     :: isplit

  integer(ip),         pointer    :: COMM_MATRIX(:)    
  integer(ip),         pointer    :: COMM_MATRIX_TMP(:)   

  integer(ip),         pointer    :: number_colors(:)
  integer(ip),         pointer    :: number_parts(:)
  integer(ip),         pointer    :: my_list_of_colors(:)
  integer(ip),         pointer    :: list_of_colors(:)
  integer(ip),         pointer    :: in_subd(:)
  integer(ip),         pointer    :: in_zone(:)

  integer(4)                      :: recvcount4
  integer(ip),         pointer    :: sendbuf(:)
  integer(ip),         pointer    :: recvbuf(:)

  real(rp) :: time1,time2

  !if( ISEQUEN ) return

  call livinf(0_ip,'PARALL: COMPUTE INTER-COLOR COMMUNICATORS',0_ip)
  ! 
  ! Nullify local pointers 
  !
  nullify(COMM_MATRIX)
  nullify(COMM_MATRIX_TMP)
  nullify(number_colors)
  nullify(number_parts)
  nullify(my_list_of_colors)
  nullify(list_of_colors) 
  nullify(in_subd) 
  nullify(in_zone)
  nullify(sendbuf)
  nullify(recvbuf)
  !
  ! Compute MZONE and MSUBD. MCODE was computed at beginning of 
  ! code, after splitting the MPI_COMM_WORLD
  !
  mzone = nzone
  msubd = nsubd

  call PAR_MAX(mzone,'IN THE WORLD')
  call PAR_MAX(msubd,'IN THE WORLD')
  mcolo = (mcode+1)*(mzone+1)*(msubd+1)-1
  m1    = mcolo+1_ip
  cz    = int(PAR_WORLD_SIZE,ip)
  ! 
  ! Allocate local memory 
  !
  allocate( COMM_MATRIX(0:mcolo) )
  COMM_MATRIX = 0
  !
  ! Allocate global memory 
  ! 
  call memory_alloca(par_memor,'I_AM_IN_COLOR',              'par_color_communicatros',I_AM_IN_COLOR,      m1,             'INITIALIZE',0_ip)            
  call memory_alloca(par_memor,'PAR_COMM_COLOR',             'par_color_communicatros',PAR_COMM_COLOR,     m1,m1,          'INITIALIZE',0_ip,0_ip)       
  call memory_alloca(par_memor,'PAR_CPU_TO_COLOR',           'par_color_communicatros',PAR_CPU_TO_COLOR,   cz,             'INITIALIZE',0_ip)            
  call memory_alloca(par_memor,'PAR_COLOR_TO_CPU',           'par_color_communicatros',PAR_COLOR_TO_CPU,   m1,             'INITIALIZE',0_ip)            
  call memory_alloca(par_memor,'PAR_COMM_COLOR_PERM',        'par_color_communicatros',PAR_COMM_COLOR_PERM,m1,m1,cz,       'INITIALIZE',0_ip,0_ip,0_ip)  
  call memory_alloca(par_memor,'PAR_COMM_WORLD_TO_CODE_PERM','par_color_communicatros',PAR_COMM_WORLD_TO_CODE_PERM,2_ip,cz,'INITIALIZE',1_ip,0_ip)       

  if( ISEQUEN ) then
     allocate( PAR_COMM_MY_CODE_ARRAY(1) )
     call PAR_INITIALIZE_COMMUNICATION_ARRAY(PAR_COMM_MY_CODE_ARRAY(1))
  end if

  allocate( PAR_COMM_COLOR_ARRAY(0:mcolo) )
  call PAR_INITIALIZE_COMMUNICATION_ARRAY(PAR_COMM_COLOR_ARRAY)
  !
  ! Check CURRENT_CODE does not exceed mcode
  !
  ierro = 0
  if( current_code > mcode ) ierro = 1
  !   call PAR_MAX(ierro,'IN THE UNIVERSE')
  call PAR_MAX_ALL(ierro)
  if( ierro == 1 ) call runend('PAR_COLOR_COMMUNICATORS: CODE NUMBER EXCEED NUMBER OF CODES')

  !----------------------------------------------------------------------
  !   
  !   Using zone and subdomain information, fill in the color arrays:
  !
  !   I_AM_IN_COLOR(icolo) = .true./.false.
  !   COMM_MATRIX(icolo)   = 1/0
  !
  !   I have zone 3,5 and subd 2; 0 is all
  !   A priori we will communicate only between zones or between subdomains
  !   Structure is prepared to do combinations
  !
  !   +---+---+---+---+---+---+
  ! 0 |   |   |   | x |   | x | 
  !   +---+---+---+---+---+---+
  ! 1 |   |   |   |   |   |   |
  !   +---+---+---+---+---+---+
  ! 2 | x |   |   |   |   |   |
  !   +---+---+---+---+---+---+
  ! 3 |   |   |   |   |   |   |
  !   +---+---+---+---+---+---+
  !     0   1   2   3   4   5   => zone
  ! 
  !----------------------------------------------------------------------

  if( PAR_MY_CODE_RANK /= 0 ) then
     !
     ! Zones and subdomains inside my code
     !
     isubd = 0
     do izone = 1,nzone
        if( I_AM_IN_ZONE(izone) ) then
           icode                = current_code
           icolo                = par_code_zone_subd_to_color(icode,izone,isubd)
           COMM_MATRIX(icolo)   = 1
        end if
     end do
     izone = 0
     do isubd = 1,nsubd
        if( I_AM_IN_SUBD(isubd) ) then
           icode                = current_code
           icolo                = par_code_zone_subd_to_color(icode,izone,isubd)
           COMM_MATRIX(icolo)   = 1
        end if
     end do
     !
     ! I am in code, of course!
     !
     icolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip)
     COMM_MATRIX(icolo) = 1
     !
     ! We are obviously in COMM_WORLD tool
     !
     COMM_MATRIX(0) = 1

  end if
  !
  ! Involve the master in all communications
  !
  allocate( COMM_MATRIX_TMP(0:mcolo) )
  COMM_MATRIX_TMP = COMM_MATRIX
  call PAR_MAX(COMM_MATRIX_TMP,'IN MY CODE')
  if( PAR_MY_CODE_RANK == 0 ) COMM_MATRIX = COMM_MATRIX_TMP
  deallocate( COMM_MATRIX_TMP )

  do icolo = 0,mcolo
     if( COMM_MATRIX(icolo) == 1 ) I_AM_IN_COLOR(icolo) = .true.
  end do

  !----------------------------------------------------------------------
  !
  ! PAR_CPU_TO_COLOR(IPART) % L(:) = List of colors of partition IPART
  ! PAR_COLOR_TO_CPU(ICOLO) % L(:) = List of partitions in color ICOLO
  !
  !----------------------------------------------------------------------
  !
  ! MY_LIST_OF_COLORS(1:NCOLO) = My list of colors
  !
  ncolo = 0
  do icolo = 0,mcolo
     if( COMM_MATRIX(icolo) == 1 ) ncolo = ncolo + 1
  end do
  allocate( my_list_of_colors(ncolo) )
  ncolo = 0
  do icolo = 0,mcolo
     if( COMM_MATRIX(icolo) == 1 ) then
        ncolo = ncolo + 1
        my_list_of_colors(ncolo) = icolo
     end if
  end do
  !
  ! NUMBER_COLORS(1:PAR_WORLD_SIZE) = Number of colors of all partitions
  !
  allocate( number_colors(0:PAR_WORLD_SIZE-1) )
  call PAR_ALLGATHER(ncolo,number_colors,1_4,'IN THE WORLD')
  !
  ! LIST_OF_COLORS(1:PAR_WORLD_SIZE) = All gather list of colors of all partitions
  ! 
  kcolo = 0
  do ipart = 0,PAR_WORLD_SIZE-1
     kcolo = kcolo + number_colors(ipart)
  end do
  allocate( list_of_colors(kcolo) )   
  call PAR_ALLGATHERV(my_list_of_colors,list_of_colors,number_colors,'IN THE WORLD')
  !
  ! PAR_CPU_TO_COLOR(IPART) % L(:) = List of colors
  !
  kcolo = 0
  do ipart = 0,PAR_WORLD_SIZE-1
     if( number_colors(ipart) > 0 ) then
        call memory_alloca(par_memor,'PAR_CPU_TO_COLOR(IPART) % L','par_color_communicatros',PAR_CPU_TO_COLOR(ipart) % l,number_colors(ipart))
        do icolo = 1,number_colors(ipart)
           kcolo = kcolo + 1
           PAR_CPU_TO_COLOR(ipart) % l(icolo) = list_of_colors(kcolo)
        end do
        call heapsorti1(2_ip,number_colors(ipart),PAR_CPU_TO_COLOR(ipart) % l)
     else
        nullify(  PAR_CPU_TO_COLOR(ipart) % l )        
     end if
  end do
  if( associated(number_colors)     ) deallocate( number_colors )
  if( associated(list_of_colors)    ) deallocate( list_of_colors )
  if( associated(my_list_of_colors) ) deallocate( my_list_of_colors )
  !
  ! PAR_COLOR_TO_CPU(ICOLO) % L(:) = List of partitions
  !
  allocate( number_parts(0:mcolo) )
  do icolo = 0,mcolo
     number_parts(icolo) = 0
  end do
  do ipart = 0,PAR_WORLD_SIZE-1
     if( associated(PAR_CPU_TO_COLOR(ipart) % l) ) then 
        do kcolo = 1,size(PAR_CPU_TO_COLOR(ipart) % l,kind=ip)
           icolo = PAR_CPU_TO_COLOR(ipart) % l(kcolo)
           number_parts(icolo) = number_parts(icolo) + 1
        end do
     end if
  end do
  do icolo = 0,mcolo
     if( number_parts(icolo) > 0 ) then
        call memory_alloca(par_memor,'PAR_COLOR_TO_CPU(IPART) % L','par_color_communicatros',PAR_COLOR_TO_CPU(icolo) % l,number_parts(icolo))
        number_parts(icolo) = 0
     else
        nullify( PAR_COLOR_TO_CPU(icolo) % l )        
     end if
  end do
  do ipart = 0,PAR_WORLD_SIZE-1
     if( associated(PAR_CPU_TO_COLOR(ipart) % l )) then 
        do kcolo = 1,size(PAR_CPU_TO_COLOR(ipart) % l,kind=ip)
           icolo = PAR_CPU_TO_COLOR(ipart) % l(kcolo)
           number_parts(icolo) = number_parts(icolo) + 1
           PAR_COLOR_TO_CPU(icolo) % l(number_parts(icolo)) = ipart
        end do
     end if
  end do
  if( associated(number_parts) ) deallocate( number_parts )
  !
  ! Up to now COMM_MATRIX is the inter-code matrix
  ! Take max over all codes to get the full one
  ! 
  call PAR_MAX(COMM_MATRIX,'IN THE WORLD')
  do icolo = 0,mcolo
     if( COMM_MATRIX(icolo) == 1 ) then
        do jcolo = icolo,mcolo
           if( COMM_MATRIX(jcolo) == 1 ) then
              PAR_COMM_COLOR(icolo,jcolo) = 1
              PAR_COMM_COLOR(jcolo,icolo) = 1
           end if
        end do
     end if
  end do
  deallocate( COMM_MATRIX )  

  !----------------------------------------------------------------------
  !
  ! Split world communicator to compute PAR_COMM_COLOR(ICOLO,JCOLO)
  ! Save permutation array PAR_COMM_COLOR_PERM(icolo,jcolo,PAR_MY_WORLD_RANK)
  !
  !----------------------------------------------------------------------
  !
  ! Inter-zone and inter-subd communicators only
  !
  do icolo = 0,mcolo
     do jcolo = icolo,mcolo
        if( PAR_COMM_COLOR(icolo,jcolo) == 1 ) then
           !
           ! Existing combinations: Should I be involved?
           !
           if( I_AM_IN_COLOR(icolo) .or. I_AM_IN_COLOR(jcolo) ) then
              icolor = 1
           else
              icolor = 0          
           end if
           isubd  = par_color_to_subd(icolo)
           izone  = par_color_to_zone(icolo)
           icode  = par_color_to_code(icolo)
           jsubd  = par_color_to_subd(jcolo)
           jzone  = par_color_to_zone(jcolo)
           jcode  = par_color_to_code(jcolo)
           isplit = .false.
           if( izone /= 0 .and. jzone /= 0 .and. isubd == 0 .and. jsubd == 0 ) isplit = .true.
           if( isubd /= 0 .and. jsubd /= 0 .and. izone == 0 .and. jzone == 0 ) isplit = .true.     
           if( icode /= 0 .and. jcode /= 0 .and. icode /= jcode .and. izone == 0 .and. jzone == 0 .and. isubd == 0 .and. jsubd == 0 ) isplit = .true.     
           if( isplit ) then
              call PAR_COMM_SPLIT(icolor,PAR_COMM_COLOR(icolo,jcolo),my_new_rank,'IN THE WORLD')
              if( my_new_rank /= -1 ) then
                 PAR_COMM_COLOR_PERM(icolo,jcolo,PAR_MY_WORLD_RANK) = my_new_rank
              else
                 PAR_COMM_COLOR(icolo,jcolo) = -1
                 PAR_COMM_COLOR_PERM(icolo,jcolo,PAR_MY_WORLD_RANK) = -1
              end if
           else
              PAR_COMM_COLOR(icolo,jcolo) = -1
              PAR_COMM_COLOR_PERM(icolo,jcolo,PAR_MY_WORLD_RANK) = -1
           end if
        else
           PAR_COMM_COLOR(icolo,jcolo) = -1
           PAR_COMM_COLOR_PERM(icolo,jcolo,PAR_MY_WORLD_RANK) = -1
        end if
        PAR_COMM_COLOR(jcolo,icolo) = PAR_COMM_COLOR(icolo,jcolo)
        PAR_COMM_COLOR_PERM(jcolo,icolo,PAR_MY_WORLD_RANK) = PAR_COMM_COLOR_PERM(icolo,jcolo,PAR_MY_WORLD_RANK)
     end do
  end do
  !
  ! COMM WORLD communicator
  !
  PAR_COMM_COLOR(0,0) = PAR_COMM_WORLD
  PAR_COMM_COLOR_PERM(0,0,PAR_MY_WORLD_RANK) = PAR_MY_WORLD_RANK
  !
  ! This code communicator
  ! At this point, PAR_COMM_COLOR_ARRAY is not fully usable. In fact, we miss
  ! the communicator for fringe geometry
  !
  icolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip)
  PAR_COMM_COLOR(icolo,icolo)                        = PAR_COMM_MY_CODE
  PAR_COMM_COLOR_PERM(icolo,icolo,PAR_MY_WORLD_RANK) = PAR_MY_CODE_RANK
  PAR_COMM_MY_CODE_ARRAY(1) % PAR_COMM_WORLD         = PAR_COMM_MY_CODE
  PAR_COMM_MY_CODE_ARRAY(1) % RANK4                  = int(PAR_MY_CODE_RANK,4)
  !call PAR_COPY_COMMUNICATION_ARRAY(PAR_COMM_MY_CODE_ARRAY,PAR_COMM_COLOR_ARRAY(icolo:))
  !
  ! Initialization of PAR_COMM_CURRENT
  !
  PAR_COMM_CURRENT = PAR_COMM_MY_CODE
  !
  ! Get the permutation of everybody
  !
  !  call PAR_SUM(PAR_COMM_COLOR_PERM,'IN THE UNIVERSE')
  call PAR_SUM_ALL(PAR_COMM_COLOR_PERM)
  !
  ! World to code rank
  !
  !allocate( recvcount4(0:PAR_WORLD_SIZE-1) )
  allocate( recvbuf(0:PAR_WORLD_SIZE*2-1) )
  allocate( sendbuf(2) )
  sendbuf(1) = current_code
  sendbuf(2) = int(PAR_MY_CODE_RANK,ip)
  !do ipart = 0,PAR_WORLD_SIZE-1
  !   recvcount4(ipart) = 2
  !end do
  recvcount4 = 2_4
  call PAR_ALLGATHER(sendbuf,recvbuf,recvcount4,'IN THE WORLD')
  ii = -1
  do ipart = 0,PAR_WORLD_SIZE-1
     ii = ii + 1
     PAR_COMM_WORLD_TO_CODE_PERM(1,ipart) = recvbuf(ii)
     ii = ii + 1
     PAR_COMM_WORLD_TO_CODE_PERM(2,ipart) = recvbuf(ii)
  end do
  deallocate( sendbuf )
  deallocate( recvbuf )
  !deallocate( recvcount4 )

  !----------------------------------------------------------------------
  !
  ! Check code number is not repeated over the different codes of the world
  !
  !----------------------------------------------------------------------

  ierro = 0
  if( PAR_MY_CODE_RANK == 0 ) then
     do ipart = 0,PAR_WORLD_SIZE-1
        if( ipart /= PAR_MY_WORLD_RANK ) then
           jcode = PAR_COMM_WORLD_TO_CODE_PERM(1,ipart)
           jpart = PAR_COMM_WORLD_TO_CODE_PERM(2,ipart)
           if( jpart == 0 .and. jcode == current_code ) then
              ierro = 1
           end if
        end if
     end do
  end if
  !  call PAR_MAX(ierro,'IN THE UNIVERSE')
  call PAR_MAX_ALL(ierro)
  if( ierro == 1 ) call runend('PAR_COLOR_COMMUNICATORS: WRONG CODE NUMBER')

  !----------------------------------------------------------------------
  !
  ! Output
  !
  !----------------------------------------------------------------------

  if( PAR_MY_CODE_RANK == 0 ) then
     ioutp(1) = mcode; ioutp(2) = mzone; ioutp(3) = msubd
     ioutp(4) = nzone; ioutp(5) = nsubd
     call outfor(66_ip,lun_outpu,'')
     do icolo = 0,mcolo
        do jcolo = icolo,mcolo
           if( PAR_COMM_COLOR(icolo,jcolo) /= -1 ) then
              isubd  = par_color_to_subd(icolo)
              izone  = par_color_to_zone(icolo)
              icode  = par_color_to_code(icolo)
              jsubd  = par_color_to_subd(jcolo)
              jzone  = par_color_to_zone(jcolo)
              jcode  = par_color_to_code(jcolo)
              ioutp(1) = icode; ioutp(2) = izone; ioutp(3) = isubd
              ioutp(4) = jcode; ioutp(5) = jzone; ioutp(6) = jsubd
              call outfor(67_ip,lun_outpu,'')
           end if
        end do
     end do
     call iofile_flush_unit(lun_outpu)
  end if

end subroutine par_color_communicators
