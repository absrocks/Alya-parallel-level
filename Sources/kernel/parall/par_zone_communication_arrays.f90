!------------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_domgra.f90
!> @author  Guillaume Houzeaux
!> @brief   Zone-wise communication array
!> @details Zone-wise communication array
!> @} 
!------------------------------------------------------------------------
subroutine par_zone_communication_arrays()
  use def_kintyp,         only :  ip,rp,comm_data_par
  use def_master,         only :  ISEQUEN,ISLAVE
  use def_master,         only :  current_code
  use def_master,         only :  current_zone
  use def_master,         only :  current_subd
  use def_master,         only :  gisca,npoi3
  use def_domain,         only :  npoin,nzone
  use def_domain,         only :  nelem,lesub,nsubd,lnods
  use def_domain,         only :  lnnod,lperi,nperi
  use mod_communications, only :  PAR_DEFINE_COMMUNICATOR
  use mod_communications, only :  PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only :  PAR_SEND_RECEIVE
  use mod_communications, only :  PAR_MAX
  use mod_communications, only :  PAR_COMM_RANK_AND_SIZE
  use mod_parall,         only :  PAR_COMM_COLOR
  use mod_parall,         only :  PAR_COMM_COLOR_ARRAY
  use mod_parall,         only :  PAR_COMM_MY_CODE_ARRAY 
  use mod_parall,         only :  PAR_COMM_COLOR_PERM
  use mod_parall,         only :  par_world_rank_of_a_code_neighbor
  use mod_parall,         only :  PAR_MY_WORLD_RANK
  use mod_parall,         only :  PAR_THIS_NODE_IS_MINE
  use mod_parall,         only :  PAR_COPY_COMMUNICATION_ARRAY 
  use mod_parall,         only :  par_code_zone_subd_to_color
  use mod_parall,         only :  PAR_SUB_COMMUNICATION_ARRAY
  use mod_messages,       only :  messages_live
  implicit none
  integer(ip)                  :: PAR_COMM_ORIGINAL
  integer(4)                   :: PAR_COMM_ORIGINAL4
  integer(4)                   :: PAR_RANK_ZONE
  integer(ip)                  :: ierro,dom_i_world,dom_i_zone,ielem
  integer(ip)                  :: icolo,ipoin,ineig,dom_i_code
  integer(ip)                  :: inode,jpoin,iperi
  type(comm_data_par), pointer :: commu

  if( ISEQUEN ) return

  call messages_live('PARALL: CREATE INTRA-ZONE COMMUNICATION ARRAY')
  ! 
  ! Allocate memory
  !
  call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_ORIGINAL4,commu)
  PAR_COMM_ORIGINAL = int(PAR_COMM_ORIGINAL4,ip)
  
  !----------------------------------------------------------------------
  !
  ! All-zone and all-subdomain communication array
  !
  !----------------------------------------------------------------------

  icolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip)   
  call PAR_COPY_COMMUNICATION_ARRAY(PAR_COMM_MY_CODE_ARRAY(1),PAR_COMM_COLOR_ARRAY(icolo),COMM_NAME='PAR_COMM_COLOR_ARRAY')

  !----------------------------------------------------------------------
  !
  ! Intra-zone communication arrays
  !
  !----------------------------------------------------------------------

  if( ISLAVE ) call memgen(1_ip,npoin,0_ip)

  do current_zone = 1,nzone

     icolo = par_code_zone_subd_to_color(current_code,current_zone,0_ip)
     PAR_RANK_ZONE = PAR_COMM_COLOR_PERM(icolo,icolo,PAR_MY_WORLD_RANK)
     !
     ! Sub-communicator only in current zone
     !            
     if( ISLAVE ) then

        do ipoin = 1,npoin
           gisca(ipoin) = 1
        end do
        call PAR_SUB_COMMUNICATION_ARRAY(commu,PAR_COMM_COLOR_ARRAY(icolo),gisca,COMM_NAME='PAR_COMM_COLOR_ARRAY')
        !
        ! Convert code neighbor to zone neighbor (dm_i_code is kfl_paral!!!)
        ! This is because the communicator was already defined using a split
        ! at the beginning (PAR_COMM_COLOR_ARRAY(icolo) % PAR_COMM_WORLD)
        ! and commu invovles the code partition numbering
        !
        do ineig = 1,PAR_COMM_COLOR_ARRAY(icolo) % nneig
           dom_i_code  = PAR_COMM_COLOR_ARRAY(icolo) % neights(ineig)
           dom_i_world = par_world_rank_of_a_code_neighbor(dom_i_code,current_code)
           dom_i_zone  = PAR_COMM_COLOR_PERM(icolo,icolo,dom_i_world)
           PAR_COMM_COLOR_ARRAY(icolo) % neights(ineig) = dom_i_zone
        end do
        do ipoin = 1,npoin
           gisca(ipoin) = 0
        end do
     end if
     !
     ! MPI communicator
     !
     PAR_COMM_COLOR_ARRAY(icolo) % PAR_COMM_WORLD = PAR_COMM_COLOR(icolo,icolo)

  end do

  if( ISLAVE ) call memgen(3_ip,npoin,0_ip)

  !----------------------------------------------------------------------
  !
  ! A short test to check inter-zone communication
  !
  !----------------------------------------------------------------------

  do current_zone = 1,nzone
     ierro = 0
     if( ISLAVE ) then
        call memgen(1_ip,npoin,0_ip)
        do ipoin = 1,npoi3
           gisca(ipoin) = 1
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM','IN MY ZONE','SYNCHRONOUS')        
        do ipoin = 1,npoin
           if( gisca(ipoin) /= 1 ) then
              ierro = ierro + 1
           end if
        end do
        call memgen(3_ip,npoin,0_ip)
     end if
     call PAR_MAX(ierro,'IN MY CODE')
     !if( ierro /= 0 ) call runend('PROBLEM WITH ZONE COMMUNICATION ARRAY')
  end do

  !----------------------------------------------------------------------
  !
  ! Very important: CURRENT_ZONE is all zones
  !
  !----------------------------------------------------------------------

  current_zone = 0

  !---------------------------------------------------------------------------------------------------!
  !
  !  ==========================
  !  DO THE SAME FOR SUBDOMAINS
  !  ==========================
  !
  !---------------------------------------------------------------------------------------------------!

  call messages_live('PARALL: CREATE INTRA-SUBDOMAIN COMMUNICATION ARRAY')
  !
  ! Allocate memory
  !
  call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_ORIGINAL4,commu)
  PAR_COMM_ORIGINAL = int(PAR_COMM_ORIGINAL4,ip)

  !----------------------------------------------------------------------
  !
  ! Inter-zone communication arrays
  !
  !----------------------------------------------------------------------

  if( ISLAVE ) call memgen(1_ip,npoin,0_ip)

  do current_subd = 1,nsubd

     icolo = par_code_zone_subd_to_color(current_code,0_ip,current_subd)
     !
     ! Sub-communicator only in current zone
     !            
     if( ISLAVE ) then

        do ielem = 1,nelem
           if( lesub(ielem) == current_subd ) then
              do inode = 1,lnnod(ielem)
                 ipoin = lnods(inode,ielem)
                 gisca(ipoin) = 1
              end do
           end if
        end do
        do iperi = 1,nperi
           ipoin = lperi(1,iperi) 
           jpoin = lperi(2,iperi)
           gisca(ipoin) = gisca(jpoin) 
        end do

        call PAR_INTERFACE_NODE_EXCHANGE(gisca,'MAX','IN MY CODE')
        call PAR_SUB_COMMUNICATION_ARRAY(commu,PAR_COMM_COLOR_ARRAY(icolo),gisca,COMM_NAME='PAR_COMM_COLOR_ARRAY')
        !
        ! Convert code neighbor to zone neighbor (dm_i_code is kfl_paral!!!)
        ! This is because the communicator was already defined using a split
        ! at the beginning (PAR_COMM_COLOR_ARRAY(icolo) % PAR_COMM_WORLD)
        ! and commu invovles the code partition numbering
        !
        do ineig = 1,PAR_COMM_COLOR_ARRAY(icolo) % nneig
           dom_i_code  = PAR_COMM_COLOR_ARRAY(icolo) % neights(ineig)
           dom_i_world = par_world_rank_of_a_code_neighbor(dom_i_code,current_code)
           dom_i_zone  = PAR_COMM_COLOR_PERM(icolo,icolo,dom_i_world)
           PAR_COMM_COLOR_ARRAY(icolo) % neights(ineig) = dom_i_zone
        end do
        do ipoin = 1,npoin
           gisca(ipoin) = 0
        end do
     end if
     !
     ! MPI communicator
     !
     PAR_COMM_COLOR_ARRAY(icolo) % PAR_COMM_WORLD = PAR_COMM_COLOR(icolo,icolo)

  end do

  if( ISLAVE ) call memgen(3_ip,npoin,0_ip)
  
  !----------------------------------------------------------------------
  !
  ! A short test to check intra-subdomain communication
  !
  !----------------------------------------------------------------------
  
  do current_subd = 1,nsubd
     ierro = 0
     if( ISLAVE ) then
        call memgen(1_ip,npoin,0_ip)
 
        do ielem = 1,nelem
           if( lesub(ielem) == current_subd ) then
              do inode = 1,lnnod(ielem)
                 ipoin = lnods(inode,ielem)
                 if( PAR_THIS_NODE_IS_MINE(ipoin) ) then
                    gisca(ipoin) =  1
                 else
                    gisca(ipoin) = -1
                 end if
              end do
           end if
        end do

        do iperi = 1,nperi
           ipoin = lperi(1,iperi)
           jpoin = lperi(2,iperi)
           if( PAR_THIS_NODE_IS_MINE(ipoin) .and. abs(gisca(jpoin)) == 1 ) gisca(ipoin) = 1
        end do

        if( npoin > 0 ) gisca = max(0_ip,gisca)
        
        call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM','IN MY SUBD','SYNCHRONOUS')        
        
        do ielem = 1,nelem
           if( lesub(ielem) == current_subd ) then
              do inode = 1,lnnod(ielem)
                 ipoin = lnods(inode,ielem)
                 if( gisca(ipoin) /= 1 ) then
                    ierro = ierro + 1
                 end if
              end do
           end if
        end do
        call memgen(3_ip,npoin,0_ip)
     end if
     !call PAR_MAX(ierro,'IN MY CODE')
     !if( ierro /= 0 ) call runend('PROBLEM WITH SUBDOMAIN COMMUNICATION ARRAY')
  end do

  !----------------------------------------------------------------------
  !
  ! Very important: CURRENT_SUBD is all subdomains
  !
  !----------------------------------------------------------------------

  current_subd = 0

end subroutine par_zone_communication_arrays


subroutine PAR_SUB_COMMUNICATION_ARRAY_1(COMM_IN,COMM_OUT,mask,COMM_NAME)

  use def_kintyp,         only : ip,rp,lg,i1p,comm_data_par,ompss_domain
  use def_domain,         only : npoin,nelem,nboun,npoin_2
  use def_domain,         only : mesh_type,npoin
  use def_domain,         only : htable_lninv_loc
  use mod_graphs,         only : graphs_number_to_linked_list
  use def_master,         only : ISEQUEN,ISLAVE,IMASTER,INOTMASTER
  use def_master,         only : npoi1,npoi2,npoi3,lninv_loc
  use def_master,         only : intost,gisca,igene,lun_outpu
  use def_master,         only : ioutp,kfl_paral
  use mod_memory,         only : memory_copy
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_size
  use mod_graphs,         only : graphs_eleele
  use mod_graphs,         only : graphs_dealep
  use mod_graphs,         only : graphs_coloring
  use mod_graphs,         only : graphs_coloring_greedy
  use mod_graphs,         only : graphs_deallocate
  use mod_htable,         only : htalid
   use mod_parall
  implicit none
    type(comm_data_par),          intent(in)    :: COMM_IN        !< Input communicator
    type(comm_data_par),          intent(inout) :: COMM_OUT       !< Output communicator
    integer(ip),                  intent(in)    :: mask(*)        !< mask (=1 to consider node)
    character(*),                 intent(in)    :: COMM_NAME
    integer(ip),         pointer                :: bound_perm(:)
    integer(ip),         pointer                :: bound_invp(:)
    integer(ip),         pointer                :: bound_size(:)
    integer(ip),         pointer                :: neights(:)
    integer(ip)                                 :: nneig,bound_dim
    integer(ip)                                 :: nneig_1,nneig_2
    integer(ip)                                 :: ineig,jj,jneig
    integer(ip)                                 :: ipoin
    character(50)                               :: my_comm_name
    
    COMM_OUT % PAR_COMM_WORLD = COMM_IN % PAR_COMM_WORLD
    COMM_OUT % RANK4          = COMM_IN % RANK4
    
       my_comm_name = trim(COMM_NAME)

    if( COMM_IN % bound_dim > 0 ) then

       nullify(bound_perm)
       nullify(bound_invp)
       nullify(bound_size)
       nullify(neights)

       allocate( bound_perm(COMM_IN % bound_dim) )
       allocate( bound_invp(COMM_IN % bound_dim) )
       allocate( bound_size(COMM_IN % nneig)     )
       allocate( neights   (COMM_IN % nneig)     )

       nneig     = 0
       bound_dim = 0

       do jj = 1,COMM_IN % bound_dim
          bound_perm(jj) = COMM_IN % bound_perm(jj)
          bound_invp(jj) = COMM_IN % bound_invp(jj)
       end do
       do ineig = 1,COMM_IN % nneig
          bound_size(ineig) = 0
       end do

       do ineig = 1,COMM_IN % nneig
          do jj = COMM_IN % bound_size(ineig),COMM_IN % bound_size(ineig+1)-1
             ipoin = COMM_IN % bound_perm(jj)
             if( mask(ipoin) > 0 ) then
                bound_size(ineig) = bound_size(ineig) + 1
             else
                bound_perm(jj) = 0
                bound_invp(jj) = 0
             end if
          end do
          bound_dim = bound_dim + bound_size(ineig)

          if( bound_size(ineig) > 0 ) nneig = nneig + 1
       end do
       !
       ! Allocate interzone communication array
       !
       COMM_OUT % nneig     = nneig
       COMM_OUT % bound_dim = bound_dim
       COMM_OUT % npoi1     = COMM_IN % npoi1
       COMM_OUT % npoi2     = COMM_IN % npoi2
       COMM_OUT % npoi3     = COMM_IN % npoi3
       COMM_OUT % npoin     = COMM_IN % npoin
       call memory_alloca(par_memor,trim(my_comm_name)//' % NEIGHTS'        ,'PAR_SUB_COMMUNICATION_ARRAY',COMM_OUT % neights,nneig)
       call memory_alloca(par_memor,trim(my_comm_name)//' % NEIGHTS_ORDERED','PAR_SUB_COMMUNICATION_ARRAY',COMM_OUT % neights_ordered,nneig)
       call memory_alloca(par_memor,trim(my_comm_name)//' % PERM_ORDERED'   ,'PAR_SUB_COMMUNICATION_ARRAY',COMM_OUT % perm_ordered,nneig)
       call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_PERM'     ,'PAR_SUB_COMMUNICATION_ARRAY',COMM_OUT % bound_perm,bound_dim)
       call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_INVP'     ,'PAR_SUB_COMMUNICATION_ARRAY',COMM_OUT % bound_invp,bound_dim)
       call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_SIZE'     ,'PAR_SUB_COMMUNICATION_ARRAY',COMM_OUT % bound_size,nneig+1)

       !allocate( COMM_OUT % neights(nneig) )
       !allocate( COMM_OUT % bound_perm(bound_dim) )
       !allocate( COMM_OUT % bound_size(nneig+1) )
       !
       ! Permutation array BOUND_PERM(1:BOUND_DIM)
       !
       jneig = 0
       bound_dim = 0
       do ineig = 1,COMM_IN % nneig
          if( bound_size(ineig) > 0 ) then

             jneig = jneig + 1
             COMM_OUT % neights(jneig)         = COMM_IN % neights(ineig)
             COMM_OUT % neights_ordered(jneig) = COMM_IN % neights_ordered(ineig)
             COMM_OUT % perm_ordered(jneig)    = COMM_IN % perm_ordered(ineig)

             do jj = COMM_IN % bound_size(ineig),COMM_IN % bound_size(ineig+1)-1
                ipoin = bound_perm(jj)
                if( ipoin > 0 ) then
                   bound_dim = bound_dim + 1
                   COMM_OUT % bound_perm(bound_dim) = ipoin
                   COMM_OUT % bound_invp(bound_dim) = ipoin
                end if
             end do

          end if
       end do
       ineig_1: do ineig = 1,COMM_IN % nneig
          if( ineig > int(COMM_IN % RANK4,ip) ) then
             nneig_1 = ineig-1
             nneig_2 = ineig
             exit ineig_1
          end if
       end do ineig_1
       !
       ! Construct linked list BOUND_SIZE(1:NNEIG+1)
       !
       COMM_OUT % bound_size(1) = 1
       jneig = 0
       do ineig = 1,COMM_IN % nneig
          if( bound_size(ineig) > 0 ) then
             jneig = jneig + 1
             COMM_OUT % bound_size(jneig+1) = &
                  COMM_OUT % bound_size(jneig) + bound_size(ineig)
          end if
       end do

       if( associated(bound_perm) ) deallocate( bound_perm )
       if( associated(bound_invp) ) deallocate( bound_invp )
       if( associated(bound_size) ) deallocate( bound_size )
       if( associated(neights)    ) deallocate( neights    )

    end if

  end subroutine PAR_SUB_COMMUNICATION_ARRAY_1
