!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_iniunk.f90
!> @author  Mariano Vazquez
!> @date    26/12/2016
!> @brief   Starting: This routine sets up the initial conditions
!> @details Starting: This routine sets up the initial conditions
!> @}
!-----------------------------------------------------------------------
subroutine exm_iniunk
   use def_master
   use def_domain
   use mod_iofile
   use def_exmedi
   use mod_maths,          only : maths_distribute_iterations
   use mod_communications, only : PAR_BROADCAST
   use mod_messages,       only : messages_live
   use mod_memory,         only : memory_initia
   use mod_memory,         only : memory_alloca
   use mod_memory,         only : memory_deallo
   use mod_exm_onecel
   use mod_exm_oneohr
   use mod_exm_onetor
   use mod_outfor
   use mod_exm_oceohr,      only : exm_oceohr, exm_oceola
   use mod_exm_sld_eccoupling, only: EXMSLD_EMEC_LAND, EXMSLD_EMEC_LAND_BIDIR
   use mod_exm_sld_eccoupling, only: cell_ca0_ecc
   use mod_exm_sld_eccoupling, only : kfl_exmsld_3Dcou_ecc
   use mod_exm_sld_eccoupling, only: state_ecc, troponin_ecc, nsdvs_ecc, cell_ca0_ecc 

   implicit none

   integer(ip), parameter               :: ohara_nstats = 3_ip       !number of values returned in ohara_stats
   integer(ip), parameter               :: torord_nstats = 3_ip


   integer(ip)                          :: kmodel_imate, kmodel_ipoin, imate
   integer(ip)                          :: ipoin, icelltype, iconc, i, j
   integer(ip)                          :: oneohr_status_tmp,onetor_status_tmp
   integer(ip)                          :: ipair, npairs,root_rank
   character(30)                        :: mesau(2)
   integer(ip), dimension(:,:), pointer :: mater_cell_pairs
   integer(ip), dimension(:),   pointer :: oneohr_status,onetor_status
   integer(ip), dimension(:),   pointer :: list_pairs
   real(rp),    dimension(7)            :: land_variables ! used temporarily
   real(rp),    dimension(7,mnode)      :: eland
   real(rp),    dimension(:,:), pointer :: land_variables_allmaterials
   integer(ip)                          :: igaus, inode, pnode, pgaus, ielem, pelty


   real(rp), dimension(:,:), pointer    :: ohara_stats  ! saves ohara stats: number of beats, tolerance, rmse
   integer(ip)                          :: ohara_nbeats ! used temporarily for printing
   real(rp)                             :: ohara_toler, ohara_rmse  ! used temporarily for printing
   character(10)                        :: ohara_stats_variablename !used for output
   logical(lg)                          :: has_land

   real(rp), dimension(:,:), pointer    :: torord_stats  ! saves ohara stats: number of beats, tolerance, rmse
   integer(ip)                          :: torord_nbeats ! used temporarily for printing
   real(rp)                             :: torord_toler, torord_rmse  ! used temporarily for printing
   character(10)                        :: torord_stats_variablename !used for output

   nullify(mater_cell_pairs,oneohr_status,onetor_status,list_pairs)
   nullify(land_variables_allmaterials)

   nullify(ohara_stats)
   nullify(torord_stats)
   call memory_initia(ohara_stats)
   call memory_initia(torord_stats)

   if( kfl_rstar == 0_ip ) then
      !
      ! Initialize
      !
      call memory_initia(elmag)
      call memory_initia(vconc)
      !
      ! check if any TT-like model is present
      !
      vcoin_exm    = 0.0_rp
      vauin_exm    = 0.0_rp
      kmodel_imate = 0_ip
      kmodel_ipoin = 0_ip

      do imate= 1,nmate
         !
         ! message: compute initial conditions for some cell models
         !
         kmodel_imate = kfl_cellmod(imate)

         if (kmodel_imate == CELL_TT2006_EXMEDI .or. kmodel_imate == CELL_OHARA_EXMEDI .or. kmodel_imate == CELL_TORORD_EXMEDI) then
            call messages_live('  COMPUTING INITIAL PHYSIOLOGICAL CONDITIONS...')
            !if (kmodel_imate == 4 .or. kmodel_imate==5 .or. kmodel_imate==6 ) then
            call messages_live('    ITERATING ONE CELL MODEL...')
            mesau(1)= intost(moneclmate_exm(1,imate)) ! beats
            mesau(2)= intost(moneclmate_exm(2,imate)) ! cyclelength
            call messages_live('    BEATS=       '//trim(mesau(1)))
            call messages_live('    CYCLELENGTH= '//trim(mesau(2)))
         end if

         call messages_live('  INITIAL CONDITIONS FOR MATERIAL= '//trim(intost(imate)))
         if (kmodel_imate == CELL_TT2006_EXMEDI) then
            call messages_live('    SUB CELL MODEL: TT HETEROGENEOUS')
         else if (kmodel_imate == CELL_OHARA_EXMEDI) then
            call messages_live('    SUB CELL MODEL: OHARA - RUDY')
         else if (kmodel_imate == CELL_TORORD_EXMEDI) then
            call messages_live('    SUB CELL MODEL: TOR-ORD')
         end if

         select case ( kfl_steadystate_variable(imate) )
            case ( EXM_CELL_STEADY_VOLTAGE )
               call messages_live("    USING VOLTAGE TO TEST CELL CONVERGENCE FOR MATERIAL.")
            case ( EXM_CELL_STEADY_CALCIUM )
               call messages_live("    USING CALCIUM TO TEST CELL CONVERGENCE FOR MATERIAL.")
            case default
               call runend("EXM_INIUNK: Unknown name of the variable to determine the cell convergence.")
         end select


      end do


      npairs = nmate*EXM_CELLTYPE_MAXID

      call memory_alloca(mem_modul(1:2,modul),'MATER_CELL_PAIRS','exm_iniunk', mater_cell_pairs, npairs, 2_ip)
      call memory_alloca(mem_modul(1:2,modul),'ONEOHR_STATUS'   ,'exm_iniunk', oneohr_status, npairs)
      call memory_alloca(mem_modul(1:2,modul),'ONETOR_STATUS'   ,'exm_iniunk', onetor_status, npairs)
      call memory_alloca(mem_modul(1:2,modul),'OHARA_STATS'     ,'exm_iniunk', ohara_stats, npairs, ohara_nstats)
      call memory_alloca(mem_modul(1:2,modul),'TORORD_STATS'    ,'exm_iniunk',torord_stats, npairs, torord_nstats)

      !save material and cell type to a linear array for paralelization
      ipair = 1_ip
      do imate= 1,nmate
         do icelltype = 1,EXM_CELLTYPE_MAXID
            mater_cell_pairs(ipair, 1_ip) = imate
            mater_cell_pairs(ipair, 2_ip) = icelltype
            oneohr_status(ipair) = -1_ip
            onetor_status(ipair) = -1_ip
            ipair = ipair + 1_ip
         end do
      end do
      !
      ! LPAIRS(IPAIR) is the list of rank owner of pair IPAIR
      !
      call memory_alloca(mem_modul(1:2,modul),'LIST_PAIRS','exm_iniunk',list_pairs,npairs)

      !this one is npairs x 6
      call memory_alloca(mem_modul(1:2,modul),'LAND_VARIABLES','exm_iniunk',land_variables_allmaterials,npairs, size(land_variables,KIND=ip))

      land_variables_allmaterials=0.0_rp
      land_variables=0.0_rp


      call maths_distribute_iterations(npart,npairs,list_pairs,DEFAULT_VALUE=kfl_paral)

      ! Calculate initial conditions for the cell model
      do ipair= 1,npairs

         if( list_pairs(ipair) == kfl_paral ) then

            imate = mater_cell_pairs(ipair, 1_ip)
            icelltype = mater_cell_pairs(ipair, 2_ip)

            kmodel_imate = kfl_cellmod(imate)

            if (kmodel_imate == CELL_TT2006_EXMEDI) then
               !this TT model needs deep revision, it may not be working correctly
               !maybe will be removed in the future
               !modifies: vminimate_exm(icelltype,imate), vauin_exm(:,icelltype,imate), vcoin_exm(:,icelltype,imate)
               call exm_onecel(imate) ! these are precalculated, just sets up the variables runs extremeley fast. Needs to be refactored to pass also celltype


            else if (kmodel_imate == CELL_OHARA_EXMEDI) then

               !
               ! Initialization of the cell model: done only by the master and then distributed to the slaves
               !
               !modifies Global vars: vminimate_exm(icelltype,imate), vauin_exm(:,icelltype,imate), vcoin_exm(:,icelltype,imate)
               !Modifies Local: elmlo_exm, viclo_exm, vcolo_exm, vaulo_exm
               !Modifies: land_variables only is land model is used and there is coupling with solidz
               call exm_oneohr(imate, icelltype, land_variables, oneohr_status_tmp, ohara_stats(ipair,:) )

               land_variables_allmaterials(ipair,:) = land_variables

               oneohr_status(ipair) = oneohr_status_tmp !die if oneohr_status != 0

            else if (kmodel_imate == CELL_TORORD_EXMEDI) then

                !
                ! Initialization of the cell model: done only by the master then distributed to the slaves
                !
                call exm_onetor(imate, icelltype, land_variables, onetor_status_tmp, torord_stats(ipair,:) )
                land_variables_allmaterials(ipair,:) = land_variables

                onetor_status(ipair) = onetor_status_tmp

            end if

         end if
      end do !ipair

      !Send oneohr_status from children to master

      !
      ! Owner of the pair broadcast the result to others
      !
      do ipair = 1,npairs
         root_rank = list_pairs(ipair)
         call PAR_BROADCAST(oneohr_status(ipair)     , ROOT_RANK=root_rank)
         call PAR_BROADCAST(onetor_status(ipair)     , ROOT_RANK=root_rank)
         call PAR_BROADCAST(mater_cell_pairs(ipair,1), ROOT_RANK=root_rank)
         call PAR_BROADCAST(mater_cell_pairs(ipair,2), ROOT_RANK=root_rank)

         do i = 1, ohara_nstats
            call PAR_BROADCAST(ohara_stats(ipair,i), ROOT_RANK=root_rank)
         end do

         do i = 1, torord_nstats
            call PAR_BROADCAST(torord_stats(ipair,i), ROOT_RANK=root_rank)
         end do

         do i = 1, size(land_variables_allmaterials, 2, KIND=ip)
            call PAR_BROADCAST(land_variables_allmaterials(ipair,i), ROOT_RANK=root_rank)
         end do
      end do

      do imate = 1,nmate
         do icelltype = 1,EXM_CELLTYPE_MAXID
            ipair     = (imate-1)*EXM_CELLTYPE_MAXID + icelltype
            root_rank = list_pairs(ipair)

            call PAR_BROADCAST(vminimate_exm(icelltype,imate), ROOT_RANK=root_rank)
            do iconc = 1,nconc_exm
               call PAR_BROADCAST(vcoin_exm(iconc, icelltype, imate), ROOT_RANK=root_rank)
            end do
            do iconc = 1,nauxi_exm
               call PAR_BROADCAST(vauin_exm(iconc, icelltype, imate), ROOT_RANK=root_rank)
            end do

         end do
      end do


      if (INOTSLAVE) then
         !Master should print the information about the cell model convergence
         call outfor(25_ip,momod(modul) % lun_outpu,'OHARA CONVERGENCE')


         do ipair= 1,npairs
            imate = mater_cell_pairs(ipair, 1_ip)
            icelltype = mater_cell_pairs(ipair, 2_ip)

            ohara_nbeats = int(ohara_stats(ipair, 1_ip))
            ohara_toler = ohara_stats(ipair, 2_ip)
            ohara_rmse = ohara_stats(ipair, 3_ip)

            torord_nbeats = int(torord_stats(ipair, 1_ip))
            torord_toler = torord_stats(ipair, 2_ip)
            torord_rmse = torord_stats(ipair, 3_ip)

            select case (kfl_steadystate_variable(imate))
            case (EXM_CELL_STEADY_VOLTAGE)
               ohara_stats_variablename = 'VOLTAGE'
               torord_stats_variablename = 'VOLTAGE'
            case (EXM_CELL_STEADY_CALCIUM)
               ohara_stats_variablename = 'CALCIUM'
               torord_stats_variablename = 'CALCIUM'
            case default
               call runend("EXM_OCEOHR or EXM_OCETOR: Unknown name of the variable to determine the cell convergence.")
            end select

            if (oneohr_status(ipair) >= 0_ip) then !if negative -- not executed
               if ( oneohr_status(ipair) == EXM_ONEOHR_ERROR_NOTCONVERGED ) then
                  call messages_live("EXM_OCEOHR: CELL MODEL FOR MATERIAL "//trim(intost(imate))//" AND CELLTYPE "//trim(intost(icelltype))//" DID NOT CONVERGE. BEATS="//trim(intost(ohara_nbeats))//", TOL="//trim(retost(ohara_toler))//", RMSE="//trim(retost(ohara_rmse)),"WARNING")
                  write(momod(modul)%lun_outpu,*) "OHARA "//trim(ohara_stats_variablename)//" NOT CONVERGED, MAT="//trim(intost(imate))//" CELLTYPE="//trim(intost(icelltype))//" BEATS="//trim(intost(ohara_nbeats))//" TOL="//trim(retost(ohara_toler))//" RMSE="//trim(retost(ohara_rmse))

                  if( kfl_ignore_steadystate_celltypes_exm(imate, icelltype) == 0_ip ) then
                     call runend("EXM_OCEOHR: STEADY STATE NOT REACHED AFTER "//trim(intost(ohara_nbeats))//" BEAT, TOL="//trim(retost(ohara_toler))//", RMSE="//trim(retost(ohara_rmse))//". TRY INCREASING THE NUMBER OF BEATS")
                  end if
               else if ( oneohr_status(ipair) == EXM_ONEOHR_ERROR_CONVERGED ) then
                  ! DO NOT ERASE THIS MESSAGE
                  call messages_live("EXM_OCEOHR: CELL MODEL FOR MATERIAL "//trim(intost(imate))//" AND CELLTYPE "//trim(intost(icelltype))//" !!!DID!!! CONVERGE. BEATS="//trim(intost(ohara_nbeats))//", TOL="//trim(retost(ohara_toler))//", RMSE="//trim(retost(ohara_rmse)) )
                  write(momod(modul)%lun_outpu,*) "OHARA "//trim(ohara_stats_variablename)//" CONVERGED, MAT="//trim(intost(imate))//" CELLTYPE="//trim(intost(icelltype))//" BEATS="//trim(intost(ohara_nbeats))//" TOL="//trim(retost(ohara_toler))//" RMSE="//trim(retost(ohara_rmse))
               else if ( oneohr_status(ipair) == EXM_ONEOHR_ERROR_NOTINITIALIZED ) then
                  call messages_live("EXM_ONEOHR: MATERIAL "//trim(intost(imate))//" CELLTYPE "//trim(intost(icelltype))//", VOLTAGES WERE NOT INITIALIZED. UNKNOWN MYOCITE TYPE. MOST LIKELY REQUESTED NONEXISTING HARDCODED CELLTYPE.","WARNING")
                  write(momod(modul)%lun_outpu,*) "OHARA NOT INITIALIZED, MAT "//trim(intost(imate))//" CELLTYPE "//trim(intost(icelltype))//" TOL "//trim(retost(ohara_toler))//" RMSE="//trim(retost(ohara_rmse))
               else
                  call runend("EXM_ONEOHR: UNHANDLED RETURN VALUE OF EXM_ONEOHR - "//trim(intost(oneohr_status(ipair))))
               end if
            end if

            if (onetor_status(ipair) >= 0_ip) then !if negative -- not executed
               if ( onetor_status(ipair) == EXM_ONETOR_ERROR_NOTCONVERGED ) then
                  call messages_live("EXM_OCETOR: CELL MODEL FOR MATERIAL "//trim(intost(imate))//" AND CELLTYPE "//trim(intost(icelltype))//" DID NOT CONVERGE. BEATS="//trim(intost(torord_nbeats))//", TOL="//trim(retost(torord_toler))//", RMSE="//trim(retost(torord_rmse)),"WARNING")
                  write(momod(modul)%lun_outpu,*) "TORORD "//trim(torord_stats_variablename)//" NOT CONVERGED, MAT="//trim(intost(imate))//" CELLTYPE="//trim(intost(icelltype))//" BEATS="//trim(intost(torord_nbeats))//" TOL="//trim(retost(torord_toler))//" RMSE="//trim(retost(torord_rmse))

                  if( kfl_ignore_steadystate_celltypes_exm(imate, icelltype) == 0_ip ) then
                     call runend("EXM_OCETOR: STEADY STATE NOT REACHED AFTER "//trim(intost(torord_nbeats))//" BEAT, TOL="//trim(retost(torord_toler))//", RMSE="//trim(retost(torord_rmse))//". TRY INCREASING THE NUMBER OF BEATS")
                  end if
               else if ( onetor_status(ipair) == EXM_ONETOR_ERROR_CONVERGED ) then
                  ! DO NOT ERASE THIS MESSAGE
                  call messages_live("EXM_OCETOR: CELL MODEL FOR MATERIAL "//trim(intost(imate))//" AND CELLTYPE "//trim(intost(icelltype))//" !!!DID!!! CONVERGE. BEATS="//trim(intost(torord_nbeats))//", TOL="//trim(retost(torord_toler))//", RMSE="//trim(retost(torord_rmse)) )
                  write(momod(modul)%lun_outpu,*) "TORORD "//trim(torord_stats_variablename)//" CONVERGED, MAT="//trim(intost(imate))//" CELLTYPE="//trim(intost(icelltype))//" BEATS="//trim(intost(torord_nbeats))//" TOL="//trim(retost(torord_toler))//" RMSE="//trim(retost(torord_rmse))
               else if ( onetor_status(ipair) == EXM_ONETOR_ERROR_NOTINITIALIZED ) then
                  call messages_live("EXM_ONETOR: MATERIAL "//trim(intost(imate))//" CELLTYPE "//trim(intost(icelltype))//", VOLTAGES WERE NOT INITIALIZED. UNKNOWN MYOCITE TYPE. MOST LIKELY REQUESTED NONEXISTING HARDCODED CELLTYPE.","WARNING")
                  write(momod(modul)%lun_outpu,*) "TORORD NOT INITIALIZED, MAT "//trim(intost(imate))//" CELLTYPE "//trim(intost(icelltype))//" TOL "//trim(retost(torord_toler))//" RMSE="//trim(retost(torord_rmse))
               else
                  call runend("EXM_ONETOR: UNHANDLED RETURN VALUE OF EXM_ONETOR - "//trim(intost(onetor_status(ipair))))
               end if
            end if
         end do
      end if ! INOTSLAVE

      !--------------------------------------------------------------
      !
      ! Distribute land parameters among cells and points
      !
      !--------------------------------------------------------------     
      if( (kfl_coupl(ID_SOLIDZ,ID_EXMEDI) >= 1_ip .or. kfl_coupl(ID_EXMEDI,ID_SOLIDZ) >=1_ip) .or. kfl_exmsld_3Dcou_ecc ) then
            has_land=.False.
            do imate = 1,nmate
               ! Save Ca0 values to automatically calibrate Hunter-McCulloch ECC model in sld_eccoup
               cell_ca0_ecc(:,imate) = vcoin_exm(1,:,imate)

               if( kfl_eccty(imate) == EXMSLD_EMEC_LAND .or. kfl_eccty(imate) == EXMSLD_EMEC_LAND_BIDIR ) then
                   has_land=.True.
               endif
            end do      

            if(has_land)then
               call messages_live("EXM_INIUNK: PASSING LAND PARAMETERS TO SOLIDZ") 
               do ipoin = 1,npoin
                  imate = nodemat(ipoin)
                  icelltype = int(celty_exm(1,ipoin))            
                  ipair     = (imate-1)*EXM_CELLTYPE_MAXID + icelltype            

                  land_variables = land_variables_allmaterials(ipair,:)

                  troponin_ecc(ipoin,1) = land_variables(3) ! CaTRPN
                  troponin_ecc(ipoin,2) = land_variables(3) ! CaTRPN
               end do

               !interpolate land values for each gauss point from points
               do ielem = 1,nelem
                  !
                  ! Get element shape function
                  !
                  pelty = ltype(ielem)
                  if( pelty > 0 ) then
                     pgaus = ngaus(pelty)
                     pnode = lnnod(ielem)
                     
                     do inode = 1,pnode
                        ipoin          = lnods(inode,ielem)
                        imate          = nodemat(ipoin)
                        icelltype      = int(celty_exm(1,ipoin))
                        ipair          = (imate-1)*EXM_CELLTYPE_MAXID + icelltype            
                        eland(:,inode) = land_variables_allmaterials(ipair,:)
                     end do

                     do igaus = 1,pgaus
                        do inode = 1,pnode
                           do i = 1,size(land_variables, KIND=ip)
                              j = (igaus - 1_ip)*nsdvs_ecc + i
                              state_ecc(j,ielem,:) = state_ecc(j,ielem,:) + eland(i,inode) * elmar(pelty)%shape(inode,igaus)
                           end do
                        end do !inode
                     end do !igauss
                  end if !pelty > 0 
               end do !ielem
         endif !has_land

      end if ! distribute parameters


      call memory_deallo(mem_modul(1:2,modul),'MATER_CELL_PAIRS','exm_iniunk',mater_cell_pairs)
      call memory_deallo(mem_modul(1:2,modul),'ONEOHR_STATUS'   ,'exm_iniunk',oneohr_status)
      call memory_deallo(mem_modul(1:2,modul),'ONETOR_STATUS' ,'exm_iniunk',onetor_status)
      call memory_deallo(mem_modul(1:2,modul),'LAND_VARIABLES',  'exm_iniunk',land_variables_allmaterials)




      call messages_live('  SETTING INITIAL CONDITIONS, NO RESTART...')
      !
      ! ELMAG
      !
      do ipoin = 1,npoin
         if (kfl_voini_exm(nodemat(ipoin)) == 1_ip) then
            elmag(ipoin,1) = voini_exm(nodemat(ipoin))    ! Hoy restart
         end if
      end do
      !
      ! FISOC
      !
      do ipoin = 1,npoin
         fisoc(ipoin) = -1.0_rp
         if( thiso_exm(1)>0.0_rp ) then !if ISHOCH AUTOSAVE
            fisoc(ipoin) = 0.0_rp !here we will save both upstrokes with positive time and downstrokes with negative time. There will never be activation at 0 time, so it's safe
         else
            fisoc(ipoin) = -1.0_rp
         end if
      end do
      !
      !  VAUXI_EXM and VCONC
      !
      do ipoin = 1,npoin
         kmodel_ipoin=kfl_cellmod(nodemat(ipoin))
         imate = nodemat(ipoin)
         if ((kmodel_ipoin == CELL_NOMOD_EXMEDI) .or. (kmodel_ipoin == CELL_FITZHUGH_EXMEDI)) then
            ! Al pedo, esta subru volo pq lo que hacia ya se hace aca
            !              call exm_inicia(kmodel_ipoin,ipoin)
         else if (kmodel_ipoin == CELL_TT2006_EXMEDI) then
            call exm_inihet(ipoin,imate)
            !call runend('EXM_INIUNK: JAZMIN, REVISAR ESTO DEL INIHET!!!')
         else if (kmodel_ipoin == CELL_OHARA_EXMEDI) then       !O'hara-Rudy
            call exm_iniohr(ipoin,imate)
         else if (kmodel_ipoin == CELL_TORORD_EXMEDI) then      ! ToR-ORd
            call exm_initor(ipoin,imate)
         else if (kmodel_ipoin == CELL_SCATRIA_EXMEDI) then     !Stem cell atria
            call exm_inisca(ipoin)
         else if (kmodel_ipoin == CELL_SCVENTRI_EXMEDI) then    !Stem cell ventricle
            call exm_iniscv(ipoin)
         end if
      end do
      !
      ! Put all components of ELMAG to initial values
      !
      call exm_updunk(ITASK_INIUNK)

   end if

end subroutine exm_iniunk
