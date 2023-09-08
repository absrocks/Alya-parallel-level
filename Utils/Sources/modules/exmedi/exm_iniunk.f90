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
  use      def_master
  use      def_domain
  use      mod_iofile
  use      def_exmedi
  use      mod_communications, only: PAR_INTERFACE_NODE_EXCHANGE, PAR_BROADCAST
  use mod_messages, only : livinf

  use mod_memory

  implicit none
  integer(ip) :: kmodel_imate, kmodel_ipoin, imate,ipoin,icomp,i,j,mat
  character(30) :: mesau(2)
  logical ::exm_oneohr


  ! auxiliary system matrix
  call memory_alloca(mem_modul(1:2,modul),'AMATR_AUXI_EXM','exm_memall',amatr_auxi_exm,size(amatr))     

  if (INOTMASTER) then
     !
     ! Initialise amatr, rhsid and friends
     !
     amatr= 0.0_rp
     rhsid= 0.0_rp     
     
     !
     ! Compute diagonal indices
     !
     if (solve_sol(1)%kfl_symme == 0) then
        call diagoi(npoin,1_ip,solve_sol(1)%kfl_symme,r_dom,c_dom,idima_exm)
     else
        call diagoi(npoin,1_ip,solve_sol(1)%kfl_symme,r_sym,c_sym,idima_exm)
     end if

  end if

  elmag    = 0.0_rp
  vcoin_exm= 0.0_rp
  vconc    = 0.0_rp
  vauxi_exm= 0.0_rp
  vicel_exm= 0.0_rp
  refhn_exm= 0.0_rp
  cell_ca0_ecc=0.0_rp

  !
  ! check if any TT-like model is present
  !
  ituss_exm= 1_ip
  kmodel_imate= 0_ip
  kmodel_ipoin= 0_ip
  do imate= 1,nmate
     kmodel_imate = kfl_cellmod(imate)

     if (INOTSLAVE) then
        !
        ! message: compute initial conditions for some cell models
        !
        if (kmodel_imate == CELL_TT2006_EXMEDI .or. kmodel_imate == CELL_OHARA_EXMEDI) then
           call livinf(42_ip,'  COMPUTING INITIAL PHYSIOLOGICAL CONDITIONS...',0_ip)
           !if (kmodel_imate == 4 .or. kmodel_imate==5 .or. kmodel_imate==6 ) then
           call livinf(42_ip,'    ITERATING ONE CELL MODEL...',0_ip)
           mesau(1)= intost(moneclmate_exm(1,imate)) ! beats
           mesau(2)= intost(moneclmate_exm(2,imate)) ! cyclelength
           call livinf(42_ip,'    BEATS=       '//trim(mesau(1)),0_ip)
           call livinf(42_ip,'    CYCLELENGTH= '//trim(mesau(2)),0_ip)
        end if
        !
        ! Initialization of the cell model: done only by the master and then distributed to the slaves
        !
        if (kmodel_imate == CELL_TT2006_EXMEDI) then
           call exm_onecel(imate)
        else if (kmodel_imate == CELL_OHARA_EXMEDI) then
           if( .NOT. exm_oneohr(imate) ) then
              call runend('EXM_INIUNK: Steady state of the cell model(s) not reached. Try increasing the number of beats in the TThet_framework')
           endif
        end if
        !
        ! message: done
        !
        call livinf(42_ip,'  INITIAL CONDITIONS COMPUTED MATERIAL= '//trim(intost(imate)),0_ip)
        if (kmodel_imate == CELL_TT2006_EXMEDI) then
           call livinf(42_ip,'    SUB CELL MODEL: TT HETEROGENEOUS',0_ip)
        else if (kmodel_imate == CELL_OHARA_EXMEDI) then
           call livinf(42_ip,'    SUB CELL MODEL: OHARA - RUDY',0_ip)
        end if

     end if
  end do
  !
  ! Master distributes cell initialized values to the slaves
  ! 
  do imate=1,nmate
     do i=1,3_ip
        call PAR_BROADCAST(vminimate_exm(i,imate))
        do j=1,nconc_exm
           call PAR_BROADCAST(vcoin_exm(j,i,imate))
        end do
        do j=1,nauxi_exm
           call PAR_BROADCAST(vauin_exm(j,i,imate))
        end do
     end do
  end do


  !
  ! Master distributes variables required for Land
  !
  if( kfl_coupl(ID_SOLIDZ,ID_EXMEDI) >= 1_ip .or. kfl_coupl(ID_EXMEDI,ID_SOLIDZ) >=1_ip ) then
     do i=1,3_ip
        do j=1,30_ip
           call PAR_BROADCAST(vaulo_exm(j,i))
        enddo
        do j=1,12_ip
           call PAR_BROADCAST(vcolo_exm(j,i))
        enddo
        do j=1,27_ip
           call PAR_BROADCAST(viclo_exm(j,i))
        enddo

        call PAR_BROADCAST(elmlo_exm(i))
     enddo

     do i=1,6_ip
        call PAR_BROADCAST(stateland(i,1,1,1))
        stateland(i,:,:,:)=stateland(i,1,1,1)
     enddo

     call PAR_BROADCAST(troponin(1))
     troponin(:)=troponin(1)
     call PAR_BROADCAST(troponin_prev(1))
     troponin_prev(:)=troponin(1)

  endif


  ! Save Ca0 values to automatically calibrate Hunter-McCulloch ECC model in sld_eccoup
  do imate= 1,nmate
     cell_ca0_ecc(:,imate)=vcoin_exm(1,:,imate)
  enddo


  if (kfl_rstar == 0_ip) then

     if (INOTSLAVE) then
        call livinf(42_ip,'  SETTING INITIAL CONDITIONS, NO RESTART...',0_ip)
     end if
     
     if(INOTMASTER) then           
        if (modac_exm < 0 .and. associated(xfiel)) then
           do icomp=1,ncomp_exm
              do ipoin=1,npoin
                 appfi_exm(ipoin) = xfiel(-modac_exm) % a(1,ipoin,1)
              end do
           end do
        else
        end if

        do ipoin=1,npoin
          if (kfl_voini_exm(nodemat(ipoin)) == 1_ip) then
            elmag(ipoin,1:ncomp_exm) = voini_exm(nodemat(ipoin))    ! initial monodomain voltage
          end if
        end do


        do ipoin=1,npoin
           kmodel_ipoin=kfl_cellmod(nodemat(ipoin)) 
           mat = nodemat(ipoin)                
           if ((kmodel_ipoin == CELL_NOMOD_EXMEDI) .or. (kmodel_ipoin == CELL_FITZHUGH_EXMEDI)) then
              ! Al pedo, esta subru volo pq lo que hacia ya se hace aca
              !              call exm_inicia(kmodel_ipoin,ipoin)
           else if (kmodel_ipoin == CELL_TT2006_EXMEDI) then
              call exm_inihet(kmodel_ipoin,ipoin,mat)
              !call runend('EXM_INIUNK: JAZMIN, REVISAR ESTO DEL INIHET!!!')
           else if (kmodel_ipoin == CELL_OHARA_EXMEDI) then    !O'hara-Rudy              
              call exm_iniohr(kmodel_ipoin,ipoin,mat)
           else if (kmodel_ipoin == CELL_SCATRIA_EXMEDI) then    !Stem cell atria             
              call exm_inisca(kmodel_ipoin,ipoin)
           else if (kmodel_ipoin == CELL_SCVENTRI_EXMEDI) then    !Stem cell ventricle           
              call exm_iniscv(kmodel_ipoin,ipoin)
           end if
        end do
     end if

  else
     ! read restart  variables from exmedi rst files
     call exm_restar(1_ip)
     
  end if


end subroutine exm_iniunk
