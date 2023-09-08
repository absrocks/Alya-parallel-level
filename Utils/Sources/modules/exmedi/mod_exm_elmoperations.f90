!-----------------------------------------------------------------------
!> @addtogroup ExmediElmoperations
!> @ingroup    Exmedi
!> @{
!> @file    mod_exm_elmoperations.f90
!> @author  Mariano Vazquez
!> @date    22/11/2016
!> @brief   Compute elemental matrix and RHS  
!> @details Compute elemental matrix and RHS
!-----------------------------------------------------------------------
module mod_exm_elmoperations
  use      def_master
  use      def_domain
  use      def_elmtyp
  use      def_exmedi
  use      mod_messages, only : livinf
  use      mod_exm_fractional_diffusion
  use mod_timings,                     only : timings_assembly

  implicit none

  private

  integer(ip), parameter :: VECTOR = 1_ip
  logical  :: flag_pseudecg
  real(rp) :: time1,time2

  !
  ! Public stuff
  !
  public :: exm_elmoperations


  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  
contains

  !-----------------------------------------------------------------------
  !> @author  Mariano Vazquez
  !> @date    22/11/2016
  !> @brief   Compute elemental matrix and RHS  
  !> @details Compute elemental matrix and RHS
  !-----------------------------------------------------------------------
  subroutine exm_elmoperations
    implicit none

    integer(ip) :: remainder,ielem_start_block,last_block,ipoin
    logical(lg), save  :: flag_compute_matrix=.true.
    integer(ip), parameter :: EXPLICIT_SCHEME=1, IMPLICIT_SCHEME=2, CRANK=2

    time1     =  0.0_rp
    time2     =  0.0_rp

    if (INOTSLAVE) then
       if (kfl_timet_exm == EXPLICIT_SCHEME) call livinf(165_ip,'  SOLVE TISSUE EXPLICITLY... ',0_ip)
       if (kfl_timet_exm == IMPLICIT_SCHEME) call livinf(165_ip,'  SOLVE TISSUE IMPLICITLY... ',0_ip)
    end if

    flag_pseudecg = .false.
    if (kcopeecg_exm .ge. 0) then
       !
       ! Compute pseudo-ECG
       !
       if (mod(ittim,nvint_exm) == 0) then
          kcopeecg_exm = kcopeecg_exm + 1
          flag_pseudecg= .true.
          pseudecg_exm   = 0.0_rp
       end if
    end if


    call cputim(time1)

    if ((kfl_comat_exm == 1 .or. ( (ittim == 1) .and. itinn(modul)==1)) .or. (flag_pseudecg) .or. (flag_compute_matrix))  then

       call livinf(165_ip,'  COMPUTE MATRIX... ',0_ip)
       flag_compute_matrix=.false.


       !
       ! Loop over blocks of elements
       !
       rhsid     = 0.0_rp
       amatr     = 0.0_rp

       remainder = mod(nelem,VECTOR)
       last_block = nelem - remainder
       elements: do ielem_start_block = 1,last_block,VECTOR
          call exm_elmope_processblock(ielem_start_block,VECTOR)
       end do elements
       ielem_start_block = last_block + 1
       if (ielem_start_block .le. nelem) call exm_elmope_processblock(ielem_start_block,remainder)
       !
       ! Store amatr as it was computed
       !
       if (kfl_timet_exm == IMPLICIT_SCHEME) amatr_auxi_exm = amatr       

       call cputim(time2)

    end if

    !
    ! CPU times
    !
    cpu_ass_sol_exm(1) = time2 - time1

       call timings_assembly(cpu_ass_sol_exm(1), TYPE_OF_ASSEMBLY='ELEMENT')

  end subroutine exm_elmoperations

  !-----------------------------------------------------------------------
  !> @author  Mariano Vazquez
  !> @date    22/11/2016
  !> @brief   Compute elemental matrix and RHS  
  !> @details Compute elemental matrix and RHS
  !-----------------------------------------------------------------------
  subroutine exm_elmope_processblock(ielem_start_block,block_size)
    use      def_master
    use      def_domain
    use      def_elmtyp
    use      def_exmedi
    implicit none

    integer(ip) :: ielem_start_block !> ielem_start_block
    integer(ip) :: block_size        !> block_size is VECTOR or remainder = mod(nelem,VECTOR)

    integer(ip) :: ivelem,kelem,kelty,mnode_block,mgaus_block

    integer(ip), parameter :: EXPLICIT_SCHEME=1, IMPLICIT_SCHEME=2, CRANK=2

    integer(ip) :: &
         ielem(VECTOR), &
         pelty(VECTOR), &
         pnode(VECTOR), &
         pevat(VECTOR), &
         pgaus(VECTOR), &
         plapl(VECTOR), &
         pmate(VECTOR), &
         noion(VECTOR), &
         vlnods(VECTOR,mnode)

    integer(ip) :: &
         inode,jnode,idime,jdime,kdime,idim2,ipoin,igaus,iroot


    real(rp)    :: &
         elrhs(VECTOR,mnode),&
         elmat(VECTOR,mnode,mnode)

    real(rp)    :: &
         eliap_iterk(VECTOR, mnode ),&
         elfib(VECTOR,ndime,mnode  ),&
         elcod(VECTOR,ndime,mnode  ),&
         eldistecg(VECTOR,mnode,nrootecg_exm),&
         elfin(VECTOR,ndime,ndime,mnode  ),&
         elfde(VECTOR            ,mnode  ),&
         eldis(VECTOR,ndime,mnode  )


    real(rp)    :: gpcar(VECTOR,ndime,mnode,mgaus),griap(VECTOR,ndime),gpsha_inode,gpsha_jnode,dvgaus
    real(rp)    :: dvolu(VECTOR,mgaus),gpdet(VECTOR,mgaus)                     ! Values at Gauss points
    real(rp)    :: hessi(VECTOR,ntens,mnode*mlapl)
    real(rp)    :: xjaci(VECTOR,ndime,ndime,mgaus),xjacm(VECTOR,ndime,ndime,mgaus)
    real(rp)    :: difin(VECTOR,ndime,ndime,mgaus), nfibe(VECTOR,ndime,mgaus),nfibe_length(VECTOR)
    real(rp)    :: difinaux(ndime,ndime,mgaus), gridi(ndime)
    real(rp)    :: gpfin(VECTOR,ndime,ndime)      , gpfde(VECTOR)

    real(rp)    :: d2sdx(ndime,ndime,ndime*mlapl)        ! that is called inside do igaus 

    real(rp)    :: crank_factor, dtimeEP, distroo, dauxi

    ! Initialisation

    eldistecg = 0.0_rp
    difinaux = 0.0_rp
    elmat = 0.0_rp
    elrhs = 0.0_rp
    griap = 0.0_rp
    pelty = 0_ip
    pnode = 0_ip
    pevat = 0_ip
    pgaus = 0_ip
    plapl = 0_ip
    pelty = 0_ip
    pmate = 0_ip
    vlnods = 0_ip
    crank_factor = 1.0_rp
    if (kfl_tisch_exm == CRANK) crank_factor = 0.5_rp

    dtimeEP = dtime

    set_up_loop: do ivelem= 1,block_size
       kelem = ivelem + ielem_start_block - 1
       kelty = ltype(kelem)
       ielem(ivelem) = kelem
       pelty(ivelem) = kelty
       pnode(ivelem) = nnode(kelty)
       pgaus(ivelem) = ngaus(kelty)
       plapl(ivelem) = llapl(kelty)
       pmate(ivelem) = 1_ip
       pevat(ivelem) = pnode(ivelem)
       if( nmate > 1_ip ) then
          pmate(ivelem) = lmate_exm(kelem)
       end if

       do inode=1,pnode(ivelem)
          ipoin= lnods(inode,kelem)
          vlnods(ivelem,inode) = ipoin

          ! gather:
          eliap_iterk(ivelem,inode) = elmag(ipoin,1)   
          elcod(ivelem,1:ndime,inode)=coord(1:ndime,ipoin)
          if( kfl_coupl(ID_SOLIDZ,ID_EXMEDI) >= 1 .or. kfl_coupl(ID_EXMEDI,ID_SOLIDZ) >=1 ) then
             if (kfl_gcoup_exm == 1) then
                eldis(ivelem,1:ndime,inode)= displ(1:ndime,ipoin,1)
                elcod(ivelem,1:ndime,inode)= elcod(ivelem,1:ndime,inode) + eldis(ivelem,1:ndime,inode)
                if (modfi_exm < 0) then
                   elfib(ivelem,1:ndime,inode)= fiber(1:ndime,ipoin)
                else if ( modfi_exm > 0 ) then
                   gfibe_exm= 0.0_rp
                   gfibe_exm(modfi_exm) = 1.0_rp
                   elfib(ivelem,1:ndime,inode)= gfibe_exm(1:ndime)
                end if
             else if (kfl_gcoup_exm == 2) then
                elfin(ivelem,1:ndime,1:ndime,inode)= gdeinv(1:ndime,1:ndime,ipoin)
                elfde(ivelem,                inode)= gdedet(                ipoin)
             end if
          end if

          if (flag_pseudecg) then
             do iroot= 1,nrootecg_exm
                distroo = 0.0_rp
                do idime= 1,ndime
                   distroo = distroo + &
                        (elcod(ivelem,idime,inode)-coordecg_exm(idime,iroot)) &
                        * (elcod(ivelem,idime,inode)-coordecg_exm(idime,iroot)) 
                end do
                if (distroo > 0) eldistecg(ivelem,inode,iroot) = 1.0_rp / sqrt(distroo)
             end do
          end if

       end do

    end do set_up_loop

    ! Block-wide maxima
    mnode_block = maxval(pnode)
    mgaus_block = maxval(pgaus)

    compute_loop: do ivelem= 1,block_size

       ! Cartesian derivatives and Jacobian 
       gauss_points_cartesian: do igaus=1,pgaus(ivelem)
          call elmder(pnode(ivelem),ndime,elmar(pelty(ivelem))%deriv(1,1,igaus),elcod(ivelem,1,1),&
               gpcar(ivelem,1,1,igaus),gpdet(ivelem,igaus),xjacm(ivelem,1,1,igaus),xjaci(ivelem,1,1,igaus))
          dvolu(ivelem,igaus)=elmar(pelty(ivelem))%weigp(igaus)*gpdet(ivelem,igaus)
          if(plapl(ivelem)==1) & 
               call elmhes(elmar(pelty(ivelem))%heslo(1,1,igaus),hessi(ivelem,1,1),ndime,pnode(ivelem),ntens, &
               xjaci(ivelem,1,1,igaus),d2sdx,elmar(pelty(ivelem))%deriv,elcod(ivelem,1,1))
       end do gauss_points_cartesian

       nfibe_length(ivelem)=0.0_rp

       ! Compute extra and intracellular conductivity matrices        
       if ( kfl_coupl(ID_SOLIDZ,ID_EXMEDI) >= 1 .or. kfl_coupl(ID_EXMEDI,ID_SOLIDZ) >=1) then          
          if (kfl_gcoup_exm == 1) then
             if (ittim == 1 .or. kfl_cellmod(pmate(ivelem)) == CELL_NOMOD_EXMEDI) then 
                ! either the first time gpfib is still zero or this is a NO_MOD exmedi material
                call exm_comcnd(ielem(ivelem),pmate(ivelem),difin(ivelem,1,1,1),noion(ivelem))
             else if (kfl_gcoup_exm == 1) then
                ! only ==1 because in the new lagrangian formulation (==2), difin is fine as it is in the reference system
                do igaus=1,pgaus(ivelem)
                   do idime=1,ndime
                      nfibe(ivelem,idime,igaus)= gpfib(idime,igaus,ielem(ivelem))
                      nfibe_length(ivelem)= nfibe_length(ivelem) + nfibe(ivelem,idime,igaus) * nfibe(ivelem,idime,igaus)
                   end do
                   nfibe_length(ivelem)= sqrt(nfibe_length(ivelem))
                   do idime=1,ndime
                      nfibe(ivelem,idime,igaus)= nfibe(ivelem,idime,igaus) / nfibe_length(ivelem)
                   end do
                   call exm_coucnd(pmate(ivelem),difin(ivelem,1,1,igaus),nfibe(ivelem,1,igaus),noion(ivelem))
                end do
             end if
          else
             ! do this when both gcoup_exm == 0 or == 2
             call exm_comcnd(ielem(ivelem),pmate(ivelem),difin(ivelem,1,1,1),noion(ivelem))
          end if
       else
          call exm_comcnd(ielem(ivelem),pmate(ivelem),difin(ivelem,1,1,1),noion(ivelem))
       end if

       !
       ! As this is a iteration increment formulation, FFD, BFD and CN needs to evaluate the RHS
       !
       gauss_points_rhs: do igaus=1,pgaus(ivelem)

          dvgaus= dvolu(ivelem,igaus)

          ! Intracellular action potential gradients 
          do inode=1,pnode(ivelem)
             griap(ivelem,1:ndime)=griap(ivelem,1:ndime) &
                  + gpcar(ivelem,1:ndime,inode,igaus) * eliap_iterk(ivelem,inode)
             gpsha_inode  = elmar(pelty(ivelem))%shape(inode,igaus)           
          end do

          !
          ! Compute diffusion
          ! 
          if ( kfl_coupl(ID_SOLIDZ,ID_EXMEDI) >= 1 .or. kfl_coupl(ID_EXMEDI,ID_SOLIDZ) >=1) then          
             if (kfl_gcoup_exm == 2) then
                ! deformation gradient inverse and determinant
                gpfin(ivelem,1:ndime,1:ndime) = 0.0_rp
                gpfde(ivelem                ) = 0.0_rp
                do inode=1,pnode(ivelem)
                   gpsha_inode  = elmar(pelty(ivelem))%shape(inode,igaus)           
                   gpfin(ivelem,1:ndime,1:ndime)=gpfin(ivelem,1:ndime,1:ndime) &
                        + gpsha_inode * elfin(ivelem,1:ndime,1:ndime,inode)
                   gpfde(ivelem)=gpfde(ivelem) &
                        + gpsha_inode * elfde(ivelem,inode)                        
                end do

                ! D^aux_ij = J Finv_ik D_kl Finv_jl (la segunda F es asi porque es la traspuesta en realidad) 
                do idime= 1,ndime
                   do jdime= 1,ndime
                      do kdime= 1,ndime
                         do idim2= 1,ndime
                            difinaux(idime,jdime,igaus)= difinaux(idime,jdime,igaus) &
                                 + gpfde(ivelem) &
                                 * gpfin(ivelem,idime,kdime) & 
                                 * difin(ivelem,kdime,idim2,igaus) &
                                 * gpfin(ivelem,jdime,idim2)
                         end do
                      end do
                   end do
                end do
             else 

                difinaux(1:ndime,1:ndime,igaus)= difin(ivelem,1:ndime,1:ndime,igaus)

             end if

          else

             difinaux(1:ndime,1:ndime,igaus)= difin(ivelem,1:ndime,1:ndime,igaus)

          end if

          !
          ! Compute pseudo ecg
          !
          if (flag_pseudecg) then
             dvgaus= dvolu(ivelem,igaus)
             do iroot= 1,nrootecg_exm
                ! grad (1 / R) 
                gridi = 0.0_rp
                do inode=1,pnode(ivelem)
                   do idime= 1,ndime
                      gridi(idime)=gridi(idime) &
                           + gpcar(ivelem,idime,inode,igaus) * eldistecg(ivelem,inode,iroot)
                   end do
                end do
                dauxi = 0.0_rp
                ! Compute (grad(V)))_i*D_ij*grad(1/R)_j
                do idime= 1,ndime
                   do jdime= 1,ndime
                      dauxi = dauxi + griap(ivelem,idime) * difin(ivelem,idime,jdime,igaus) * gridi(jdime) 
                   end do
                end do
                pseudecg_exm(iroot)= pseudecg_exm(iroot) + dauxi * dvgaus
             end do
          end if

       end do gauss_points_rhs

       if (kfl_timet_exm == IMPLICIT_SCHEME) then
          !
          ! 
          !
          gauss_points_matrix: do igaus=1,pgaus(ivelem)

             dvgaus= dvolu(ivelem,igaus)

             do inode= 1,pnode(ivelem)
                gpsha_inode  = elmar(pelty(ivelem))%shape(inode,igaus)           
                do jnode= 1,pnode(ivelem)
                   do idime= 1,ndime
                      do jdime= 1,ndime
                         elmat(ivelem,inode,jnode)= elmat(ivelem,inode,jnode) &
                              + dvgaus * gpcar(ivelem,idime,inode,igaus) & 
                              * difinaux(idime,jdime,igaus) * gpcar(ivelem,jdime,jnode,igaus) * crank_factor
                      end do
                   end do
                end do
             end do
          end do gauss_points_matrix

       end if

       !
       ! Matrix assembly
       !
       if (kfl_timet_exm == IMPLICIT_SCHEME) call assmat(&
            solve(1)%ndofn,pnode(ivelem),pevat(ivelem),solve(1)%nunkn,&
            solve(1)%kfl_algso,ielem(ivelem),lnods(1,ielem(ivelem)),elmat(ivelem,1,1),amatr)

!!$       !
!!$       ! RHS assembly
!!$       !
!!$       if (kfl_comat_exm == 1) call assrhs(&
!!$            1_ip,pnode(ivelem),lnods(1,ielem(ivelem)),elrhs(ivelem,1),rhsid)


    end do compute_loop



  end subroutine exm_elmope_processblock

end module mod_exm_elmoperations
!> @}
