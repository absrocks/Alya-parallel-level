!------------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_assbou.f90
!> @date    21/06/2012
!> @author  Mariano Vazquez
!> @brief   Boundary conditions assembly for implicit scheme
!> @details Boundary conditions assembly for implicit schemes \n
!!          Conditions are always set on physical variables.
!> @}
!------------------------------------------------------------------------
subroutine nsa_assbou(pdofn,pnode,pevat,kfl_elfix,lnode,&
     elrhs,elmat,elunk,elphy,elhcp,elwme,elbve,elsou)
  !-----------------------------------------------------------------------
  !
  ! BOUNDARY CONDITIONS:
  !
  !  
  !-----------------------------------------------------------------------
  use      def_master
  use      def_nastal, only: cvcoe_nsa,ncomp_nsa,ndofn_nsa,kfl_fixrs_nsa,jacrot_du_dq_nsa,&
       jacrot_dq_du_nsa,runiv_nsa,kfl_linea_nsa
  use      def_domain, only: ndime,mnode,skcos,exnor,lpoty,lnods
  implicit none

  integer(ip)    :: pdofn,pnode,pevat,inode,jnode,ipoin,ibopo,idime,jdime,ievat,jevat
  integer(ip)    :: kfl_elfix(pdofn,pnode),lnode(pnode),kinfl,iroty,idofn,jdofn,nofix,nofix_vel
  real(rp)       :: elrhs(pevat),elmat(pevat,pevat), elsou(pevat),adiag, xvalu, xmach,&
       elunk(ndime+2,mnode,ncomp_nsa),elbve(ndime+2,mnode),velsq,xvelo(ndime),xadve(ndime),&
       xpres,xdens,xtemp,xener,xmome(ndime), elphy(ndofn_nsa,mnode),&
       rotgl(ndime,ndime),rotlg(ndime,ndime),rotqu_aux(ndofn_nsa,ndofn_nsa),&
       rotuq_aux(ndofn_nsa,ndofn_nsa),&
       elhcp(mnode),elwme(mnode),&
       rgacv,rgasc,xhecv,matri_debu(ndofn_nsa,ndofn_nsa),matri_debu_ndime(ndime,ndime)
  !
  ! Boundary conditions: Dirichlet (essential B.Cs)
  !
  ! PROGRAMMED BOUNDARY DIRICHLET BOUNDARY CONDITIOS:
  !
  ! 00000 : free, no boundary, (U,rho,E)
  ! 11111 : all prescribed, (U,rho,E)
  ! 50000 : checking inflow/outflow 
  ! 20000 : velocity local base prescribed, (u,rho,T)
  ! 20100 : velocity local base prescribed, (u,rho,t) "the famous condition 20100"
  ! 11101 : no slip NS, (u,rho,T) or subsonic inflow 
  ! 11103 : no slip NS, (u,rho,T) 
  ! 11100 : no slip NS, (u,rho,T) but no prescription on T 
  ! 00010 : subsonic outflow
  ! 11110 : dani 1
  ! 11102 : subsonic inflow, prescribing p instead of T (u,rho,p)
  ! 00002 : subsonic outflow, prescribing p instead of rho (u,rho,p)
  !

  do inode = 1,pnode
     ipoin = lnode(inode)
     ibopo = lpoty(ipoin)

     !
     ! Check if this node has some kind of boundary condition
     !
     nofix = 0
     do idofn= 1,pdofn
        if (kfl_elfix(idofn,inode) .gt. 0) then
           nofix= nofix + 1
        end if
     end do
     nofix_vel= 0
     do idime= 1,ndime
        if (kfl_elfix(idime,inode) .gt. 0) then
           nofix_vel= nofix_vel + 1
        end if
     end do

     if (nofix > 0 .and. ibopo > 0) then 

        ! this is a boundary node with a condition

        rgasc = runiv_nsa / elwme(inode)
        xhecv = elhcp(inode) - rgasc                      !Cv is computed from R & Cp
        rgacv = runiv_nsa / elwme(inode) / xhecv          !R / Cv

        xdens= elunk(ndime+1,inode,ITER_K)
        xener= elunk(ndime+2,inode,ITER_K)
        xpres= elphy(ndime+1,inode)
        xtemp= elphy(ndime+2,inode)
        velsq= 0.0_rp
        do idime= 1,ndime
           xvelo(idime) = elphy(idime,inode)
           xadve(idime) = elphy(idime,inode)
           xmome(idime) = elunk(idime,inode,ITER_K)
           velsq = velsq + xvelo(idime)*xvelo(idime)
        end do
        ! When coupled with alefor, substract the mesh velocity velom to the advection velocity xadve
        if( kfl_coupl(ID_NASTAL,ID_ALEFOR) /= 0 ) then  
           do idime=1,ndime
              xadve(idime) = xadve(idime) - velom(idime,ipoin)
           end do
        end if


        if(kfl_elfix(1,inode)==5) then         
           !
           ! 1. Check wether inflow or outflow, sub or supersonic
           !

           if( solve(1) % kfl_iffix .ne. 0 ) call runend("NSA_ASSBOU: BC 50000 NOT PREPARED FOR ZEROFIX, TO BE PROGRAMMED.")

           call nsa_chkinf(&
                kinfl,xmach,ipoin,xadve(1:ndime),xvelo(1:ndime),xpres,xdens)

           if (kinfl==1) then
              ! inflow
              ! velocity prescribed
              do idime=1,ndime
                 kfl_elfix(idime,inode)= 1  
              end do
              kfl_elfix(ndime+2,inode)= 1
              if (xmach >= 1.0_rp) then 
                 ! supersonic, so fix also the density
                 kfl_elfix(ndime+1,inode)  = 1
              end if
           else
              ! outflow
              if (xmach < 1.0_rp) then 
                 ! subsonic, so fix only the density
                 kfl_elfix(ndime+1,inode)  = 1
              end if
           end if

        end if

        !
        ! Change variables from conservative to physical 
        !
        ! 1. Compute base change jacobians dU/dQ and dQ/dU
        !
        ! These are the boundary conditions options so far:
        !
        ! (a) u,rho,T (so we use U,rho,E as unknowns) : no jacobians computed -> nofix_vel  > 0 and nofix == ndofn_nsa
        ! (b) u, T    : jacobians computed    -> nofix_vel  > 0 and nofix  < ndofn_nsa
        ! (c) rho     : no jacobians computed -> nofix_vel == 0 and nofix  < ndofn_nsa
        ! (d) u       : jacobians computed, but only for the velocity local frame -> nofix_vel > 0 and nofix < ndofn_nsa,
        !                                                                            but kfl_elfix(ndime+1,inode) == 0
        ! (e) u, rho  : no jacobians computed    -> nofix_vel  > 0 and nofix  < ndofn_nsa
        ! (f) u, p    : jacobians computed       -> nofix_vel  > 0 and nofix  < ndofn_nsa
        ! (g) p       : jacobians computed       -> nofix_vel == 0 and nofix  < ndofn_nsa

        if ((nofix < ndofn_nsa .and. kfl_elfix(ndime+1,inode) == 0) .or. kfl_elfix(ndime+2,inode) == 2) then 
           ! (b), (d) and (f)
           ! Initialize rotation matrices. As this is done over ibopos but within an element loop,
           ! the rotation matrices are computed and stored more than once, redundantly. 
           ! Pero bueno, la vida es asín. 
           do idofn= 1,ndofn_nsa
              do jdofn= 1,ndofn_nsa
                 jacrot_dq_du_nsa(idofn,jdofn,ibopo)=0.0_rp
                 jacrot_du_dq_nsa(idofn,jdofn,ibopo)=0.0_rp
              end do
              jacrot_dq_du_nsa(idofn,idofn,ibopo)=1.0_rp
              jacrot_du_dq_nsa(idofn,idofn,ibopo)=1.0_rp              
           end do

           ! Compute rotqu and rotuq 

           if (kfl_elfix(ndime+2,inode) == 2) then   ! U=(U, rho, E) and Q=(u, rho, p)

              do idime=1,ndime
                 jacrot_du_dq_nsa(  idime,  idime,ibopo) =   xdens         ! dU_i / du_i 
                 jacrot_du_dq_nsa(  idime,ndime+1,ibopo) =   xvelo(idime)  ! dU_i / drho
                 jacrot_du_dq_nsa(ndime+2,  idime,ibopo) =   xmome(idime)  ! dE / du_i   
                 
                 jacrot_dq_du_nsa(  idime,  idime,ibopo) =   1.0_rp / xdens            ! du / dU
                 jacrot_dq_du_nsa(  idime,ndime+1,ibopo) = - xvelo(idime) / xdens      ! du / drho           
                 jacrot_dq_du_nsa(ndime+2,  idime,ibopo) = - xvelo(idime) * rgacv      ! dp / dU_i                 
              end do

              jacrot_du_dq_nsa(ndime+2,ndime+1,ibopo) = 0.5_rp * velsq                 ! dE / drho
              jacrot_du_dq_nsa(ndime+2,ndime+2,ibopo) = 1.0 / rgacv                    ! dE / dp
              jacrot_du_dq_nsa(ndime+1,ndime+1,ibopo) = 1.0_rp                         ! drho / drho
              
              jacrot_dq_du_nsa(ndime+2,ndime+1,ibopo) = rgacv * 0.5_rp * velsq         ! dp / drho
              jacrot_dq_du_nsa(ndime+2,ndime+2,ibopo) = rgacv                          ! dp / dE
              jacrot_dq_du_nsa(ndime+1,ndime+1,ibopo) = 1.0_rp                         ! drho / drho 
              
           else if (kfl_elfix(ndime+1,inode) == 0) then    ! U=(U, rho, E) and Q=(u, rho, T)
              ! u and T prescribed
              do idime=1,ndime
                 jacrot_du_dq_nsa(  idime,  idime,ibopo) =   xdens         ! dU_i / du_i = d (rho u_i) / du_i 
                 jacrot_du_dq_nsa(  idime,ndime+1,ibopo) =   xvelo(idime)  ! dU_i / drho = d (rho u_i) / drho 
                 jacrot_du_dq_nsa(ndime+2,  idime,ibopo) =   xmome(idime)  ! dE / du_i   = d (rho c_v T + rho 0.5 u^2) / du_i     
                 
                 jacrot_dq_du_nsa(  idime,  idime,ibopo) =   1.0_rp / xdens
                 jacrot_dq_du_nsa(  idime,ndime+1,ibopo) = - xvelo(idime)/ xdens                   
                 jacrot_dq_du_nsa(ndime+2,  idime,ibopo) = - xvelo(idime) / xhecv / xdens 
                 
              end do
              jacrot_du_dq_nsa(ndime+2,ndime+1,ibopo) = xtemp * xhecv + 0.5_rp * velsq    ! dE / drho = d (rho c_v T + rho 0.5 u^2) / drho
              jacrot_du_dq_nsa(ndime+2,ndime+2,ibopo) = xdens * xhecv                     ! dE / dT   = d (rho c_v T + rho 0.5 u^2) / dT
              jacrot_du_dq_nsa(ndime+1,ndime+1,ibopo) = 1.0_rp                            ! drho / drho 
              
              jacrot_dq_du_nsa(ndime+2,ndime+1,ibopo) = (velsq - xener / xdens) / xhecv / xdens
              jacrot_dq_du_nsa(ndime+2,ndime+2,ibopo) = 1.0_rp / xhecv / xdens 
              jacrot_dq_du_nsa(ndime+1,ndime+1,ibopo) = 1.0_rp                            ! drho / drho 

           end if

           !           if( kfl_elfix(ndime+2,inode) > 0 ) then                       ! (b), temperature prescribed             
           !              do idime=1,ndime
           !                 jacrot_du_dq_nsa(ndime+2,  idime,ibopo) =    xmome(idime) / xdens
           !                 jacrot_dq_du_nsa(  idime,ndime+2,ibopo) =   -xmome(idime) / xhecv / xdens / xdens
           !              end do
           !              jacrot_du_dq_nsa(ndime+2,ndime+1,ibopo) = xtemp * xhecv - 0.5_rp * velsq
           !              jacrot_du_dq_nsa(ndime+2,ndime+2,ibopo) = xdens * xhecv
           !              jacrot_dq_du_nsa(ndime+1,ndime+2,ibopo) = (velsq - xener / xdens) / xhecv / xdens 
           !              jacrot_dq_du_nsa(ndime+2,ndime+2,ibopo) = xener /  xhecv / xdens 
           !           end if


           !
           ! 2. Check the specific prescription
           !
           if( kfl_elfix(1,inode) == 2 ) then                             ! compute local frame
              !
              ! 1. First, check local base velocity prescriptions
              !
              !   rotau(1:ndime,1:ndime,1)  is  GL , goes in rotuq
              !   rotau(1:ndime,1:ndime,2)  is  LG , goes in rotqu
              !
              !   Choose the proper local basis
              !
              iroty=kfl_fixrs_nsa(ibopo)
              if( iroty == -1 ) then                                    ! Tangent system
                 do idime=1,ndime
                    do jdime= 1,ndime
                       rotlg(idime,jdime)= exnor(idime,jdime,ibopo)
                       rotgl(idime,jdime)= exnor(jdime,idime,ibopo)
                    end do
                 end do

              else if( iroty >= 1 ) then                                ! Given system
                 do idime=1,ndime
                    do jdime= 1,ndime
                       rotlg(idime,jdime)= skcos(idime,jdime,iroty)
                       rotgl(idime,jdime)= skcos(jdime,idime,iroty)
                    end do
                 end do
                 !           else if( iroty == -2 ) then                               ! Given system
              else if( iroty == -3 ) then                               ! Geometrical normal
                 do idime=1,ndime
                    do jdime= 1,ndime
                       rotlg(idime,jdime)= skcos(idime,jdime,ibopo)
                       rotgl(idime,jdime)= skcos(jdime,idime,ibopo)
                    end do
                 end do
              end if


              !
              ! 2. Correct jacrot_du_dq_nsa and jacrot_dq_du_nsa accordingly: 
              !    F   = dU/dQ *   LG   in jacrot_du_dq_nsa(...)  
              !    F^⁻1=  GL  * dQ/dU   in jacrot_dq_du_nsa(...)

              rotqu_aux= 0.0_rp
              rotuq_aux= 0.0_rp
              do idime=1,ndime
                 do idofn= 1,ndofn_nsa
                    do jdime=1,ndime
                       rotqu_aux(idofn,idime)= &
                            rotqu_aux(idofn,idime) + jacrot_du_dq_nsa(idofn,jdime,ibopo) * rotlg(jdime,idime) 
                       rotuq_aux(idime,idofn)= &
                            rotuq_aux(idime,idofn) + rotgl(idime,jdime) * jacrot_dq_du_nsa(jdime,idofn,ibopo)
                    end do
                 end do
              end do

              do idime=1,ndime
                 do idofn= 1,ndofn_nsa
                    jacrot_du_dq_nsa(idofn,idime,ibopo)=rotqu_aux(idofn,idime) 
                    jacrot_dq_du_nsa(idime,idofn,ibopo)=rotuq_aux(idime,idofn) 
                 end do
              end do

           end if

           !
           ! Rotate elmat and elrhs
           !
           call nsa_rotsys(1_ip,&
                inode,pnode,ndofn_nsa,pevat,elmat,elrhs,&
                jacrot_du_dq_nsa(1,1,ibopo),jacrot_dq_du_nsa(1,1,ibopo),kfl_linea_nsa)

           !           if( kfl_elfix(1,inode) == 2) then
           !              matri_debu=0.0_rp
           !              matri_debu_ndime=0.0_rp
           !              matri_debu=matmul(jacrot_du_dq_nsa(:,:,ibopo),jacrot_dq_du_nsa(:,:,ibopo))
           !              matri_debu_ndime=matmul(rotgl(:,:),rotlg(:,:))
           !              matri_debu_ndime=0.0_rp
           !           end if

        end if


        if( solve(1) % kfl_iffix == 0 ) then


           !
           ! Correct boundary conditions
           !
           !
           !matrix and rhsid modified now to account for the boundary conditions
           !

           ! Momentum dof's
           do idime = 1,ndime
              if( kfl_elfix(idime,inode) == 1 .or. kfl_elfix(idime,inode) == 2) then
                 ! either local frame or cartesian frame, both are equally treated
                 ievat = (inode-1)*pdofn + idime
                 adiag = elmat(ievat,ievat)
                 ! xvalu is the velocity

!!$                 xvalu = elbve(idime,inode)                                   ! (b) and (d), velocity
!!$                 if (nofix == ndofn_nsa) xvalu= xvalu *  elbve(ndime+1,inode) ! (a), all
                 do jnode = 1,pnode 
                    do jdofn = 1,ndofn_nsa
                       jevat = (jnode-1)*pdofn + jdofn
                       elmat(ievat,jevat) = 0.0_rp
!!$                       elrhs(jevat) = elrhs(jevat) - elmat(jevat,ievat) * xvalu
                       elmat(jevat,ievat) = 0.0_rp
                    end do
                 end do
                 elmat(ievat,ievat)        = adiag
!!$                 elrhs(ievat)              = adiag * xvalu
                 elrhs(ievat)              = 0.0_rp
                 !              if (kfl_linea_nsa == 2) elsou(ievat) = adiag * xvalu
              end if
           end do


           ! Continuity dof
           if( kfl_elfix(ndime+1,inode) == 1 ) then
              ievat = (inode-1)*pdofn + ndime + 1
              adiag = elmat(ievat,ievat)
!!$              xvalu = elbve(ndime+1,inode)             ! (a) and (c), continuity
              do jnode = 1,pnode 
                 do jdofn = 1,ndofn_nsa
                    jevat = (jnode-1)*pdofn + jdofn
                    elmat(ievat,jevat) = 0.0_rp
!!$                    elrhs(jevat) = elrhs(jevat) - elmat(jevat,ievat) * xvalu
                    elmat(jevat,ievat) = 0.0_rp
                 end do
              end do
              elmat(ievat,ievat)        = adiag
!!$              elrhs(ievat)              = adiag * xvalu
              elrhs(ievat)              = 0.0_rp
              !           if (kfl_linea_nsa == 2) elsou(ievat) = adiag * xvalu
           end if

           ! Energy or Temperature
           if( kfl_elfix(ndime+2,inode) > 0) then
              ievat = (inode-1)*pdofn + ndime + 2
              adiag = elmat(ievat,ievat)
              ! xvalu is the temperature, which is always the physical quantity fixed
!!$              xvalu = elbve(ndime+2,inode)   ! (b), temperature 
!!$              if (nofix == ndofn_nsa) then   ! (a), all (i.e. energy)
!!$                 !
!!$                 !when all dof's are fixed, conservative quantities are used as boundary conditions
!!$                 !
!!$                 xvalu= elbve(ndime+1,inode) * (xvalu * xhecv + 0.5_rp * velsq) ! total energy
!!$              end if
              do jnode = 1,pnode 
                 do jdofn = 1,ndofn_nsa
                    jevat = (jnode-1)*pdofn + jdofn
                    elmat(ievat,jevat) = 0.0_rp
!!$                    elrhs(jevat) = elrhs(jevat) - elmat(jevat,ievat) * xvalu
                    elmat(jevat,ievat) = 0.0_rp
                 end do
              end do
              elmat(ievat,ievat)        = adiag
!!$              elrhs(ievat)              = adiag * xvalu
              elrhs(ievat)              = 0.0_rp
              !           if (kfl_linea_nsa == 2) elsou(ievat) = adiag * xvalu
           end if

        end if

     end if   ! nofix

  end do


end subroutine nsa_assbou
