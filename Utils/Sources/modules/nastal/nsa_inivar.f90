!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_inivar.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   
!> @details 
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_inivar(itask)
  use def_parame
  use def_master
  use def_domain
  use def_nastal
  use def_solver
  implicit none
  integer(ip), intent(in) :: itask  !< Who is calling
  integer(4)              :: istat
  integer(ip)             :: &
       ielty,pgaus,igaus,jgaus,idime,jdime,ifreq,jfreq,kfreq,frequ(3)
  real(rp)                :: disga,arfou,refle

  select case(itask)

  case(0)

     !
     ! Wopos definition
     !
     
     postp(1) % wopos( 1, 1) = 'VELOC'
     postp(1) % wopos( 1, 2) = 'MOMEN'
     postp(1) % wopos( 1, 3) = 'VSHOM'  ! Schock speed from RH for momentum equation
     if(ndime == 3) postp(1) % wopos( 1, 4) = 'VORTI'

     postp(1) % wopos( 2, 1) = 'VECTO'
     postp(1) % wopos( 2, 2) = 'VECTO'
     postp(1) % wopos( 2, 3) = 'VECTO'
     if(ndime == 3) postp(1) % wopos( 2, 4) = 'VECTO'

     postp(1) % wopos( 1,20) = 'PRESS'
     postp(1) % wopos( 1,21) = 'DENSI'
     postp(1) % wopos( 1,22) = 'TEMPE'
     postp(1) % wopos( 1,23) = 'ENERG'
     postp(1) % wopos( 1,24) = 'VISCO'
     postp(1) % wopos( 1,25) = 'MACHN'
     postp(1) % wopos( 1,26) = 'ENTHA'
     postp(1) % wopos( 1,27) = 'TEMCE'  ! temperature in celsius
     postp(1) % wopos( 1,28) = 'BAROT'  ! barotropic coefficient
     if (ndime == 2)  postp(1) % wopos( 1,29) = 'VORTI'  ! in 2D problems, vorticity is a scalar
     postp(1) % wopos( 1,30) = 'MODVE'  ! velocity module
     postp(1) % wopos( 1,31) = 'MODVO'  ! vorticity module
     postp(1) % wopos( 1,32) = 'LAMB2'  ! velocity gradient 2nd invariant: 2ND_INVARIANT_GRAD_VELOCITY  
     postp(1) % wopos( 1,33) = 'XVELO'  ! X-velocity
     postp(1) % wopos( 1,34) = 'YVELO'  ! Y-velocity
     postp(1) % wopos( 1,35) = 'ZVELO'  ! Z-velocity
     postp(1) % wopos( 1,36) = 'TAUTO'  ! total tau
     postp(1) % wopos( 1,37) = 'TAUVI'  ! viscosity tau
     postp(1) % wopos( 1,38) = 'TAUCO'  ! convection tau
     postp(1) % wopos( 1,39) = 'TAUSO'  ! acoustic tau
     postp(1) % wopos( 1,40) = 'DPHYD'  ! p - p_hydrostatic
     postp(1) % wopos( 1,41) = 'DDHYD'  ! rho - rho_hydrostatic
     postp(1) % wopos( 1,42) = 'THBAS'  ! tempe - tempe_reference
     postp(1) % wopos( 1,43) = 'SGSMX'  ! umome sgs
     postp(1) % wopos( 1,44) = 'SGSMY'  ! vmome sgs
     postp(1) % wopos( 1,45) = 'SGSDE'  ! density sgs
     postp(1) % wopos( 1,46) = 'SGSEN'  ! energy or theta sgs
     postp(1) % wopos( 1,47) = 'SGSMO'  ! sgs of momentum as a vector
     postp(1) % wopos( 1,48) = 'XMOME'  ! X-momentum
     postp(1) % wopos( 1,49) = 'YMOME'  ! Y-momentum
     postp(1) % wopos( 1,50) = 'ZMOME'  ! Z-momentum
     postp(1) % wopos( 1,51) = 'FREQU'  ! frequency of nondiagonal VMS
     postp(1) % wopos( 1,52) = 'VSHOC'  ! Schock speed from RH for continuity equation
     postp(1) % wopos( 1,53) = 'VSHOE'  ! Schock speed from RH for energy equation
     postp(1) % wopos( 1,54) = 'NUTUR'  ! Turbulent viscosity
!
     postp(1) % wopos( 1,55) = 'AVVEL'  ! Time-averaged velocity vector
     postp(1) % wopos( 1,56) = 'AVVE2'  ! Time-averaged normal stress Vx*Vx, Vy*Vy & Vz*Vz
     postp(1) % wopos( 1,57) = 'AVVXY'  ! Time-averaged shear stress Vx*Vy, Vx*Vz & Vy*Vz
     postp(1) % wopos( 1,58) = 'AVMOM'  ! Time-averaged momentum ro*Vx, ro*Vy & ro*Vz
     postp(1) % wopos( 1,59) = 'AVPRE'  ! Time-averaged pressure
     postp(1) % wopos( 1,60) = 'AVTEM'  ! Time-averaged temperature     
     postp(1) % wopos( 1,63) = 'RHSMO'
     postp(1) % wopos( 1,64) = 'RHSDE'
     postp(1) % wopos( 1,65) = 'RHSEN'
     postp(1) % wopos( 1,66) = 'DTMOM'
     postp(1) % wopos( 1,67) = 'DTDEN'
     postp(1) % wopos( 1,68) = 'DTENE'
     postp(1) % wopos( 1,69) = 'PCOEF'  ! Pressure coefficient
     postp(1) % wopos( 1,70) = 'GROUP'  ! groups

     
     postp(1) % wopos( 2,20) = 'SCALA'
     postp(1) % wopos( 2,21) = 'SCALA'
     postp(1) % wopos( 2,22) = 'SCALA'
     postp(1) % wopos( 2,23) = 'SCALA'
     postp(1) % wopos( 2,24) = 'SCALA'
     postp(1) % wopos( 2,25) = 'SCALA'
     postp(1) % wopos( 2,26) = 'SCALA'
     postp(1) % wopos( 2,27) = 'SCALA'
     postp(1) % wopos( 2,28) = 'SCALA'
     if (ndime == 2)  postp(1) % wopos( 2,29) = 'SCALA'
     postp(1) % wopos( 2,30) = 'SCALA'  ! velocity module
     postp(1) % wopos( 2,31) = 'SCALA'  ! vorticity module
     postp(1) % wopos( 2,32) = 'SCALA'  ! velocity gradient 2nd invariant: 2ND_INVARIANT_GRAD_VELOCITY  
     postp(1) % wopos( 2,33) = 'SCALA'  ! X-velocity 
     postp(1) % wopos( 2,34) = 'SCALA'  ! Y-velocity 
     postp(1) % wopos( 2,35) = 'SCALA'  ! Z-velocity 
     postp(1) % wopos( 2,36) = 'SCALA'  ! total tau
     postp(1) % wopos( 2,37) = 'SCALA'  ! viscosity tau
     postp(1) % wopos( 2,38) = 'SCALA'  ! convection tau
     postp(1) % wopos( 2,39) = 'SCALA'  ! acoustic tau
     postp(1) % wopos( 2,40) = 'SCALA'  ! acoustic tau
     postp(1) % wopos( 2,41) = 'SCALA'  ! acoustic tau
     postp(1) % wopos( 2,42) = 'SCALA'  ! acoustic tau
     postp(1) % wopos( 2,43) = 'SCALA'  ! umomentum sgs
     postp(1) % wopos( 2,44) = 'SCALA'  ! vmomentum sgs
     postp(1) % wopos( 2,45) = 'SCALA'  ! density sgs
     postp(1) % wopos( 2,46) = 'SCALA'  ! energy or theta sgs
     postp(1) % wopos( 2,47) = 'VECTO'  ! sgs of momentum as a vector
     postp(1) % wopos( 2,48) = 'SCALA'  ! X-momentum
     postp(1) % wopos( 2,49) = 'SCALA'  ! Y-momentum
     postp(1) % wopos( 2,50) = 'SCALA'  ! Z-momentum
     postp(1) % wopos( 2,51) = 'SCALA'  ! frequency of nondiagonal VMS
     postp(1) % wopos( 2,52) = 'SCALA'  ! VSHOC
     postp(1) % wopos( 2,53) = 'SCALA'  ! VSHOE
     postp(1) % wopos( 2,54) = 'SCALA'  ! Turbulent viscosity
!
     postp(1) % wopos( 2,55) = 'VECTO'  ! Time-averaged velocity vector
     postp(1) % wopos( 2,56) = 'VECTO'  ! Time-averaged normal stress Vx*Vx, Vy*Vy & Vz*Vz
     postp(1) % wopos( 2,57) = 'VECTO'  ! Time-averaged shear  stress Vx*Vy, Vx*Vz & Vy*Vz
     postp(1) % wopos( 2,58) = 'VECTO'  ! Time-averaged momentum ro*Vx, ro*Vy & ro*Vz
     postp(1) % wopos( 2,59) = 'SCALA'  ! Time-averaged pressure
     postp(1) % wopos( 2,60) = 'SCALA'  ! Time-averaged temperature

    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------! 
    !../../kernel/defmod/def_kintyp.f90. nvarp = postprocess variables  
    postp(1) % wopos( 1,61) = 'CHA01'  !> LODI      
    postp(1) % wopos( 2,61) = 'VECTO'

    postp(1) % wopos( 1,62) = 'CHA02'  !> SCHLIEREN SCHLIEREN  
    postp(1) % wopos( 2,62) = 'VECTO'

    postp(1) % wopos( 1,63) = 'CHA03'  !> SCHLIEREN SCHLIEREN  
    postp(1) % wopos( 2,63) = 'VECTO'

    postp(1) % wopos( 1,64) = 'CHA04'  !> SCHLIEREN SCHLIEREN  
    postp(1) % wopos( 2,64) = 'VECTO'

    postp(1) % wopos( 1,65) = 'CHA05'  !> SCHLIEREN SCHLIEREN  
    postp(1) % wopos( 2,65) = 'VECTO'

    !-----------------------------------------------------------------------! 
!     postp(1) % wopos( 2,63) = 'SCALA'
!     postp(1) % wopos( 2,64) = 'SCALA'
!     postp(1) % wopos( 2,65) = 'SCALA'

     postp(1) % wopos( 2,66) = 'SCALA'
     postp(1) % wopos( 2,67) = 'SCALA'
     postp(1) % wopos( 2,68) = 'SCALA'

     postp(1) % wopos( 2,69) = 'SCALA'

     postp(1) % wopos( 2,70) = 'SCALA'

     postp(1) % wopos( 1,36) = '2DPRO'
     postp(1) % wopos( 2,36) = 'DLINE'
     
     !
     ! Sets variables
     !
     postp(1) % woese (  1) = 'VELOC'
     postp(1) % woese (  2) = 'VORTI'
     
     postp(1) % wonse (  1) = 'VELOX'
     postp(1) % wonse (  2) = 'VELOY'
     postp(1) % wonse (  3) = 'VELOZ'    
     postp(1) % wonse (  4) = 'PRESS'

     postp(1) % wobse (  1) = 'MASS '     ! Mass rho*u.n
     postp(1) % wobse (  2) = 'DENSI'     ! Density rho.n
     postp(1) % wobse (  3) = 'UBULK'     ! Bulk velocity rho*u.n / rho.n
     !
     ! Witness variables 
     !
     postp(1) % wowit (1)     = 'VELOX'   ! x-velocity
     postp(1) % wowit (2)     = 'VELOY'   ! y-velocity
     postp(1) % wowit (3)     = 'VELOZ'   ! z-velocity
     postp(1) % wowit (4)     = 'PRESS'   ! Pressure
     postp(1) % wowit (5)     = 'TEMPE'   ! Temperature
     !
     ! Unknowns element-wise dimensions
     !
     ndofn_nsa = ndime+2
     ndof2_nsa = ndofn_nsa*ndofn_nsa
     nevat_nsa = ndofn_nsa*mnode
     nevab_nsa = ndime*mnode
     nflub_nsa = 3        

     ndtdf_nsa = 3                                     ! Forced compressible flow (old kfl_foreg_nsa = 0)
     nkeep_nsa = ndofn_nsa + 2 + 3                     ! Base state (i.e. for KEEPHYDRO)

     ! Derived dimensions according to the general algorithm type
     !

     nunkn_nsa = ndofn_nsa *  npoin

     !
     ! Solvers
     !     
     call soldef(-1_ip)             ! only one problem to solve
     solve(1) % wprob     = 'PHYSICAL_VARIABLES'
     solve(1) % ndofn     = ndofn_nsa
     solve(1) % kfl_solve = 1                 ! Output flag

     ! In the case of implicit, these values will be defined in nsa_reanut, when calling reasol
     ! These are default values, corresponding to explicit with global time step
     solve(1) % kfl_algso = SOL_SOLVER_RICHARDSON
     solve(1) % kfl_preco = SOL_CLOSE_MASS_MATRIX   
!!solve(1) % kfl_preco = SOL_MASS_MATRIX

     ! Default: when imposing bc, modify matrix and rhs from at element level
     solve(1) % kfl_iffix = 0

     kfl_refre_nsa = 0

     izone_nsa= lzone(ID_NASTAL)

     resou_nsa = 1.0_rp
     resin_nsa = 1.0_rp

     !
     ! Nullify pointers
     !
     nullify(avvel_nsa)
     nullify(avve2_nsa)
     nullify(avvxy_nsa)
     nullify(avmom_nsa)
     nullify(avpre_nsa)
     nullify(avtem_nsa)

 case(1)   

     if (ndime==2) then
        nindx_nsa(1,1)= 1
        nindx_nsa(2,2)= 2
        nindx_nsa(1,2)= 3
        nindx_nsa(2,1)= 3
     else if (ndime==3) then        
        nindx_nsa(1,1)= 1
        nindx_nsa(2,2)= 2
        nindx_nsa(3,3)= 3
        nindx_nsa(1,2)= 4
        nindx_nsa(1,3)= 5
        nindx_nsa(2,3)= 6
        nindx_nsa(2,1)= 4
        nindx_nsa(3,1)= 5
        nindx_nsa(3,2)= 6
     end if

     !
     ! Number of time advance components 
     !
     if(kfl_timei_nsa==1) then
        kfl_timei=1
!!        kfl_tiacc_nsa = 1             ! in nastal this is the default "time accuracy" 

        ! kfl_timul_nsa has dissapeared
        ncomp_nsa=2+kfl_tiacc_nsa  ! needs previous time step fields to compute time residual           


!!$        if(kfl_timul_nsa==0) then
!!$           !
!!$           ! low-order mono-step
!!$           !
!!$           ncomp_nsa=2+kfl_tiacc_nsa  ! needs previous time step fields to compute time residual           
!!$        else if(kfl_timul_nsa==1) then
!!$           !
!!$           ! high-order mono-step (RK, FRK, CN)
!!$           !
!!$           ncomp_nsa=2+kfl_tiacc_nsa  ! needs previous time step fields to compute time residual           
!!$
!!$           call nsa_partim            ! set the scheme parameters
!!$
!!$        else if(kfl_timul_nsa==2) then
!!$           !
!!$           ! BDF scheme, variables-based multi-step
!!$           !
!!$           ncomp_nsa=2+kfl_tiacc_nsa
!!$           call runend ('NSA_INIVAR: TO BE PROGRAMMED IN IMPLICIT SCHEMES')
!!$
!!$           call nsa_partim            ! set the scheme parameters
!!$
!!$        else if(kfl_timul_nsa==3) then
!!$           !
!!$           ! residuals-based multi-step
!!$           !
!!$           ncomp_nsa=2+kfl_tiacc_nsa  ! needs previous time step fields to compute time residual           
!!$           call nsa_partim            ! set the scheme parameters
!!$
!!$        end if
     else

        call runend('NSA_INIVAR: ONLY TRANSIENT PROBLEMS')

     end if
     !
     ! Time accuracy: save original value
     !
     kfl_tiaor_nsa=kfl_tiacc_nsa

  case(2)   

     !
     ! Cosine - sine coefficients matrices
     !
     
!     write(6,*)

!!$     do ielty = 1,nelty
!!$        if( lexis(ielty) == 1 ) then 
!!$           pgaus=ngaus(ielty)
!!$           do igaus= 1,pgaus
!!$              do jgaus=1,pgaus
!!$                 do ifreq= - mfreq_nsa , mfreq_nsa
!!$                    frequ(1) = ifreq 
!!$                    do jfreq= - mfreq_nsa , mfreq_nsa 
!!$                       frequ(2) = jfreq 
!!$                       do kfreq= - mfreq_nsa , mfreq_nsa
!!$                          frequ(3) = kfreq 
!!$                          arfou= 0.0_rp
!!$                          do idime=1,ndime
!!$                             disga= elmar(ielty)%posgp(idime,igaus) - elmar(ielty)%posgp(idime,jgaus) 
!!$                             arfou= arfou + pi * real(frequ(idime)) * disga
!!$                          end do
!!$
!!$                          cosma_nsa(igaus,jgaus,&
!!$                               mfreq_nsa+1+frequ(1),mfreq_nsa+1+frequ(2),mfreq_nsa+1+frequ(3),ielty) &
!!$                               = cos(arfou) 
!!$                          sinma_nsa(igaus,jgaus,&
!!$                               mfreq_nsa+1+frequ(1),mfreq_nsa+1+frequ(2),mfreq_nsa+1+frequ(3),ielty) &
!!$                               = sin(arfou) 
!!$print*
!!$print*, 'igaus jgaus ifreq jfreq kfreq cosma sinma', igaus,jgaus,ifreq,jfreq,kfreq,cos(arfou),sin(arfou)
!!$                        
!!$                       end do
!!$                    end do
!!$                 end do
!!$              end do
!!$           end do
!!$        end if
!!$     end do


     do ielty = 1,nelty
        if( lexis(ielty) == 1 ) then 
           pgaus=ngaus(ielty)
           do igaus= 1,pgaus
              do jgaus=1,pgaus
                 do ifreq= - frmax_nsa , frmax_nsa
                    frequ(1) = ifreq 
                    do jfreq= - frmax_nsa , frmax_nsa 
                       frequ(2) = jfreq 
                       do kfreq= - frmax_nsa , frmax_nsa
                          frequ(3) = kfreq 
                          arfou= 0.0_rp
                          do idime=1,ndime
                             disga= elmar(ielty)%posgp(idime,igaus) - elmar(ielty)%posgp(idime,jgaus) 
                             arfou= arfou + 2.0_rp * pi * real(frequ(idime)) * disga / hnatu(ielty)
                          end do

                          cosma_nsa(igaus,jgaus,&
                               frmax_nsa+1+frequ(1),frmax_nsa+1+frequ(2),frmax_nsa+1+frequ(3),ielty) &
                               = cos(arfou) 
                          sinma_nsa(igaus,jgaus,&
                               frmax_nsa+1+frequ(1),frmax_nsa+1+frequ(2),frmax_nsa+1+frequ(3),ielty) &
                               = sin(arfou)                        
                      end do
                    end do
                 end do
              end do
           end do
        end if
     end do

  end select
  
end subroutine nsa_inivar
