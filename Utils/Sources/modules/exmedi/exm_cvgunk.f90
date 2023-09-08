!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_cvgunk.f90
!> @date    16/11/1966
!> @author  Mariano Vazquez
!> @brief   Unknowns convergence
!> @details Unknowns convergence
!> @}
!------------------------------------------------------------------------
subroutine exm_cvgunk(itask)
!-----------------------------------------------------------------------
!****f* Exmedi/exm_cvgunk
! NAME 
!    exm_cvgunk
! DESCRIPTION
!    This routine performs several convergence checks for the 
!    incompressible NS equations
! USES
!    exm_endite (itask=1,2) (TO BE UNPROGRAMMED)
!    exm_endste (itask=3)
! USED BY
!    Exmedi
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_exmedi
  use mod_communications , only: PAR_SUM  
  implicit none
  integer(ip), intent(in) :: itask   !> itask number
  integer(ip), save       :: ipass=0
  integer(ip)             :: idofn,iroot
  real(rp)                :: numer,denom,riexm(4),cpnew,cpdif
!
! Initializations
!
  numer = 0.0_rp
  denom = 0.0_rp


  select case(itask)


  case(1)
     
     !
     ! Check convergence of the inner (sub)iterations:
     ! || u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||
     !
     
     riexm=0.0_rp

     call residu(kfl_normc_exm,one,1_ip,unkno,elmag    ,one,  one,one,1.0_rp,riexm(1))
     
     riexm(4) = riexm(1) + riexm(2)

     if(riexm(4)<cotol_exm.or.itinn(modul)>=miinn_exm) kfl_goite_exm = 0

!!!     call minmax(one,npoin,tempe,vamxm_nsa(1,3),vamxm_nsa(2,3))

     if(itinn(modul) == 1) then
        do idofn= 1,4
           resin_first_exm(idofn) = 1.0_rp
           if (riexm(idofn) > zeexm) then ! to avoid division by zero
              resin_first_exm(idofn) = riexm(idofn)
           end if
        end do
     end if

     if (kfl_refre_exm(itask) == 0) then
        do idofn=1,4
           if (riexm(idofn) > zeexm) then ! to avoid division by zero
              resin_exm(idofn) = riexm(idofn)
           else
              resin_exm(idofn) = zeexm
           end if
        end do
        kfl_refre_exm(itask) = 1
     end if

     !
     ! check stopping condition for subiterations
     !
     if (kfl_adres_exm == 0) then
        if(riexm(4)<cotol_exm.or.itinn(modul)>=miinn_exm) kfl_goite_exm = 0
     else
        if(itinn(modul)>=miinn_exm) then
           kfl_goite_exm = 0
        else
           if (riexm(4) < cotol_exm) kfl_goite_exm = 0                     ! firstly, absolute check
           if (riexm(4)/resin_first_exm(4) < corat_exm) kfl_goite_exm = 0  ! then, relative check
        end if
     end if


     !
     ! Sum up slaves contribution for volelec_exm and pseudecg_exm
     !
!!!     call pararr('SUM',0_ip,1_ip,volelec_exm)
     call PAR_SUM(volelec_exm,'IN MY CODE')
     call PAR_SUM(nrootecg_exm,pseudecg_exm,'IN MY CODE')


     if(INOTSLAVE) then

        if(ipass==0.and.kfl_rstar/=2) then
           ipass=1
           write(momod(modul)%lun_conve,200)
           write(         lun_vinte_exm,300)
        end if

        ! cpold_exm is initialized when defined in def_exmedi, it is meaningful only for the master and sequential
        ! cases
        call cputim(cpnew)
        cpdif = cpnew - cpold_exm
        cpold_exm = cpnew
!
!       Write convergence files
!
        if(dtinv_exm.eq.0.0_rp) dtinv_exm=1.0_rp ! In it=0 gives NaN as parameter

        write(momod(modul)%lun_conve,1000) &
             ittim,&
             itcou,&
             itinn(modul),&
             last_iters_exm,&
             cutim,&
             riexm(4),riexm(1),riexm(2),riexm(3),&
             1.0_rp/dtinv_exm,&
             cpdif,&
             volelec_exm
        flush(momod(modul)%lun_conve)
!!        write(lun_maxmi_nsa,1000) ittim,itcou,itinn(modul),&
!!             vamxm_nsa(1,1),vamxm_nsa(2,1),vamxm_nsa(1,2),vamxm_nsa(2,2),&
!!             vamxm_nsa(1,3),vamxm_nsa(2,3),vamxm_nsa(1,4),vamxm_nsa(2,4)
!!        flush(lun_maxmi_nsa)

     end if


!
! Check convergence of the outer iterations in the norm selected by the user:
! || u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||
!
  case(2)


     call residu(kfl_normc_exm,one,1_ip,unkno,elmag    ,one,  one,one,1.0_rp,riexm(1))
     if (kfl_cemod_exm == 2) call residu(kfl_normc_exm,one,1_ip,unkno,elmag    ,one,one+1,one,1.0_rp,riexm(2))
!!     call residu(kfl_normc_exm,one,ndofn_exm,unkno,refhn_exm,one,  one,one,1.0_rp,riexm(3))

     if (mod(ittim,nvint_exm) == 0) then
        write(         lun_vinte_exm,2000) &
             ittim,&
             cutim,&
             volelec_exm,&
             pseudecg_exm(1:nrootecg_exm)
     end if

!
! Check residual of the time evolution, always in the L2 norm:
! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||
!

  end select
!
! Formats
!
200 format('# ',' --|  Exmedi Convergence File '       ,/,&
         & '# ',' --|  Columns displayed:' ,/,&
         & '# ',' --|  1. Time Step   2. Global Iteration   3. Sub Iteration   4. Solver iterations'  ,/,&
         & '# ',' --|  5. Current time  ',/,&
         & '# ',' --|  6. Total Res.  7. Intra. Res.        8. Extra. Res.       9. Recov. Res.' ,/,&
         & '# ',' --|  10. Critical Time Step  11. CPU time        12. Vol. integral Phi' ,/,&
         & '# ',' --|  ') 

300 format('# ',' --|  Exmedi Volume Integrals File '       ,/,&
         & '# ',' --|  Columns displayed:' ,/,&
         & '# ',' --|  1. Time Step   '  ,/,&
         & '# ',' --|  2. Current time  ',/,&
         & '# ',' --|  3. Vol. integral Phi' ,/,&
         & '# ',' --|  4:nroot. Pseudo ECG' ,/,&
         & '# ',' --|  ') 

1000 format(4x,4(i9,2x),15(2x,e13.6))
2000 format(4x,1(i9,2x),100(2x,e13.6))
        
end subroutine exm_cvgunk


