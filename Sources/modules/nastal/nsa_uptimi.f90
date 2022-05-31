!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_uptimi.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   This routine manages the high order time integration schemes
!> @details This routine manages the high order time integration schemes
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_uptimi(itask)
  use      def_master
  use      def_parame
  use      def_domain
  use      def_nastal
  implicit none
  integer(ip), intent(in) :: itask   !< Who is calling
  integer(ip) :: ipoin,idofn,ievat,jevat,idime,jtinn,ktinn,ifrap
  real(rp)    :: ufact,rauxi,rauti

  

  if( IMASTER ) return

  select case(itask)

  case(one)     ! -->> SOLITE is calling

     if (itinn(modul)>1) return

!!$     if (kfl_timul_nsa /= 0) then
!!$        do jtinn= 1,nromp_nsa+1
!!$           do ipoin= 1,npoin
!!$              do idofn= 1,ndofn_nsa
!!$                 ievat = (ipoin-1)*ndofn_nsa + idofn
!!$                 do ifrap= 1,nfrap_nsa
!!$                    rhsax_nsa(ievat,jtinn,ifrap) =  0.0_rp            ! A-residual
!!$                    rhsax_nsa(ievat,jtinn,ifrap) =  0.0_rp            ! b-residual
!!$                 end do
!!$              end do
!!$           end do
!!$        end do
!!$     end if
     
  case(two)     ! -->> UPFMOM is calling


  case(three)   ! -->> UPFMOM is calling
     

  case(four)    ! -->> UPCONS is calling

     jtinn = itinn(modul)
     

     if(kfl_tisch_nsa == 41) then
        !
        ! RK-3:  we store the following variables that will be used by the RK solver:
        !
        !
        !  sum(alpha_ik * U(k) <--- 
        !
        !
        do ipoin= 1,npoin
           do idofn= 1,ndofn_nsa
              ievat = (ipoin-1)*ndofn_nsa + idofn
              !
              !    Storing unknown at each jtinn iteration
              !
              unkit_nsa(ievat,jtinn) = unkno(ievat) !u(jtinn -1)


              !
              !    Adding up beta-residual
              !
              rhsid(ievat) = parkb_nsa(jtinn,1) * rhsid(ievat)


              !
              !    Adding up alpha-residual
              !
              unkna_nsa(ievat) = 0.0_rp
              do ktinn = 1, jtinn
                 unkna_nsa(ievat) = unkna_nsa(ievat) + parka_nsa(jtinn,ktinn,1) * unkit_nsa(ievat,ktinn)
              end do
              

           end do
        end do

     else
        !
        ! Heun 3 (maybe, ask Mariano)
        !
        do ipoin= 1,npoin
           do idofn= 1,ndofn_nsa
              ievat = (ipoin-1)*ndofn_nsa + idofn
              !
              !    Storing each jtinn residual
              !
              rhsax_nsa(ievat,jtinn,1)= rhsid(ievat  )
              !
              !    Adding up b-residual
              !
              rhsax_nsa(ievat,nromp_nsa+1,1)=  &
                   rhsax_nsa(ievat,nromp_nsa+1,1) + parkb_nsa(jtinn,1) * rhsid(ievat  )

              if (jtinn == miinn_nsa) then
                 !
                 !    Last step: rauxi is the final b-residual
                 !
                 rauxi= rhsax_nsa(ievat,nromp_nsa+1,1)
              else
                 !
                 !    Adding up A-residual
                 !
                 rauxi= 0.0_rp
                 do ktinn= 1, jtinn
                    rauxi=  rauxi + parka_nsa(ktinn,jtinn,1) * rhsax_nsa(ievat,ktinn,1)
                 end do
              end if
              !
              !    Residual assignement
              !
              rhsid(ievat  )= rauxi
           end do
        end do
     end if
  end select

end subroutine nsa_uptimi
