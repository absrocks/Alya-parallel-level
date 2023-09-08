!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_coupli.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   ALEFOR coupling
!> @details ALEFOR coupling:
!>    FMALE coupling at n+1.
!>    In Parallel, some checks are necessary because of the following:
!>    - Different subdomains can find the host element of the previous
!>      position of node IPOIN, because of tolerance: GISCA(IPOIN)=IELEM
!>    - They agree (according to the rank) on who computes it. The one
!>      that owns the node has GISCA(IPOIN)=ILEME; the toher put
!>      GISCA(IPOIN)=0
!>    - Finally, the previous value u^n is set to zero when GISCA(IPOIN)=0
!>      and then the value are exchanged in the usual way between slaves
!> @} 
!-----------------------------------------------------------------------
subroutine ale_coupli(itask)
  use def_master
  use def_kermod
  use def_domain
  use def_elmtyp
  use def_alefor
  use mod_intpol

  implicit none
  integer(ip), intent(in) :: itask !> where the subroutine is called 
  integer(ip)             :: idime,ipoin,inode,jpoin
  real(rp)                :: coor1(ndime),coor2(ndime,npoin_2),tempo
  real(rp),pointer        :: shapl(:)

  select case ( itask )

  case ( ITASK_ENDITE )

     !-------------------------------------------------------------------
     !
     ! ENDITE
     !
     !-------------------------------------------------------------------
     !
     ! Coupling with modules
     !
     if(    (kfl_coupl(ID_NASTIN,ID_ALEFOR) /= 0 .or. kfl_coupl(ID_TEMPER,ID_ALEFOR) /= 0) & 
          .and. (associated(lnint) .and. ittim >= 1) ) then

        !call memgen(0_ip,ndime,npoin_2)
        !do ipoin = 1,npoin
        !   do idime = 1,ndime
        !      gevec(idime,ipoin) = dispm(idime,ipoin,1)
        !   end do
        !end do

        !call pararr('FRI',NPOIN_TYPE,ndime*npoin_2,gevec)

        !do ipoin = 1,npoin_2
        !   if (lntib(ipoin) < 0) then
        !      do idime = 1,ndime
        !         coor2(idime,ipoin) = coord(idime,ipoin) + gevec(idime,ipoin)           
        !         coord(idime,ipoin) = coord(idime,ipoin) 
        !      end do
        !   end if
        !end do
        !call memgen(2_ip,ndime,npoin_2)

        if( kfl_coupl(ID_NASTIN,ID_ALEFOR) /= 0 ) then                       
           call memgen(0_ip,ndime,npoin)
        elseif( kfl_coupl(ID_TEMPER,ID_ALEFOR) /= 0 ) then
           call memgen(0_ip,npoin,0_ip)
        end if

        do ipoin = 1,npoin
           if (lntib(ipoin) <= 0 .and. lpoty(ipoin) == 0) then
              do idime = 1,ndime
                coor1(idime) = coord(idime,ipoin) + dispm(idime,ipoin,1)
              end do

              allocate(shapl(lnin2(ipoin) % limit+1))
              call krigin(coord,coor1,lnin2(ipoin) % limit,lnin2(ipoin) % lnode,shapl)
              !
              ! NASTIN
              !
              if( kfl_coupl(ID_NASTIN,ID_ALEFOR) /= 0 ) then
                 do idime = 1,ndime
                    gevec(idime,ipoin) = 0.0_rp
                 end do                 
                 !if (ipoin <= npoi1 .or. (ipoin >= npoi2 .and. ipoin <= npoi3)) then
                 !   do idime = 1,ndime
                 !      gevec(idime,ipoin) = gevec(idime,ipoin) + shapl(lnin2(ipoin) % limit+1) &
                 !           * veloc(idime,ipoin,3)
                 !   end do                    
                 !end if
                 !
                 ! We take into account the parallel reduction in nodes that share more than one subdomain
                 !
                 do inode = 1,lnin2(ipoin) %limit
                    jpoin = lnin2(ipoin) %lnode(inode)                    
                    if (jpoin <= npoin) then
                       do idime = 1,ndime
                          gevec(idime,ipoin) = gevec(idime,ipoin) + lnin2(ipoin) % shapl(inode) & 
                               * shapl(inode) * veloc(idime,jpoin,3)
                       end do
                    end if
                 end do
              end if
              !
              ! TEMPER
              !
              if( kfl_coupl(ID_TEMPER,ID_ALEFOR) /= 0 ) then                
                 gesca(ipoin) = 0.0_rp                
                 !if (ipoin <= npoi1 .or. (ipoin >= npoi2 .and. ipoin <= npoi3)) then
                 !   gesca(ipoin) = gesca(ipoin) + shapl(lnin2(ipoin) % limit+1) &
                 !        * tempe(ipoin,3)
                 !end if
                 !
                 ! We take into account the parallel reduction in nodes that share more than one subdomain
                 !
                 do inode = 1,lnin2(ipoin) % limit
                    jpoin = lnin2(ipoin) %lnode(inode)        
                    if (jpoin <= npoin) then                                        
                       gesca(ipoin) = gesca(ipoin) + lnin2(ipoin) % shapl(inode) &
                            * shapl(inode) * tempe(jpoin,3)
                    end if
                 end do
              end if
              deallocate(shapl)
           end if
        end do
        !
        ! Recover the original values for coord
        !
        !do ipoin = 1,npoin_2
        !   if (lntib(ipoin) < 0) then
        !      do idime = 1,ndime                      
        !         coord(idime,ipoin) = coor2(idime,ipoin)
        !      end do
        !   end if
        !end do
        !
        ! Interchange information between subdomains
        !      
        if( kfl_coupl(ID_NASTIN,ID_ALEFOR) /= 0 ) then
           call pararr('SLX',NPOIN_TYPE,(ndime)*npoin,gevec)
        end if
        if( kfl_coupl(ID_TEMPER,ID_ALEFOR) /= 0 ) then
           call pararr('SLX',NPOIN_TYPE,npoin,gesca)
        end if



        do ipoin = 1,npoin
           if (lntib(ipoin) == 0 .and. lpoty(ipoin) == 0) then
              !
              ! NASTIN
              !
              if( kfl_coupl(ID_NASTIN,ID_ALEFOR) /= 0 ) then
                 do idime = 1,ndime
                    veloc(idime,ipoin,3) = gevec(idime,ipoin)
                 end do
              end if
              !
              ! TEMPER
              !
              if( kfl_coupl(ID_TEMPER,ID_ALEFOR) /= 0 ) then
                 tempe(ipoin,3) = gesca(ipoin)
              end if
           end if
        end do

        if( kfl_coupl(ID_NASTIN,ID_ALEFOR) /= 0 ) then
           call memgen(2_ip,ndime,npoin)
        end if
        if( kfl_coupl(ID_TEMPER,ID_ALEFOR) /= 0 ) then
           call memgen(2_ip,npoin,0_ip)
        end if

        do ipoin = 1,npoin
           if (lntib(ipoin) < 0_ip .and. lnti2(ipoin) > 0 ) then
              print *,"--| ALYA     THERE IS A NEW NODE: ",ipoin
           end if
        end do
        
     end if
  end select

!===============================================================| contains |===!
contains
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date
  !> @brief
  !> @details
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_ale_concou(n_iters_ale) ! unused??
  !< For global iterations in coupling schemes
  use def_alefor, only: coord_ale
  use def_master, only: dispm
  implicit none
  integer(ip) :: counter = 0_ip !!??
  integer(ip) :: n_iters_ale
  !
  if(kfl_gocou == 1) then
    !
    if(INOTMASTER) then                                 !< ale_updunk.case(6_ip)
      counter = counter + 1_ip
      !
      if(counter == n_iters_ale) then
        counter = 0_ip
        return
      end if
      !
      dispm(    1:ndime,1:npoin,1) = dispm(    1:ndime,1:npoin,3)
      coord_ale(1:ndime,1:npoin,1) = coord_ale(1:ndime,1:npoin,3)
      coord(    1:ndime,1:npoin  ) = coord_ale(1:ndime,1:npoin,1)
    end if
    !

  endif
  !
  if(.not.INOTMASTER) then
    print *, "[commdom_alya_ale_concou] "
  endif
  !
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date
  !> @brief
  !> @details
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_ale_doiter( kfl_bvess_ale, kfl_solve_ale )
  ! def_alefor.f90:
  !   kfl_solve_ale                                                 !< ALE solved
  ! ale_doiter.f90:
  !   kfl_bvess_ale = 1_ip                                          !< FSI
  implicit none
  integer(ip), intent(out) :: kfl_bvess_ale, kfl_solve_ale
     kfl_bvess_ale = 1
     kfl_solve_ale = 1
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date
  !> @brief
  !> @details
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_ale_outvar()
  use def_master, only: dispm, gevec
  implicit none
      gevec => dispm(:,:,1)
  end subroutine
  !-----------------------------------------------------------------------||---!

!===============================================================| contains |===!
end subroutine ale_coupli
 
