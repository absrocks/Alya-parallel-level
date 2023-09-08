!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_iniunk.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   This routine sets up the initial condition for the mesh velocity
!> @details This routine sets up the initial condition for the mesh velocity
!> @} 
!-----------------------------------------------------------------------
subroutine ale_iniunk()
  use def_master
  use def_domain
  use def_alefor
  use def_kermod, only : kfl_adj_prob,sens_mesh
  
  implicit none

  integer(ip)             :: kpoin,idime,iimbo,ipoin

  if(kfl_rigid_ale == 1) call ale_inirbo()
  !
  ! Reading mesh sensitivities for optimization
  !
  if( INOTMASTER ) then
    if (kfl_adj_prob == 1_ip) then
      do idime = 1, ndime
        do ipoin=1,npoin
          sens_mesh(idime,ipoin) = xfiel(-kfl_sensi_ale) %a(idime,ipoin,1)
        end do
      enddo
    end if
  endif
  
  if( kfl_rstar == 0 ) then  

     if (moddi_ale < 0) then
        do ipoin=1,npoin
           do idime=1,ndime
              coord(idime,ipoin) = coord(idime,ipoin) + xfiel(-moddi_ale)%a(idime,ipoin,1)
!              coord_ale(idime,ipoin,1) = coord_ale(idime,ipoin,1) + xfiel(-moddi_ale)%a(idime,ipoin)
!              dispm(idime,ipoin,1:3) = xfiel(-moddi_ale)%a(idime,ipoin)
           end do
        end do
     end if

     !
     ! If ALE is only used for smoothing purpose
     !
     if( coupling('ALEFOR','SOLIDZ') == 0 ) then
        call ale_doiter() 
     end if

  else
     
     call ale_restar(1_ip)

     if(kfl_rigid_ale == 1) then

        ! call ale_restar(1_ip)
        call ale_parall(2_ip)

        if( INOTMASTER ) then
           !
           ! obtain  % new coord using  coord & dispm
           !
           do ipoin = 1,npoin
              do idime = 1,ndime    
                 coord(idime,ipoin) = coord(idime,ipoin) + dispm(idime,ipoin,1)
                 coord_ale(idime,ipoin,1) = coord(idime,ipoin)
              end do
           end do
           !
           ! reobtain  % cooib using  coord & dispm
           !
           do iimbo = 1,nrbod
              do kpoin = 1,rbbou(iimbo) % npoib
                 ipoin = rbbou(iimbo) % lninv(kpoin)
                 do idime = 1,ndime    
                    !                 rbbou(iimbo) % cooib(idime,kpoin) = coord(idime,ipoin) + dispm(idime,ipoin,1)
  
                    rbbou(iimbo) % cooib(idime,kpoin) = coord(idime,ipoin)
                 end do
              end do
           end do
        end if
        !
        ! something similar to   call nsi_updunk(11_ip)                  ! VELOC(:,:,1) & VELOC(:,:,2)   <= VELOC(:,:,nprev_nsi)
        !      
        call ale_updunk(31_ip)

     end if
  end if

  call mescek(1_ip)

end subroutine ale_iniunk
