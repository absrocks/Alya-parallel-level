  !----------------------------------------------------------------------
  !> @addtogroup Nastin
  !> @{
  !> @file    nsi_fsiexch.f90
  !> @author  J.C. Cajas
  !> @date    17/04/2014
  !> @brief   Send force to Solidz
  !> @details Send force to Solidz using the coupling structures and functions.
  !> @} 
  !----------------------------------------------------------------------

subroutine nsa_fsiexch(itask)
  use def_master
  use def_kermod
  use def_domain
  use def_elmtyp
  use def_nastal
  use mod_ker_proper 
  use def_coupli,        only :  coupling_type
  use mod_couplings,     only :  COU_INTERPOLATE_NODAL_VALUES

  implicit none
  integer(ip), intent(in)     :: itask
  integer(ip)                 :: idofn
  integer(ip)                 :: ipoin
  integer(ip)                 :: jpoin
  integer(ip)                 :: idime
  integer(ip)                 :: jdime
  integer(ip)                 :: izdom
  integer(ip)                 :: jzdom

  real(rp),    pointer        :: xvalu(:,:)
  real(rp),    pointer        :: svalu(:,:)
  real(rp),    pointer        :: force(:,:)
  real(rp)                    :: foref

  nullify(xvalu) 
  nullify(force)

  if ( INOTMASTER ) then

     allocate(xvalu(ndime,npoin))
     allocate(svalu(ndime,npoin))
     allocate(force(ndime,npoin))

  else

     allocate(xvalu(1_ip,1_ip))
     allocate(svalu(1_ip,1_ip))
     allocate(force(1_ip,1_ip))

  end if

  if ( itask == 1_ip ) then ! Send force to Solidz

     ! Arguments for the coupling function: 
     ! COU_INTERPOLATE_NODAL_VALUES(coupling label, number of dimensions of the variable to interpolate, 
     ! array to store the results, variable to interpolate )
     ! The setup of the communicators, coordinates and everything else is done in COU_INITIALIZE_COUPLING
     force = 0_rp
     if ( INOTMASTER ) then
!         call nsa_force
        do ipoin = 1,npoin 
           if (kfl_codno(1,ipoin) == 16_ip) then
              do idime = 1,ndime
                 
                !  force(idime,ipoin) = intfo_nsi(ipoin) % bu(idime)
                !  jzdom = 0
                !  do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                    
                !     jpoin = c_dom(izdom)
                !     jzdom = jzdom + 1
                !     do jdime = 1,ndime

                !        force(idime,ipoin) = force(idime,ipoin) - intfo_nsi(ipoin) % Auu(jdime,idime,jzdom) * veloc(jdime,jpoin,1) 
                       
                !     end do
                !     force(idime,ipoin) = force(idime,ipoin) - intfo_nsi(ipoin) % Aup(idime,jzdom) * press(jpoin,1)

                !  end do
                 force(idime,ipoin) = 1.0e-3_rp

              end do
           end if
        end do

     end if
!      print*, "DEBUG: force= ", force(:,1:10)         
     ! force
!      print*, "DEBUG: wet points= ", coupling_type(2_ip) % geome % npoin_wet 
     call COU_INTERPOLATE_NODAL_VALUES(1_ip,ndime,xvalu,force)  

  end if

  if ( itask == 2 ) then ! Receive displacement from Solidz

     xvalu = 0_rp
     ! displacement
     call COU_INTERPOLATE_NODAL_VALUES(2_ip,ndime,xvalu,svalu)

     if ( INOTMASTER ) then

        bvess_ale = 0_rp
        if( associated(xvalu) ) then
           
           do ipoin = 1, npoin
              
!               if (kfl_codno(1,ipoin) == 16_ip) then
              
                 do idime = 1, ndime

                    bvess_ale(idime,ipoin) = xvalu(idime,ipoin) !- coord(idime,ipoin)

                 end do
! print*, "DEBUG: bvess_ale = ", bvess_ale(:,ipoin)
!               end if

           end do

        end if

     end if

  end if

  if( associated(xvalu) ) deallocate( xvalu )
  if( associated(svalu) ) deallocate( svalu )
  if( associated(force) ) deallocate( force )

end subroutine nsa_fsiexch
!>
