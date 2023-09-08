subroutine ale_updbcs()
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_updbcs
  ! NAME 
  !    ale_updbcs
  ! DESCRIPTION
  !    This routine sets up the initial condition for the mesh velocity
  ! USED BY
  !    ale_begste
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_alefor
  use def_kermod
  use mod_ker_space_time_function

  implicit none
  integer(ip) :: ipoin,ifunc,kk,idime
  real   (rp) :: vefun(3), funtime(3),funold,xx


  if( kfl_conbc_ale == 0 .and. INOTMASTER ) then
     !
     ! Space/Time functions
     !
     if( number_space_time_function > 0 ) then
        do ipoin = 1,npoin
           if( kfl_funno_ale(ipoin) < 0 ) then 
              ifunc = -kfl_funno_ale(ipoin)            
              call ker_space_time_function(&
                 ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,vefun(1))
              call ker_space_time_function(&
                 ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim-dtime,vefun(2))
              bvess_ale(1:ndime,ipoin) = bvess_ref(1:ndime,ipoin)*(vefun(1) - vefun(2))
           end if
        end do
     end if

     if( kexist_tran_fiel > 0 ) then
     !
     ! Transient fields
     !
        do ipoin = 1,npoin
           if( kfl_fixno_ale(1,ipoin) /= 0 ) then
              if( kfl_funno_ale(ipoin) /= 0 ) then 
                 ifunc = kfl_funno_ale(ipoin)
                 if (ifunc  > interval_funno ) then 
                    !
                    ! BC is associated to a displacement field, meaning that it is a transient one
                    !
                    ifunc = ifunc - interval_funno  ! correcting the real ifunc value, which is the associated field number
                    kk = k_tran_fiel(ifunc)
                    xx = x_tran_fiel(ifunc)
                    do idime = 1,ndime
                       funtime(idime) = xfiel(ifunc) % a(idime,ipoin,kk) * xx + &
                            xfiel(ifunc) % a(idime,ipoin,kk+1) * (1.0_rp-xx)
                    end do

!!$           if (ipoin .eq. 309) then
!!$              write(6666,100) cutim,xfiel(ifunc) % a(1,ipoin,kk),xfiel(ifunc) % a(1,ipoin,kk+1),funtime(1)
!!$           end if
!!$
!!$           if (ipoin .eq. 12) then
!!$              write(6667,100) cutim,xfiel(ifunc) % a(1,ipoin,kk),xfiel(ifunc) % a(1,ipoin,kk+1),funtime(1)
!!$           end if
                    ! done only for 1st order movement schemes                    
                    do idime = 1,ndime
                       bvess_ale(idime,ipoin) = funtime(idime)
                    end do

                 else
                    call ale_funcre(&
                         funpa_ale(ifunc) % a,   &
                         kfl_funty_ale(ifunc,2), &
                         kfl_funty_ale(ifunc,1), &
                         cutim,                  &
                         coord(1,ipoin),         &
                         coord_ale(1,ipoin,2),   &
                         bvess_ale(1,ipoin)      )
                 end if

              end if
           end if
        end do

     end if

  end if

end subroutine ale_updbcs
