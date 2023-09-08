 subroutine nsi_inipre()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_inipre
  ! NAME
  !    nsi_inipre
  ! DESCRIPTION
  !    This routine solves the initial pressure
  ! USES
  !    nsi_ifconf
  !    nsi_solmon
  !    nsi_solbgs
  !    nsi_rotunk
  ! USED BY
  !    nsi_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use def_solver, only : SOL_MATRIX_HAS_CHANGED
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: izmat,izdom,ipoin,iboun,inodb,ibopo
  integer(ip) :: jbopo,izdod,jpoin,jzdom
  real(rp)    :: dtinv_nsi_save
  real(rp)    :: Qd

  call livinf(81_ip,'INITIAL PRESSURE USING LAPLACIAN',0_ip)

  ivari_nsi =  ivari_nsi_cont
  solve_sol => solve(2:)
  solve_sol %  kfl_assem = SOL_MATRIX_HAS_CHANGED

  if( INOTMASTER ) then
     !
     ! Initialize matrix
     !
     do izmat = 1,solve_sol(1) % nzmat
        lapla_nsi(izmat) = 0.0_rp
     end do
     do ipoin = 1,npoin
        rhsid(ipoin) = 0.0_rp
        unkno(ipoin) = 0.0_rp
     end do
     !
     ! Boundary values
     !
     call memgen(zero,nbopo,zero)

     do iboun = 1,nboun
        if( kfl_fixbo_nsi(iboun) == 2 .or. kfl_fixbo_nsi(iboun) == 6 .or. kfl_fixbo_nsi(iboun) == 20 ) then
           do inodb = 1,nnode(ltypb(iboun))
              ipoin = lnodb(inodb,iboun)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                    gesca(ibopo) = bvnat_nsi(1,iboun,1)
                 end if
              end if
           end do
        end if
     end do
     !
     ! Define some variables
     !
     dtinv_nsi_save = dtinv_nsi
     !
     ! Assemble equation
     !
     if( kfl_vector_nsi == 0 ) then
        !
        ! Classical version
        !
        call nsi_elmope_omp(6_ip)
     else
        !
        ! Vectorized version
        !
        call nsi_elmope_all(6_ip)
     end if
     !
     ! Recover original variables
     !
     dtinv_nsi = dtinv_nsi_save 
     !
     ! Impose b.c.
     !
     if( solve(2)%kfl_symme == 1 ) then

        call runend('NIS_INIPRE: NOT CODED: CHECK IT')
        do ipoin = 1,npoin
           do izdom = r_sym(ipoin),r_sym(ipoin+1) - 2
              jpoin = c_sym(izdom)
              jbopo = lpoty(jpoin)
              if( jbopo /= 0 ) then
                 if( kfl_fixpr_nsi(1,jpoin) > 0 ) then
                    rhsid(ipoin) = rhsid(ipoin) - lapla_nsi(izdom) * gesca(jbopo)
                    lapla_nsi(izdom) = 0.0_rp
                 end if
              end if
           end do
        end do

        do ipoin = 1,npoin

           ibopo = lpoty(ipoin)
           if( ibopo > 0 ) then

              if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                 !
                 ! IZDOD: Diagonal
                 !
                 izdod = r_sym(ipoin+1) - 1
                 Qd = lapla_nsi(izdod)
                 if( abs(Qd) < zeror ) Qd = 1.0_rp
                 !
                 ! Set line to zero
                 !
                 do izdom = r_sym(ipoin),r_sym(ipoin+1) - 1
                    lapla_nsi(izdom) = 0.0_rp
                 end do
                 !
                 ! Prescribe value
                 !
                 lapla_nsi(izdod) = Qd
                 rhsid(ipoin)     = Qd * gesca(ibopo)
                 unkno(ipoin)     = gesca(ibopo)

              end if

           end if

        end do

     else

        do jpoin = 1,npoin
           do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
              ipoin = c_dom(jzdom)
              ibopo = lpoty(ipoin)
              if ( ibopo /= 0 .and. ipoin /= jpoin ) then
                 if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                    rhsid(jpoin) = rhsid(jpoin) - lapla_nsi(jzdom) * gesca(ibopo)
                    lapla_nsi(jzdom) = 0.0_rp
                 end if
              end if
           end do
        end do

        do ipoin = 1,npoin

           ibopo = lpoty(ipoin)
           if( ibopo > 0 ) then

              if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                 !
                 ! IZDOD: Diagonal
                 !
                 izdod = r_dom(ipoin) - 1
                 jpoin = 0
                 do while( jpoin /= ipoin )
                    izdod = izdod + 1
                    jpoin = c_dom(izdod)
                 end do
                 Qd = lapla_nsi(izdod)
                 if( abs(Qd) < zeror ) Qd = 1.0_rp
                 !
                 ! Set line to zero
                 !
                 do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                    lapla_nsi(izdom) = 0.0_rp
                 end do
                 !
                 ! Prescribe value
                 !
                 lapla_nsi(izdod) = Qd
                 rhsid(ipoin)     = Qd * gesca(ibopo)
                 unkno(ipoin)     = gesca(ibopo)

              end if

           end if
        end do

     end if
     call memgen(two,nbopo,zero)
  end if
  !
  ! Solve system
  !
  call solver(rhsid,unkno,lapla_nsi,pmatr)
  !
  ! Tell solver that matrix will be changed
  !
  solve_sol %  kfl_assem = SOL_MATRIX_HAS_CHANGED
  !
  ! Update pressure
  !
  if( INOTMASTER ) then
     do ipoin = 1,npoin
        press(ipoin,1) = unkno(ipoin)
        press(ipoin,2) = unkno(ipoin)
        press(ipoin,3) = unkno(ipoin)
     end do
  end if

end subroutine nsi_inipre
