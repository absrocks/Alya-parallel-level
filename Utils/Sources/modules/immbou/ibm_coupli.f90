subroutine ibm_coupli(itask)
  !-----------------------------------------------------------------------
  !****f* immbou/ibm_coupli
  ! NAME
  !    ibm_coupli
  ! DESCRIPTION
  !    This routines couples iimbou
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_domain
  use def_immbou
  use mod_kdtree

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,idime,kpoin,iimbo
  integer(ip)             :: dummi
  real(rp)                :: dummr,dumma(ndime),propo(ndime)
  real(rp)                :: x(3),v(3),coor1(ndime)

  select case ( itask )

  case ( ITASK_ENDITE )

     !-------------------------------------------------------------------
     !
     ! ENDITE
     !
     !-------------------------------------------------------------------

     if( kfl_coupl(ID_ALEFOR,ID_IMMBOU) == 1 ) then
        !
        ! Coupling with Alefor
        ! For fixed mesh ALE, fixity is negative so that the mesh velocity is inverted
        !
        if( INOTMASTER ) then
           if( associated(kfl_fixno_ale) ) then

              do ipoin = 1,npoin
                 if (lpoty(ipoin) /= 0) then
                    do idime = 1,ndime
                       kfl_fixno_ale(idime,ipoin) = 3
                       bvess_ale(idime,ipoin)     = 0.0_rp
                    end do
                 else
                    do idime = 1,ndime
                       kfl_fixno_ale(idime,ipoin) = -1
                    end do
                 end if
              end do

              !do iimbo = 1,nimbo
              !   call kdtree(&
              !        0_ip,mnoib,imbou(iimbo) % npoib,imbou(iimbo) % nboib,           &
              !        imbou(iimbo) % cooi2,imbou(iimbo) % lnoib,imbou(iimbo) % ltyib, &
              !        imbou(iimbo) % fabox,imbou(iimbo) % bobox,imbou(iimbo) % sabox, &
              !        imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist, & 
              !        imbou(iimbo) % lnele )                            
              !end do
              do ipoin = 1,npoin
                 !
                 ! Capture the mesh velocity
                 !
                 if (lntib(ipoin) < 0_ip) then
                    iimbo    = abs(lntib(ipoin))
                    !if (lnint(ipoin) % limit > 1) then
                    !   call dpopar(&
                    !        1_ip,coord(1:ndime,ipoin),&
                    !        imbou(iimbo) % npoib,mnoib,imbou(iimbo) % nboib,1.0e10_rp,&
                    !        imbou(iimbo) % ltyib,imbou(iimbo) % lnoib,imbou(iimbo) % cooi2,&
                    !        dummr,dumma,propo,dummi,&
                    !        imbou(iimbo) % fabox,imbou(iimbo) % sabox,imbou(iimbo) % blink,&
                    !        imbou(iimbo) % struc,imbou(iimbo) % ldist,imbou(iimbo) % lnele)
                    !elseif (lnint(ipoin) %limit == 1) then
                    !   do idime = 1,ndime
                    !      coor1(idime) =  coord(idime,ipoin) - 2.0_rp*(coord(idime,lnint(ipoin)%lnode(1)) - coord(idime,ipoin))
                    !   end do
                    !   call faceli(&
                    !        imbou(iimbo) % sabox,imbou(iimbo) % blink,imbou(iimbo) %ltyib,imbou(iimbo) %lnoib,imbou(iimbo) %cooi2, & 
                    !        coor1,coord(1,lnint(ipoin)%lnode(1)),coord(1,ipoin),mnoib,imbou(iimbo) %nboib,propo,dummi,dummr)              
                    !end if
                    !x(1)     = propo(1) - imbou(iimbo) % posil(    1,1)
                    !x(2)     = propo(2) - imbou(iimbo) % posil(    2,1)
                    !x(3)     = 0.0_rp
                    !if ( ndime == 3)  x(ndime) = propo(ndime) - imbou(iimbo) % posil(ndime,1)
                    !v(1)     = 0.0_rp
                    !v(2)     = 0.0_rp
                    !v(3)     = 0.0_rp
                    !call vecpro(imbou(iimbo)%veloa,x,v,3_ip)            
                    !v(1)     = v(1)     + imbou(iimbo) % velol(    1,1)
                    !v(2)     = v(2)     + imbou(iimbo) % velol(    2,1)
                    !if ( ndime == 3) v(3) = v(3) + imbou(iimbo) % velol(3,1)

                    !v(1) =  sin(pi*propo(1)-0.7_rp)*sin(pi*propo(2)+0.2_rp)
                    !v(2) =  cos(pi*propo(1)-0.7_rp)*cos(pi*propo(2)+0.2_rp)
                    !v(1) =  propo(1)
                    !v(2) = -propo(2)
                    do idime = 1,ndime               
                       kfl_fixno_ale(idime,ipoin)  = 3
                       !bvess_ale(idime,ipoin)      = propo(idime) - coord(idime,ipoin)
                       !veloc(idime,ipoin,3)        = v(idime) 
                       bvess_ale(idime,ipoin)      = -(imbou(iimbo) % posil(idime,1) - imbou(iimbo) % posil(idime,2))
                    end do

                 end if
              end do
              !do iimbo = 1,nimbo
              !   call kdtree(&
              !        0_ip,mnoib,imbou(iimbo) % npoib,imbou(iimbo) % nboib,           &
              !        imbou(iimbo) % cooib,imbou(iimbo) % lnoib,imbou(iimbo) % ltyib, &
              !        imbou(iimbo) % fabox,imbou(iimbo) % bobox,imbou(iimbo) % sabox, &
              !        imbou(iimbo) % blink,imbou(iimbo) % struc,imbou(iimbo) % ldist, & 
              !        imbou(iimbo) % lnele )                                     
              !end do
           end if
        end if

     else if( kfl_coupl(ID_ALEFOR,ID_IMMBOU) == 2 ) then
        !
        ! Coupling with Alefor
        ! Pure ALE
        !
        if( INOTMASTER ) then
           if( associated(kfl_fixno_ale) ) then
              if ( kfl_mvext_ibm == 1 ) then ! set x & y components of displacement = disp COG
                 do ipoin = 1,npoin
                    do idime = 1,2
                       if (kfl_fixno_ale(idime,ipoin) == 1) bvess_ale(idime,ipoin) = imbou(1) % posil(idime,1) & 
                            - imbou(1) % posil(idime,2)    !delta x_cog
                    end do
                 end do
              end if
              do iimbo = 1,nimbo
                 if( imbou(iimbo) % kfl_typeb >= 1 ) then
                    do kpoin = 1,imbou(iimbo) % npoib
                       ipoin = imbou(iimbo) % lninv(kpoin)
                       do idime = 1,ndime    
                          bvess_ale(idime,ipoin)     = imbou(iimbo) % cooib(idime,kpoin) - coord(idime,ipoin)
                          kfl_fixno_ale(idime,ipoin) = 1
                       end do
                    end do
                 end if
              end do
           end if
        end if

     end if

     if( kfl_coupl(ID_NASTIN,ID_IMMBOU) /= 0 ) then
        !
        ! Coupling with Nastin: Cut elements defined in ibm_insout
        !
        if( INOTMASTER ) then
           do iimbo = 1,nimbo
              if( imbou(iimbo) % kfl_typeb == -2 ) then
                 continue
              end if
           end do
        end if
     end if

  end select

end subroutine ibm_coupli
