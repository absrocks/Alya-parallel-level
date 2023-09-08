subroutine got_updunk(itask)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_updunk
  ! NAME 
  !    got_updunk
  ! DESCRIPTION
  !    This routine performs several types of updates for the droplet velocity and
  !    droplet concentration.
  ! USED BY
  !    got_begste (itask=1)
  !    got_begite (itask=2)
  !    got_endite (itask=3, inner loop) 
  !    got_endite (itask=4, outer loop) 
  !    got_endste (itask=5)
  !    got_endste (itask=6)               Updates Droplet concentration
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_gotita
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,itotv,idime,icomp,ielem,igaus
  integer(ip)             :: pelty,pgaus,ndofn,ibopo
  real(rp)                :: rela1_got,worve(3),worma(3),udotn

  if(kfl_paral/=0) then

     select case (itask)

     case(1) 
        !
        ! Assign u(n,0,*) <-- u(n-1,*,*), initial guess for outer iterations
        !
        icomp=min(3,ncomp_got) 
        if(kfl_probl_got<=2) then
           do ipoin=1,npoin
              do idime=1,ndime
                 vdrop(idime,ipoin,2) = vdrop(idime,ipoin,icomp)
              end do
           end do
        end if
        if(kfl_probl_got/=2) then
           do ipoin=1,npoin
              cdrop(ipoin,2) = cdrop(ipoin,icomp)
           end do
        end if

     case(2)
        !
        ! Assign u(n,i,0) <-- u(n,i-1,*), initial guess for inner iterations
        !
        if(kfl_probl_got<=2) then
           do ipoin=1,npoin
              do idime=1,ndime
                 vdrop(idime,ipoin,1) = vdrop(idime,ipoin,2)
              end do
           end do
        end if
        if(kfl_probl_got/=2) then          
           do ipoin=1,npoin
              cdrop(ipoin,1)=cdrop(ipoin,2)
           end do
        end if
        if(kfl_algor_got==4) then
           !
           ! Block Gauss-Seidel: UNKNO=[ u1 v1 u2 v2 ... p1 p2 ... ]
           !
           ndofn=npoin*ndime
           do ipoin=1,npoin
              itotv=(ipoin-1)*ndime
              do idime=1,ndime
                 itotv=itotv+1
                 unkno(itotv)=vdrop(idime,ipoin,2)
              end do
              unkno(ndofn+ipoin)=cdrop(ipoin,2)
           end do
        else
           !
           ! Monolithic: UNKNO=[ u1 v1 p1 u2 v2 p2 ... ]
           !
           if(kfl_probl_got==1) then
              do ipoin=1,npoin
                 itotv=(ipoin-1)*ndofn_got(3)
                 do idime=1,ndime
                    itotv=itotv+1
                    unkno(itotv)=vdrop(idime,ipoin,2)
                 end do
                 unkno(itotv+1)=cdrop(ipoin,2)
              end do
           else if(kfl_probl_got==2) then
              do ipoin=1,npoin
                 itotv=(ipoin-1)*ndime
                 do idime=1,ndime
                    itotv=itotv+1
                    unkno(itotv)=vdrop(idime,ipoin,2)
                 end do
              end do
           else if(kfl_probl_got==3) then
              do ipoin=1,npoin
                 unkno(ipoin)=cdrop(ipoin,2)
              end do              
           end if
        end if

     case(3)
        !
        ! Assign u(n,i,j-1) <-- u(n,i,j), update of the droplet velocity and droplet concentration
        !
        if(postp(1)%npp_stepi(9)>0) then
           do ipoin=1,npoin
              cdold_got(ipoin)=cdrop(ipoin,1)
           end do
        end if
        if(postp(1)%npp_stepi(10)>0) then
           do ipoin=1,npoin
              do idime=1,ndime
                 vdold_got(idime,ipoin)=vdrop(idime,ipoin,1)
              end do
           end do
        end if
        !
        ! Put normal component of Droplet velocity when u.n<0
        !
        goto 10
        do ipoin=1,npoin
           ibopo=lpoty(ipoin)
           if(ibopo/=0) then
              if(  .not.(coord(1,ipoin)<xmini_got(1).or.coord(1,ipoin)>xmaxi_got(1).or.&
                   &     coord(2,ipoin)<xmini_got(2).or.coord(2,ipoin)>xmaxi_got(2))) then
                 udotn=0.0_rp
                 do idime=1,ndime
                    udotn=udotn+vdrop(idime,ipoin,1)*exnor(idime,1,ibopo)
                 end do
                 if(udotn<0.0_rp) then
                    itotv=(ipoin-1)*ndofn_got(3)
                    do idime=1,ndime 
                       itotv=itotv+1
                       worve(idime)=unkno(itotv)
                    end do
                    call mbvatb(worma,exnor(1,1,ibopo),worve,ndime,ndime) ! Global to local
                    worma(1)=0.0_rp                                       ! u.n=0
                    call mbvab0(worve,exnor(1,1,ibopo),worma,ndime,ndime) ! Local to global
                    itotv=(ipoin-1)*ndofn_got(3)
                    do idime=1,ndime 
                       itotv=itotv+1
                       unkno(itotv)=worve(idime)
                    end do
                 end if
              end if
           end if
        end do
10      continue

        rela1_got=1.0_rp-relax_got
        if(kfl_algor_got==4) then
           !
           ! Block Gauss-Seidel: UNKNO=[ u1 v1 u2 v2 ... p1 p2 ... ]
           !
           ndofn=npoin*ndime
           do ipoin=1,npoin
              itotv=(ipoin-1)*ndime
              do idime=1,ndime
                 itotv=itotv+1
                 vdrop(idime,ipoin,1) = relax_got*unkno(itotv)&
                      &               + rela1_got*vdrop(idime,ipoin,1)
              end do
              cdrop(ipoin,1)=relax_got*unkno(ndofn+ipoin)&
                   +rela1_got*cdrop(ipoin,1)
              cdrop(ipoin,1)=max(cdrop(ipoin,1),almax_got*cutof_got)              
           end do
        else
           !
           ! Monolithic: UNKNO=[ u1 v1 p1 u2 v2 p2 ... ]
           !
           if(kfl_probl_got==1) then 
              do ipoin=1,npoin
                 itotv=(ipoin-1)*ndofn_got(3)
                 do idime=1,ndime
                    itotv=itotv+1
                    vdrop(idime,ipoin,1)=relax_got*unkno(itotv)&
                         +rela1_got*vdrop(idime,ipoin,1)
                 end do
                 cdrop(ipoin,1)=relax_got*unkno(itotv+1)&
                      +rela1_got*cdrop(ipoin,1)
                 cdrop(ipoin,1)=max(cdrop(ipoin,1),almax_got*cutof_got)
              end do
           else if(kfl_probl_got==2) then 
              do ipoin=1,npoin
                 itotv=(ipoin-1)*ndime
                 do idime=1,ndime
                    itotv=itotv+1
                    vdrop(idime,ipoin,1)=relax_got*unkno(itotv)&
                         +rela1_got*vdrop(idime,ipoin,1)
                 end do
              end do
           else if(kfl_probl_got==3) then 
              do ipoin=1,npoin
                 cdrop(ipoin,1)=relax_got*unkno(ipoin)&
                      +rela1_got*cdrop(ipoin,1)
                 cdrop(ipoin,1)=max(cdrop(ipoin,1),almax_got*cutof_got)
              end do
           end if
        end if

     case(4)
        !
        ! Assign u(n,i-1,*) <-- u(n,i,*)
        !        
        ! If droplet velocity residual is required, save previous interation in veold_got
        !
        if(kfl_probl_got==1) then
           do ipoin=1,npoin
              do idime=1,ndime
                 vdrop(idime,ipoin,2)=vdrop(idime,ipoin,1)
              end do
              cdrop(ipoin,2)=cdrop(ipoin,1)
           end do
        else if(kfl_probl_got==2) then
           do ipoin=1,npoin
              do idime=1,ndime
                 vdrop(idime,ipoin,2)=vdrop(idime,ipoin,1)
              end do
           end do
        else if(kfl_probl_got==3) then
           do ipoin=1,npoin
              cdrop(ipoin,2)=cdrop(ipoin,1)
           end do
        end if

     case(5)
        !
        ! u(n-1,*,*) <-- u(n,*,*)
        !        
	if(kfl_tiacc_got==2) then
           !
           ! Crank-Nicolson method 
           !
           do ipoin=1,npoin
              do idime=1,ndime
                 vdrop(idime,ipoin,1) = 2.0_rp*vdrop(idime,ipoin,1)&
                      -vdrop(idime,ipoin,3)
              end do
              cdrop(ipoin,1) = 2.0_rp*cdrop(ipoin,1)-cdrop(ipoin,3)
           end do
        end if
        if(kfl_probl_got==1) then
           do ipoin=1,npoin
              do idime=1,ndime
                 vdrop(idime,ipoin,3) = vdrop(idime,ipoin,1)
              end do
              cdrop(ipoin,3) = cdrop(ipoin,1)
           end do
        else if(kfl_probl_got==2) then
           do ipoin=1,npoin
              do idime=1,ndime
                 vdrop(idime,ipoin,3) = vdrop(idime,ipoin,1)
              end do
           end do
        else if(kfl_probl_got==3) then
           do ipoin=1,npoin
              cdrop(ipoin,3) = cdrop(ipoin,1)
           end do
        end if

        if(kfl_sgsti_got==1) then
           !
           ! Time tracking of the subscales
           !
           if(kfl_tiacc_got==2) then
              do ielem=1,nelem
                 pelty=ltype(ielem)
                 pgaus=ngaus(pelty)
                 do igaus=1,pgaus
                    do idime=1,ndime
                       vesgs(ielem)%a(idime,igaus,1) = &
                            2.0_rp*vesgs(ielem)%a(idime,igaus,1)&
                            -vesgs(ielem)%a(idime,igaus,2)
                    end do
                 end do
              end do
           end if
           do ielem=1,nelem
              pelty=ltype(ielem)
              pgaus=ngaus(pelty)
              do igaus=1,pgaus
                 do idime=1,ndime
                    vesgs(ielem)%a(idime,igaus,2) = vesgs(ielem)%a(idime,igaus,1)
                 end do
              end do
           end do
        end if

     case(11)
        !
        ! Assign u(n,i,*)  <-- u(n-1,*,*), initial guess after reading restart
        !        u'(n,i,*) <-- u'(n-1,*,*)
        !
        icomp=min(3,ncomp_got)
        if(kfl_probl_got==1) then 
           do ipoin=1,npoin
              do idime=1,ndime           
                 vdrop(idime,ipoin,1) = vdrop(idime,ipoin,icomp) 
              end do
              cdrop(ipoin,1) = cdrop(ipoin,icomp) 
           end do
        else if(kfl_probl_got==2) then 
           do ipoin=1,npoin
              do idime=1,ndime           
                 vdrop(idime,ipoin,1) = vdrop(idime,ipoin,icomp) 
              end do
           end do
        else if(kfl_probl_got==3) then 
           do ipoin=1,npoin
              cdrop(ipoin,1) = cdrop(ipoin,icomp) 
           end do
        end if
        if(kfl_sgsti_got==1) then
           do ielem=1,nelem
              pelty=ltype(ielem)
              pgaus=ngaus(pelty)
              do igaus=1,pgaus
                 do idime=1,ndime
                    vesgs(ielem)%a(idime,igaus,1)=vesgs(ielem)%a(idime,igaus,2)
                 end do
              end do
           end do
        end if

     end select

  end if

end subroutine got_updunk

