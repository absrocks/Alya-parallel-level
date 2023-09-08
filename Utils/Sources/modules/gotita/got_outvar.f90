subroutine got_outvar(ivari)
  !------------------------------------------------------------------------
  !****f* Gotita/got_output
  ! NAME 
  !    got_output
  ! DESCRIPTION
  !    Output a postprocess variable
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    got_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_gotita
  use mod_postpr
  use mod_gradie
  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip)             :: pelty,porde
  integer(ip)             :: ibopo,ipoin,idime,jdime,pnode,ielem,inode
  real(rp)                :: vdiff(3),vnorm,ReD,CDReD,vegot(3),xter2
  real(rp)                :: xterm,dtcri,gpvel(3),gpvdr(3),gpcdr,rutim
  real(rp)                :: gppor,gpdif,gpvno,chale(2)
  character(5)            :: wopos_tmp(2)
  !
  ! Define postprocess variable
  !

  select case (ivari)  

  case(0)

     return

  case(1)
     !
     ! Droplet concentration
     !
     gesca => cdrop(:,1) 

  case(2)
     !
     ! Droplet velocity
     !
     if(pos_cutof_got==0.0_rp) then
        gevec => vdrop(:,:,1) 
     else
        if(kfl_paral/=0) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              if(cdrop(ipoin,1)<pos_cutof_got*almax_got) then
                 do idime=1,ndime
                    gevec(idime,ipoin)=0.0_rp
                 end do
              else
                 do idime=1,ndime
                    gevec(idime,ipoin)=vdrop(idime,ipoin,1)
                 end do
              end if
           end do
        end if
     end if

  case(3)
     !
     ! Droplet Momentum
     !
     if(kfl_paral/=0) then
        call memgen(zero,ndime,npoin)
        do ipoin=1,npoin
           do idime=1,ndime
              gevec(idime,ipoin)=cdrop(ipoin,1)*vdrop(idime,ipoin,1)
           end do
        end do
     end if

  case(6)
     !
     ! Drag
     !
     if(kfl_paral/=0) then
        call memgen(zero,npoin,0_ip)
        do ipoin=1,npoin
           if(kfl_velfu_got==0) then
              do idime=1,ndime
                 vdiff(idime)=veloc(idime,ipoin,1)-vdrop(idime,ipoin,1)
              end do
           else if(kfl_velfu_got==-1) then
              do idime=1,ndime
                 vdiff(idime)=veloc_got(idime,ipoin)-vdrop(idime,ipoin,1)
              end do
           else
              call got_veloci(ipoin,coord(1,ipoin),vegot)
              do idime=1,ndime
                 vdiff(idime)=vegot(idime)-vdrop(idime,ipoin,1)
              end do
           end if
           call vecnor(vdiff,ndime,vnorm,2_ip)
           Red=deair_got*ddrop_got*veair_got*vnorm/muair_got
           call got_dragco(ReD,CDReD)
           gesca(ipoin)=CDReD/(24.0_rp*kfact_got)
        end do
     end if

  case(7)
     !
     ! Air velocity
     !
     if(kfl_velfu_got==0) then
        gevec => veloc(:,:,1) 
     else if(kfl_velfu_got==-1) then
        gevec => veloc_got
     else
        if(kfl_paral/=0) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              call got_veloci(ipoin,coord(1,ipoin),gevec(1,ipoin))
           end do
        end if
     end if

  case(8)
     !
     ! Terms of the equation
     !
     wopos_tmp(2)=postp(1)%wopos(2,ivari)
     if(kfl_paral/=0) then
        call memgen(zero,-ndime,npoin)
        call memgen(zero, npoin, 0_ip)
        call memgen(zero, ndime,npoin)
     end if
     !
     ! TERM1: |alpha*du/(theta*dt)|
     !

     !
     ! TERM2: |alpha*(u.grad)u|
     !     
     call gradie(vdrop(:,:,1),geten)
     do ipoin=1,npoin
        gesca(ipoin)=0.0_rp
        do idime=1,ndime
           xterm=0.0_rp
           do jdime=1,ndime
              xterm=xterm+vdrop(jdime,ipoin,1)*geten(jdime,idime,ipoin)
           end do
           xterm=cdrop(ipoin,1)*xterm
           gesca(ipoin)=gesca(ipoin)+xterm*xterm
        end do
        gesca(ipoin)=sqrt(gesca(ipoin))
     end do
     wopos_tmp(1)=postp(1)%wopos(1,ivari)(1:4)//'2'
     call postpr(gesca,wopos_tmp,ittim,cutim)
     !
     ! TERM3: |alpha*sig*u|
     !
     do ipoin=1,npoin
        gesca(ipoin)=0.0_rp
        if(kfl_velfu_got==0) then
           do idime=1,ndime
              vegot(idime)=veloc(idime,ipoin,1)
           end do
        else if(kfl_velfu_got==-1) then
           do idime=1,ndime
              vegot(idime)=veloc_got(idime,ipoin)
           end do
        else
           call got_veloci(ipoin,coord(1,ipoin),vegot)
        end if
        do idime=1,ndime
           vdiff(idime)=vegot(idime)-vdrop(idime,ipoin,1)
        end do
        call vecnor(vdiff,ndime,vnorm,2_ip)
        Red=deair_got*ddrop_got*veair_got*vnorm/muair_got
        call got_dragco(ReD,CDReD)
        xterm=CDReD/(24.0_rp*kfact_got)
        do idime=1,ndime
           gesca(ipoin)=gesca(ipoin)&
                +(xterm*vdrop(idime,ipoin,1)*cdrop(ipoin,1))**2
        end do
        gesca(ipoin)=sqrt(gesca(ipoin))
     end do
     wopos_tmp(1)=postp(1)%wopos(1,ivari)(1:4)//'3'
     call postpr(gesca,wopos_tmp,ittim,cutim)
     !
     ! TERM4: |alpha*sig*ua|
     !
     do ipoin=1,npoin
        gesca(ipoin)=0.0_rp
        if(kfl_velfu_got==0) then
           do idime=1,ndime
              vegot(idime)=veloc(idime,ipoin,1)
           end do
        else if(kfl_velfu_got==-1) then
           do idime=1,ndime
              vegot(idime)=veloc_got(idime,ipoin)
           end do
        else
           call got_veloci(ipoin,coord(1,ipoin),vegot)
        end if
        do idime=1,ndime
           vdiff(idime)=vegot(idime)-vdrop(idime,ipoin,1)
        end do
        call vecnor(vdiff,ndime,vnorm,2_ip)
        Red=deair_got*ddrop_got*veair_got*vnorm/muair_got
        call got_dragco(ReD,CDReD)
        xterm=CDReD/(24.0_rp*kfact_got)
        do idime=1,ndime
           gesca(ipoin)=gesca(ipoin)&
                +(xterm*vegot(idime)*cdrop(ipoin,1))**2
        end do
        gesca(ipoin)=sqrt(gesca(ipoin))
     end do
     wopos_tmp(1)=postp(1)%wopos(1,ivari)(1:4)//'4'
     call postpr(gesca,wopos_tmp,ittim,cutim)
     !
     ! TERM6: |alpha*div(u)|
     !
     do ipoin=1,npoin
        xterm=0.0_rp
        do idime=1,ndime
           xterm=xterm+geten(idime,idime,ipoin)
        end do
        gesca(ipoin)=abs(xterm*cdrop(ipoin,1))
     end do
     wopos_tmp(1)=postp(1)%wopos(1,ivari)(1:4)//'6'
     call postpr(gesca,wopos_tmp,ittim,cutim)
     !
     ! TERM7: |grad(alpha).u|
     !
     call gradie(cdrop(:,1),gevec)
     do ipoin=1,npoin
        xterm=0.0_rp
        do idime=1,ndime
           xterm=xterm+gevec(idime,ipoin)*vdrop(idime,ipoin,1)
        end do
        gesca(ipoin)=abs(xterm)
     end do
     wopos_tmp(1)=postp(1)%wopos(1,ivari)(1:4)//'7'
     call postpr(gesca,wopos_tmp,ittim,cutim)

  case(9)
     !
     ! Water volume fraction Residual
     !
     call memgen(zero, npoin, 0_ip)
     do ipoin=1,npoin
        xterm=abs(cdrop(ipoin,1))
        if(xterm==0.0_rp) xterm=1.0_rp
        gesca(ipoin)=abs(cdrop(ipoin,1)-cdold_got(ipoin))/xterm
     end do

  case(10)
     !
     ! Droplet velocity residual
     !
     call memgen(zero, npoin, 0_ip)
     do ipoin=1,npoin
        xterm=0.0_rp
        xter2=0.0_rp
        do idime=1,ndime
           xterm=xterm+(vdrop(idime,ipoin,1)-vdold_got(idime,ipoin))**2
           xter2=xter2+(vdrop(idime,ipoin,1))**2
        end do
        if(xter2==0.0_rp) xter2=1.0_rp
        gesca(ipoin)=sqrt(xterm/xter2)
     end do

  case(11)
     !
     ! Fixity
     !
     call memgen(zero,npoin,0_ip)
     do ipoin=1,npoin
        gesca(ipoin)=real(kfl_fixno_got(ipoin))
     end do

  case(12)
     !
     ! Time step
     !
     if(associated(vmasc)) then
        call memgen(zero, 2_ip,npoin)
        call memgen(zero,npoin, 0_ip)
        do ielem=1,nelem
           pelty=ltype(ielem)
           pnode=nnode(pelty)
           porde=lorde(pelty)
           call got_elmgat(&
                2_ip,1_ip,pnode,lnods(1,ielem),elvdr_got,&
                elcdr_got,elcod_got,elvel_got,eldif_got)
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              call got_elmtss(&
                   pnode,porde,elcod_got,elvel_got,elvdr_got,elcdr_got,&
                   eldif_got,elmar(pelty)%shacg,elmar(pelty)%dercg,&
                   elmar(pelty)%weicg,dtcri)
              gevec(1,ipoin)=gevec(1,ipoin)+vmasc(ipoin)*dtcri
              gevec(2,ipoin)=gevec(2,ipoin)+vmasc(ipoin)
           end do
        end do
        do ipoin=1,npoin
           gesca(ipoin)=safet_got*gevec(1,ipoin)/gevec(2,ipoin)
        end do
     end if

  case(13)
     !
     ! Diffusion
     !
     if(kfl_diffu_got==1.and.associated(vmasc)) then
        call memgen(zero, 2_ip,npoin)
        call memgen(zero,npoin, 0_ip)
        do ielem=1,nelem
           pelty=ltype(ielem)
           pnode=nnode(pelty)
           call got_elmgat(&
                2_ip,1_ip,pnode,lnods(1,ielem),elvdr_got,&
                elcdr_got,elcod_got,elvel_got,eldif_got)
           call got_elmave(&
                pnode,lnods(1,ielem),elcdr_got,elvdr_got,elvel_got,&
                gpcdr,gpvdr,gpvel)
           call vecnor(gpvdr,ndime,gpvno,2_ip)
           call got_elmpro(&
                pnode,1_ip,1_ip,1_ip,eldif_got,elmar(pelty)%shape,&
                gpvdr,gpvel,gpvno,gpcdr,chale,gppor,gpdif)
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              gevec(1,ipoin)=gevec(1,ipoin)+vmasc(ipoin)*gpdif
              gevec(2,ipoin)=gevec(2,ipoin)+vmasc(ipoin)
           end do
        end do
        do ipoin=1,npoin
           gesca(ipoin)=gevec(1,ipoin)/gevec(2,ipoin)
        end do
     else if(kfl_diffu_got==2) then
        gesca => diffm_got
     end if

  case(14)
     !
     ! Collection efficiency beta
     !
     !ibopo =  1
     !call memgen(zero,nbopo,zero)
     !do ipoin=1,npoin
     !   ibopo=lpoty(ipoin)
     !   if(ibopo/=0) then
     !      if(kfl_fixno_got(ipoin)==0) then
     !         xterm=0.0_rp
     !         do idime=1,ndime
     !            xterm=xterm+vdrop(idime,ipoin,1)*exnor(idime,1,ibopo)
     !         end do
     !         if(xterm<0.0_rp) then
     !            gesca(ibopo)=cdrop(ipoin,1)*xterm
     !         else
     !            gesca(ibopo)=0.0_rp
     !         end if
     !      else
     !         gesca(ibopo)=0.0_rp
     !      end if
     !   end if
     !end do

     if(kfl_paral/=0) then
        call memgen(zero,ndime,npoin)
        do ipoin=1,npoin
           if(  .not.(coord(1,ipoin)<xmini_got(1).or.coord(1,ipoin)>xmaxi_got(1).or.&
                &     coord(2,ipoin)<xmini_got(2).or.coord(2,ipoin)>xmaxi_got(2).or.&
                &     coord(ndime,ipoin)<xmini_got(ndime).or.coord(ndime,ipoin)>xmaxi_got(ndime))) then
              do idime=1,ndime
                 gevec(idime,ipoin)=0.0_rp 
              end do
              ibopo=lpoty(ipoin)
              if(ibopo/=0.and.kfl_fixno_got(ipoin)==0) then
                 xterm=0.0_rp
                 do idime=1,ndime
                    xterm=xterm+vdrop(idime,ipoin,1)*exnor(idime,1,ibopo)
                 end do
                 if(xterm>0.0_rp) then
                    do idime=1,ndime
                       gevec(idime,ipoin)=-cdrop(ipoin,1)*xterm*exnor(idime,1,ibopo)
                    end do
                 end if
              end if
           end if
        end do
     end if

  end select

  call outvar(&
       ivari,&
       ittim,rutim,postp(1) % wopos(1,ivari))

end subroutine got_outvar

subroutine got_elmave(&
     pnode,lnods,elcdr,elvdr,elvel,&
     gpcdr,gpvdr,gpvel)
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  implicit none
  integer(ip), intent(in)  :: pnode,lnods(pnode)
  real(rp),    intent(in)  :: elcdr(pnode),elvdr(ndime,pnode)
  real(rp),    intent(in)  :: elvel(ndime,pnode)
  real(rp),    intent(out) :: gpcdr,gpvdr(ndime),gpvel(ndime)
  integer(ip)              :: inode,idime,ipoin
  real(rp)                 :: rnode
  
  gpcdr=0.0_rp
  do idime=1,ndime
     gpvdr(idime)=0.0_rp
     gpvel(idime)=0.0_rp
  end do
  
  do inode=1,pnode
     ipoin=lnods(inode)
     gpcdr=gpcdr+elcdr(inode)
     do idime=1,ndime
        gpvdr(idime)=gpvdr(idime)+elvdr(idime,inode)
        gpvel(idime)=gpvel(idime)+elvel(idime,inode)
     end do
  end do
  
  rnode=1.0_rp/real(pnode)
  gpcdr=gpcdr*rnode
  do idime=1,ndime
     gpvdr(idime)=gpvdr(idime)*rnode
     gpvel(idime)=gpvel(idime)*rnode
  end do

end subroutine got_elmave
