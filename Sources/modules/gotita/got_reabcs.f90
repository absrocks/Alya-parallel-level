subroutine got_reabcs
  !------------------------------------------------------------------------
  !****f* Gotita/got_reabcs
  ! NAME 
  !    got_reabcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions for GOTITA.
  !  
  !    For conditions on nodes, bvess_got(idime,ipoin,1) contains the
  !    droplet velocity value. 
  !    The different codes for kfl_fixno_got(ipoin) are:
  !    = 1 ... Dirichlet
  !    = 0 ... Free or initial
  !    = 8 ... Interpolated Dirichlet
  ! USES
  !    got_bcntoe
  !    ecoute
  !    memchk
  !    runend
  !    got_bounod
  !    got_autbcs
  ! USED BY
  !    got_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_gotita 
  implicit none
  integer(ip) :: idofn,ipoin,idime,ibopo,inibc
  real(rp)    :: bvmax(4),udotn,veldr(3),vnorm,vangl,alpha

  if(kfl_paral<=0) then

     bvmax(1)=-1.0_rp
     bvmax(2)=-1.0_rp
     bvmax(3)=-1.0_rp
     bvmax(4)=-1.0_rp
     !
     ! Allocate memory for the vectors needed to define the BC's 
     !
     call got_membcs(1_ip)
     !
     ! Reach the nodal-wise section
     !
     call ecoute('got_reabcs')
     do while(words(1)/='BOUND')
        call ecoute('got_reabcs')
     end do
     !
     ! Initial condition
     !
     inibc=1
     if(exists('AIRVE')) then
        inibc=1
     else if(exists('CONST')) then
        inibc=2
     else if(exists('ZEROI')) then
        inibc=3
     else if(exists('AOAIN')) then
        inibc=4
     end if
     !
     ! OJO: QUITAR
     !
     !do ipoin=1,npoin
     !   bvess_got(1,ipoin) = 40.0_rp
     !   bvess_got(2,ipoin) = 0.0_rp
     !   bvess_got(3,ipoin) = 0.0_rp
     !   bvess_got(4,ipoin) = 0.001_rp        
     !   ibopo=lpoty(ipoin)
     !   if(ibopo/=0) then
     !      if(  coord(1,ipoin)<-2.99999_rp.or.&
     !           coord(2,ipoin)<-0.99999_rp.or.&
     !           coord(2,ipoin)> 0.99999_rp.or.&
     !           coord(3,ipoin)< 0.00001_rp.or.&
     !           coord(3,ipoin)> 1.49999_rp) then
     !         kfl_fixno_got(ipoin) = 1
     !      end if
     !   end if
     !end do
     !return
     !
     ! If boundary conditions are computed automatically
     !
     if(exists('AUTOM')) then
        vangl = getrea('ANGLE', 0.0_rp,'#Angle of attack')
        vnorm = getrea('VELOC', 1.0_rp,'#Velocity norm')
        alpha = getrea('ALPHA', 1.0_rp,'#Water volume fraction')
        vangl=vangl/180.0_rp*pi
        veldr(1)=cos(vangl)
        veldr(2)=sin(vangl)
        veldr(3)=0.0_rp
        do ipoin=1,npoin
           ibopo=lpoty(ipoin)
           if(ibopo>0) then
              if(  coord(1,ipoin)<xmini_got(1).or.coord(1,ipoin)>xmaxi_got(1).or.&
                   coord(2,ipoin)<xmini_got(2).or.coord(2,ipoin)>xmaxi_got(2).or.&
                   coord(ndime,ipoin)<xmini_got(ndime).or.coord(ndime,ipoin)>xmaxi_got(ndime)) then
                 udotn=0.0_rp
                 do idime=1,ndime
                    udotn=udotn+veldr(idime)*exnor(idime,1,ibopo)
                 end do
                 if(udotn<=1.0e-12) then
                    kfl_fixno_got(ipoin) = 1
                    if(kfl_probl_got<=2) then
                       do idime=1,ndime
                          bvess_got(idime,ipoin)=vnorm*veldr(idime)
                       end do
                    end if
                    if(kfl_probl_got==1) then
                       bvess_got(ndofn_got(3),ipoin)=alpha
                    end if
                    if(kfl_probl_got==3) then
                       bvess_got(1,ipoin)=alpha
                    end if
                 end if
              end if
           end if
        end do
        almax_got=alpha
        vemax_got=vnorm

        call ecoute('got_reabcs')
        do while(words(1)/='ENDBO')
           if(words(1)=='ONNOD') then
              !
              ! On nodes
              !
              call ecoute('got_reabcs')
              do while(words(1)/='ENDON')
                 ipoin                = int(param(1))
                 kfl_fixno_got(ipoin) = int(param(2))
                 if(kfl_fixno_got(ipoin)==1) then
                    if(kfl_probl_got==1) then
                       bvess_got(ndofn_got(3),ipoin)=alpha
                       do idime=1,ndime
                          bvess_got(idime,ipoin)=vnorm*veldr(idime)
                       end do
                    else if(kfl_probl_got==2) then
                       do idime=1,ndime
                          bvess_got(idime,ipoin)=vnorm*veldr(idime)
                       end do
                    else if(kfl_probl_got==3) then
                       bvess_got(1,ipoin)=alpha
                    end if
                 end if
                 call ecoute('got_reabcs')
              end do
           end if
           call ecoute('got_reabcs')
        end do
        !
        ! Initial conditions
        !
        if(kfl_probl_got==1) then
           if(inibc==4) then
              do ipoin=1,npoin
                 if(  .not.(coord(1,ipoin)<xmini_got(1).or.coord(1,ipoin)>xmaxi_got(1).or.&
                      &     coord(2,ipoin)<xmini_got(2).or.coord(2,ipoin)>xmaxi_got(2).or.&
                      &     coord(ndime,ipoin)<xmini_got(ndime).or.coord(ndime,ipoin)>xmaxi_got(ndime))) then                 
                    bvess_got(ndofn_got(3),ipoin)=alpha
                 else
                    bvess_got(ndofn_got(3),ipoin)=0.0_rp
                 end if
              end do
           else
              do ipoin=1,npoin
                 if(kfl_fixno_got(ipoin)==0)&
                    bvess_got(ndofn_got(3),ipoin)=alpha
              end do
           end if
        else if(kfl_probl_got==3) then
           if(inibc==4) then
              do ipoin=1,npoin
                 if(  .not.(coord(1,ipoin)<xmini_got(1).or.coord(1,ipoin)>xmaxi_got(1).or.&
                      &     coord(2,ipoin)<xmini_got(2).or.coord(2,ipoin)>xmaxi_got(2).or.&
                      &     coord(ndime,ipoin)<xmini_got(ndime).or.coord(ndime,ipoin)>xmaxi_got(ndime))) then 
                    bvess_got(1,ipoin)=0.0_rp
                 else
                    bvess_got(1,ipoin)=alpha
                 end if
              end do
           else
              do ipoin=1,npoin
                 bvess_got(1,ipoin)=alpha
              end do
           end if
        end if

        if(kfl_probl_got<=2) then
           if(inibc==1.and.associated(veloc_got)) then
              do ipoin=1,npoin
                 if(kfl_fixno_got(ipoin)==0) then
                    do idime=1,ndime
                       bvess_got(idime,ipoin)=veloc_got(idime,ipoin)
                    end do
                 end if
              end do
           else if(inibc==2) then
              do ipoin=1,npoin
                 if(kfl_fixno_got(ipoin)==0) then
                    do idime=1,ndime
                       bvess_got(idime,ipoin)=vnorm*veldr(idime)
                    end do
                 end if
              end do
           else if(inibc==3) then
              do ipoin=1,npoin
                 if(kfl_fixno_got(ipoin)==0) then
                    do idime=1,ndime
                       bvess_got(idime,ipoin)=0.0_rp
                    end do
                 end if
              end do
           else if(inibc==4) then
              do ipoin=1,npoin                 
                 if(  .not.(coord(1,ipoin)<xmini_got(1).or.coord(1,ipoin)>xmaxi_got(1).or.&
                      &     coord(2,ipoin)<xmini_got(2).or.coord(2,ipoin)>xmaxi_got(2).or.&
                      &     coord(ndime,ipoin)<xmini_got(ndime).or.coord(ndime,ipoin)>xmaxi_got(ndime))) then
                    do idime=1,ndime
                       bvess_got(idime,ipoin)=0.0_rp
                    end do
                 else
                    do idime=1,ndime
                       bvess_got(idime,ipoin)=vnorm*veldr(idime)
                    end do                    
                 end if
              end do
           end if
        end if
        
     else
        !
        ! Loop over nodes and or boundaries
        !
        call ecoute('got_reabcs')
        do while(words(1)/='ENDBO')
           if(words(1)=='ONNOD') then
              !
              ! On nodes
              !
              call ecoute('got_reabcs')
              do while(words(1)/='ENDON')
                 ipoin                = int(param(1))
                 kfl_fixno_got(ipoin) = int(param(2))
                 if(kfl_probl_got==1) then
                    do idofn=1,ndofn_got(3)
                       bvess_got(idofn,ipoin)=param(2+idofn)
                       if(bvess_got(idofn,ipoin)>bvmax(idofn)) bvmax(idofn)=bvess_got(idofn,ipoin)
                    end do
                 else if(kfl_probl_got==2) then
                    do idofn=1,ndime
                       bvess_got(idofn,ipoin)=param(2+idofn)
                       if(bvess_got(idofn,ipoin)>bvmax(idofn)) bvmax(idofn)=bvess_got(idofn,ipoin)
                    end do
                 else  if(kfl_probl_got==3) then
                    bvess_got(idofn,ipoin)=param(2+1)                    
                 end if
                 call ecoute('got_reabcs')
              end do
           end if
           call ecoute('got_reabcs')
        end do
        if(kfl_probl_got==1) then
           do ipoin=1,npoin
              if(kfl_fixno_got(ipoin)==0) then
                 do idofn=1,ndofn_got(3)
                    bvess_got(idofn,ipoin)=bvmax(idofn)
                 end do
              end if
           end do
        else if(kfl_probl_got==2) then
           do ipoin=1,npoin
              if(kfl_fixno_got(ipoin)==0) then
                 do idofn=1,ndime
                    bvess_got(idofn,ipoin)=bvmax(idofn)
                 end do
              end if
           end do
        else if(kfl_probl_got==3) then
           do ipoin=1,npoin
              if(kfl_fixno_got(ipoin)==0) then
                 bvess_got(1,ipoin)=bvmax(1)
              end if
           end do
        end if
     end if

  end if

end subroutine got_reabcs
