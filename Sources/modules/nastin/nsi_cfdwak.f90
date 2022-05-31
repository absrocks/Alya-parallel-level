!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_cfdwak.f90
!> @author  Guillaume Houzeaux
!> @brief   Read physical data
!> @details Read physical data
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_cfdwak(itask)
  use def_kintyp
  use def_parame,  only : pi,zero
  use def_elmtyp , only : ELFEM
  use def_inpout
  use def_master
  use def_nastin
  use def_domain
  use mod_memory
  use def_kermod
  use mod_elmgeo, only  : elmgeo_natural_coordinates 
  use mod_elsest, only  : elsest_host_element
  use mod_ker_proper
!  use def_postpr
  use mod_iofile
  use mod_outfor, only : outfor
  implicit none
  integer(ip), intent(in)    :: itask
  integer(ip)                :: ielem,pelty,pnode,pgaus,inode
  integer(ip)                :: idime,igaus,pmate,kmate,ipoin
  integer(ip)                :: kmate_wake,imate,dummi,ifoun
  integer(ip)                :: ptopo,jelse, ntabl
  integer(ip)                :: icoun, ipass
  real(rp)                   :: elcod(ndime,mnode), elvel(mnode, mnode)
  real(rp)                   :: xjacm(9),gpdet,gpvol
  real(rp)                   :: gpcod(3),lmini,lmaxi,dista, gpvel(3)
  real(rp)                   :: A_disk,d_ref,n_disk(3),V_disk,d_disk
  real(rp)                   :: shapf(64),deriv(3,64),coloc(3), radius
  real(rp)                   :: integ, inte2, delta, r_hub
  integer(ip), pointer, save :: lelem_material_nsi(:) => null()
  real(rp),    pointer       :: xvolu(:),xxcog(:,:),xxref(:,:)
  real(rp),    pointer       :: xvelo(:,:), xdensi(:), xxvel(:,:)
  real(rp)                   :: vehub(3), veave(3), Uinft, Cp, Ct, uhubn, uaven, power, yaw, xcowt(3)
  character(len=3)           :: messa
  character(7)               :: statu
  character(11)              :: forma
  character(6)               :: posit
  !
  ! Initialization
  !
  if( kfl_force_nsi /= 1 ) return
  kmate = max(nmate,1_ip)

  select case ( itask )

  case ( 1_ip )
     
     !----------------------------------------------------------------
     !
     ! Initialize CFD Wake (for material force term)
     !
     !---------------------------------------------------------------
     !
     ! Allocate memory
     !
     nullify(xvolu)
     nullify(xxcog)
     nullify(xxref)
     call memory_alloca(mem_modul(1:2,modul),'XVOLU','nsi_cfdwak',xvolu,kmate)
     call memory_alloca(mem_modul(1:2,modul),'XXVOG','nsi_cfdwak',xxcog,3_ip,kmate)
     call memory_alloca(mem_modul(1:2,modul),'XXREF','nsi_cfdwak',xxref,3_ip,kmate)
     !
     ! KMATE_WAKE = Number of wake materials
     !
     kmate_wake = 0
     do pmate = 1,kmate
        if( lforc_material_nsi(pmate) == 1 ) kmate_wake = kmate_wake + 1
     end do
     !
     ! Compute volume
     !
     if( INOTMASTER ) then

        do ielem = 1,nelem
           pelty = ltype(ielem)  
           pmate = 1
           if( pelty > 0 ) then
              pnode = nnode(pelty) 
              pgaus = ngaus(pelty)
              if( nmate > 1 ) pmate = lmate(ielem)
              if( lforc_material_nsi(pmate) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    do idime = 1,ndime
                       elcod(idime,inode) = coord(idime,ipoin)
                    end do
                 end do
                 do igaus = 1,pgaus
                    gpcod(1) = 0.0_rp
                    gpcod(2) = 0.0_rp
                    gpcod(3) = 0.0_rp
                    do inode = 1,pnode
                       do idime = 1,ndime
                          gpcod(idime) = gpcod(idime) + elmar(pelty) % shape(inode,igaus) * elcod(idime,inode)
                       end do
                    end do
                    call jacdet(&
                         ndime,pnode,elcod,elmar(pelty) % deriv(1,1,igaus),&
                         xjacm,gpdet)
                    gpvol = elmar(pelty) % weigp(igaus) * gpdet
                    xvolu(pmate) = xvolu(pmate) + gpvol
                    do idime = 1,ndime
                       xxcog(idime,pmate) = xxcog(idime,pmate) + gpcod(idime) * gpvol
                    end do
                 end do
              end if
           end if
        end do

     end if !INOTMASTER
     !
     ! Parall
     !
     call pararr('SUM',0_ip,kmate,xvolu)
     call pararr('SUM',0_ip,kmate*3_ip,xxcog)
     do pmate = 1,kmate
        if( lforc_material_nsi(pmate) == 1 ) then
           xxcog(1,pmate) = xxcog(1,pmate) / xvolu(pmate)
           xxcog(2,pmate) = xxcog(2,pmate) / xvolu(pmate)
           xxcog(3,pmate) = xxcog(3,pmate) / xvolu(pmate)
        end if
     end do
     !
     ! Output information
     !
     call outfor(59_ip,momod(modul) % lun_outpu,' ')
     !
     ! Save results in force parameters
     !
     do pmate = 1,kmate
        if( lforc_material_nsi(pmate) == 1 ) then
           xforc_material_nsi( 6,pmate) = xxcog(1,pmate)
           xforc_material_nsi( 7,pmate) = xxcog(2,pmate)
           xforc_material_nsi( 8,pmate) = xxcog(3,pmate)
           xforc_material_nsi(15,pmate) = xvolu(pmate)
           d_ref                        = xforc_material_nsi( 2,pmate)
           n_disk(1)                    = xforc_material_nsi( 3,pmate) 
           n_disk(2)                    = xforc_material_nsi( 4,pmate) 
           n_disk(3)                    = xforc_material_nsi( 5,pmate) 
           V_disk                       = xforc_material_nsi(15,pmate)
           d_disk                       = xforc_material_nsi( 1,pmate)
           
           xforc_material_nsi(12,pmate) = xforc_material_nsi(6,pmate) + d_ref  * n_disk(1)
           xforc_material_nsi(13,pmate) = xforc_material_nsi(7,pmate) + d_ref  * n_disk(2)
           xforc_material_nsi(14,pmate) = xforc_material_nsi(8,pmate) + d_ref  * n_disk(3)
           xforc_material_nsi(16,pmate) = 0.7_rp  ! initial or default Ct
           ioutp( 1) = pmate
           routp( 1) = V_disk/d_disk
           routp( 2) = d_ref     
           routp( 3) = n_disk(1)     
           routp( 4) = n_disk(2)     
           routp( 5) = n_disk(3)     
           routp( 6) = xforc_material_nsi( 6,pmate)
           routp( 7) = xforc_material_nsi( 7,pmate)
           routp( 8) = xforc_material_nsi( 8,pmate)
           routp( 9) = xforc_material_nsi(12,pmate)
           routp(10) = xforc_material_nsi(13,pmate)
           routp(11) = xforc_material_nsi(14,pmate)
           routp(12) = V_disk
           routp(13) = d_disk

           call outfor(60_ip,momod(modul) % lun_outpu,' ')
           
           ! if nonuniform force, dimensionalization of force distribution  
           if (ntabr_nsi(pmate).gt.0) then
              ! inquire Area and rotor radius
              A_disk = V_disk/d_disk
              radius = sqrt(A_disk/pi)
              
              ntabl =  ntabr_nsi(pmate)  ! table dimension

              ! check that data is OK
              if (abs(radiu_nsi(ntabl,pmate)-1.0_rp).gt.1.0e-2_rp) call runend('nsi_cfdwak:Total radius in rotat table should be 1.0')
              
              
              r_hub= 0.061_rp ! hub radius, prevents force going to infinite

             
              ipoin= 1
              do while (radiu_nsi(ipoin, pmate).lt.r_hub) 
                 ipoin =ipoin+1
              end do
              ! redefinition of hub radius, closest poin in table wich is greater than previous r_hub
              r_hub = radiu_nsi(ipoin, pmate)
              if (ipoin.gt.1) then 
                 do inode = 1, ipoin -1 ! points inside hub  ! scale force with radius inside the hub
                    forcn_nsi(inode, pmate) = forcn_nsi(ipoin, pmate)* radiu_nsi(inode, pmate)/r_hub
                    forct_nsi(inode, pmate) = forct_nsi(ipoin, pmate)* radiu_nsi(inode, pmate)/r_hub
                 end do
            !  else if (ipoin.eq.ntabl) then
            !     do inode = 1, ipoin -1
            !        forcn_nsi(inode, pmate) =  radiu_nsi(inode, pmate)/r_hub
            !        forct_nsi(inode, pmate) =  radiu_nsi(inode, pmate)/r_hub
            !     end do
              end if
              ! dimensionalize radius from 0 to Rotor radius
              radiu_nsi(1:ntabl, pmate) =  radiu_nsi(1:ntabl, pmate)*radius

              ! scaling of tangencial and normal forces such their integral is
              !              normal force q_n  such integral \int_0^R (q_n* dr) = A
              !              tangencial force q_t such integral \int_0^R (q_t *r *dr) = A

              ! in this way total axial force is Fn= 0.5 rho Ct U_\infty^2 Area, and
              ! Power = Torque * \Omega = 0.5 * rho*Cp U_infty^3* Area

              integ = 0.0_rp !integral of normal force q_n
              inte2 = 0.0_rp !integral of tangential force by radius q_t*r

              ! assuming linear behavior until first radial node, (according to interpolation in nsi_forcent)
              delta = radiu_nsi(1, pmate)
              integ = integ + delta*forcn_nsi(1, pmate)  

              ! exact integral considering tang force is proportional to radial position until first radial node, that is, q_t = a*r.
              ! The integral of \int_0^r1 (q_t r dr) = 1/3*a*r1^3 
              
              inte2 = inte2 + (2.0_rp/3.0_rp)*delta*forct_nsi(1, pmate)*radiu_nsi(1, pmate)  
             !  
              do inode = 1, ntabl - 1 ! number of elements
                 delta = (radiu_nsi(inode+1, pmate) - radiu_nsi(inode, pmate))
                 integ = integ + delta*(forcn_nsi(inode+1, pmate) + forcn_nsi(inode, pmate) )
                 inte2 = inte2 + delta*(forct_nsi(inode+1, pmate)*radiu_nsi(inode+1, pmate) + forct_nsi(inode, pmate)*radiu_nsi(inode, pmate) )
              end do
!             total integral 
              integ = 0.5_rp*integ
              inte2 = 0.5_rp*inte2
!             scaling of forces
              if (abs(inte2).lt.1.0e-7_rp) then
                 forct_nsi(1:ntabl, pmate) =  0.0_rp
                 forcn_nsi(1:ntabl, pmate) =  forcn_nsi(1:ntabl, pmate)* A_disk/integ
              else                 
                 forcn_nsi(1:ntabl, pmate) =  forcn_nsi(1:ntabl, pmate)* A_disk/integ
                 forct_nsi(1:ntabl, pmate) =  forct_nsi(1:ntabl, pmate)* A_disk/inte2
              end if
           end if
        end if
     end do
     !
     ! Deallocate memory
     !
     if( IMASTER ) call nsi_memphy(-5_ip)
     call memory_deallo(mem_modul(1:2,modul),'XVOLU','nsi_cfdwak',xvolu)
     call memory_deallo(mem_modul(1:2,modul),'XXVOG','nsi_cfdwak',xxcog)
     call memory_deallo(mem_modul(1:2,modul),'XXREF','nsi_cfdwak',xxref)
     !
     ! Look for host element of x_ref points
     !
     if( INOTMASTER ) then
        call memory_alloca(mem_modul(1:2,modul),'LELEM_MATERIAL_NSI','nsi_cfdwak',lelem_material_nsi,kmate)
        imate = 0
        do pmate = 1,kmate
           if( lforc_material_nsi(pmate) == 1 ) then
            
              jelse = ielse(14)
              ielse(14) = ELFEM

              call elsest_host_element(&
                ielse,relse,1_ip,meshe(ndivi),xforc_material_nsi(12:,pmate),ielem,&
                shapf,deriv,coloc,dista)

              !call elsest(&
              !     2_ip,1_ip,ielse,mnode,ndime,npoin,nelem,nnode(1:),&
              !     lnods,ltype,ltopo,coord,xforc_material_nsi(12,pmate),relse,&
              !     ielem,shapf,deriv,coloc,lelch)

              if( ielem > 0 ) then
                 imate = imate + 1
                 lelem_material_nsi(pmate) = ielem
              end if
              ielse(14)=jelse
           end if
          
        end do 
     end if !inotmaster
     call parari('SUM',0_ip,1_ip,imate)
     if( imate /= kmate_wake .and. INOTMASTER ) then
        npari =  kmate
        pari1 => lelem_material_nsi
        call par_lagran(3_ip) 
        if( INOTMASTER .and. pard1 == 1 ) then 
           do pmate = 1,kmate
              lelem_material_nsi(pmate) = max(lelem_material_nsi(pmate),0_ip)
           end do
        end if
     end if
     !
     ! open files related with wind turbines
     !
     if (INOTSLAVE) then
        if(kfl_rstar==2) then 
           statu='unknown'
           forma='formatted'
           posit='append'
        else
           statu='unknown'
           forma='formatted'
           posit='asis'
        end if
        
          do pmate = 1,kmate
             if( lforc_material_nsi(pmate) == 1 ) then
                if ((pmate-1).lt.10) then
                   messa(1:2)='00'
                   write(messa(3:3), '(i1)') pmate -1
                else if((pmate-1).lt.100) then
                   messa(1:1)='0'
                   write(messa(2:3), '(i2)') pmate -1
                else if ((pmate-1).lt.1000) then
                   write(messa(1:3), '(i3)') pmate -1
!                else 
!                   write(messa(1:4), '(i4)') pmate -1
                end if
                
                fil_windt_nsi = adjustl(trim(namda))//'.windturbine_n'//adjustl(trim(messa))
                call iofile(zero,lun_windt_nsi+pmate-2,trim(fil_windt_nsi),'NASTIN TURBINE',statu,forma,posit)
                if (kfl_rstar/=2) write (lun_windt_nsi + pmate -2, 100) 
             end if
          end do
       end if
100    format('# 1:Step  2:Time(s)   3:ThrustCoeff   4:PowerCoeff   5:Vel Infty   6:Vel Hub   7:Vel Averaged   8:Power(kW)   9:Yaw (degrees) 10: Xcoor  11: Ycoor ')
  case ( 2_ip )

     !----------------------------------------------------------------
     !
     ! Compute velocity at reference position and calculate Ct
     !
     !----------------------------------------------------------------

     nullify(xvelo)
     call memory_alloca(mem_modul(1:2,modul),'XVELO','nsi_cfdwak',xvelo,3_ip,kmate)
     nullify(xdensi)
     call memory_alloca(mem_modul(1:2,modul),'XDENSI','nsi_cfdwak',xdensi,kmate)
     nullify(xxvel)
     call memory_alloca(mem_modul(1:2,modul),'XXVEL','nsi_cfdwak',xxvel,3_ip,kmate)
   
     if( INOTMASTER ) then
        lmini = -relse(1)
        lmaxi =  relse(1) + 1.0_rp

        do pmate = 1,kmate
           if( lforc_material_nsi(pmate) == 1 ) then
              ielem = lelem_material_nsi(pmate) 
              if( ielem > 0 ) then !only one process has ielem >0
                 pelty = abs(ltype(ielem))
                 pnode = nnode(pelty)
                 ptopo = ltopo(pelty)
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    do idime = 1,ndime
                       elcod(idime,inode) = coord(idime,ipoin)
                    end do
                 end do
                 call elmgeo_natural_coordinates(      &
                      ndime,pelty,pnode,elcod,shapf,   &
                      deriv,xforc_material_nsi(12:,pmate),coloc,ifoun)

                 !call elsest_chkelm(&
                 !     ndime,ptopo,pnode,elcod,shapf,deriv,&
                 !     xforc_material_nsi(12,pmate),coloc,ifoun,lmini,lmaxi)

                 if( ifoun == 0 ) then
                    call runend('NSI_CFDWAK: Element not found for center of disc')
                 else
                    do inode = 1,pnode                 
                       ipoin = lnods(inode,ielem)
                       do idime = 1,ndime
                          xvelo(idime,pmate) = xvelo(idime,pmate) + shapf(inode) * veloc(idime,ipoin,1)
                       end do                    
                    end do
                    ! loads density
                    call ker_proper('DENSI','COG  ',dummi,ielem,xdensi(pmate:pmate),dummi,dummi,shapf,deriv)
                 end if
              end if
           end if
        end do
        ! calculates average velocities
        do ielem = 1,nelem
           pelty = ltype(ielem)  
           pmate = 1
           if( pelty > 0 ) then
              pnode = nnode(pelty) 
              pgaus = ngaus(pelty)
              if( nmate > 1 ) pmate = lmate(ielem)
              if( lforc_material_nsi(pmate) == 1 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    do idime = 1,ndime
                       elcod(idime,inode) = coord(idime,ipoin)
                       elvel(idime,inode) = veloc(idime, ipoin,1)
                    end do
                 end do
                 do igaus = 1,pgaus
                    gpvel(1:3)=0.0_rp
                    do inode = 1,pnode
                       do idime = 1,ndime
                          gpvel(idime) = gpvel(idime) + elmar(pelty) % shape(inode,igaus) * elvel(idime,inode)
                       end do
                    end do
                    call jacdet(&
                         ndime,pnode,elcod,elmar(pelty) % deriv(1,1,igaus),&
                         xjacm,gpdet)
                    gpvol = elmar(pelty) % weigp(igaus) * gpdet
                    do idime = 1,ndime
                       xxvel(idime,pmate) = xxvel(idime,pmate) + gpvel(idime) * gpvol
                    end do
                 end do
              end if
           end if
        end do
        
     end if

     call pararr('SUM',0_ip,kmate*3_ip,xxvel)
     call pararr('SUM',0_ip,kmate*3_ip,xvelo)
     call pararr('SUM',0_ip,kmate,xdensi)
     do pmate = 1,kmate
        if( lforc_material_nsi(pmate) == 1 ) &  ! averaged velocity over the disc volume
             xxvel(1:3,pmate) = xxvel(1:3,pmate) /  xforc_material_nsi( 15,pmate)
     end do

     ! only valid for incompressible flows
     ipass =0 
     ifoun =0 ! turbine found
     icoun =0 !
     do pmate=1, kmate
        if( lforc_material_nsi(pmate) == 1 ) then
           if (abs(xdensi(pmate)).lt.0.0001_rp) then ! wind turbine was not found
              if (INOTSLAVE) write(*, '(a,i4,a)') 'nsi_cfdwak: Wind turbine #', pmate-1, ' not found'
              if (ifoun/=0) then !a wind turbine has been found befor
                 xdensi(pmate) = xdensi(ifoun)
                 if (icoun/=0) xdensi(icoun) = xdensi(ifoun)
              else ! any wind turbine has already been found 
                 icoun = pmate
              end if
              
           else if (ipass==0) then!  a wind turbine was found
              ipass= 1
              ifoun = pmate
           end if
        end if
            
     end do
        
     if( INOTMASTER ) then
        do pmate = 1,kmate
           if( lforc_material_nsi(pmate) == 1 ) then
              xforc_material_nsi( 9:11,pmate) = xvelo(1:3,pmate)
              xforc_material_nsi(19:21,pmate) = xxvel(1:3, pmate) ! averaged velocity
              ! calculates Ct and power
              if (veave_nsi(1,2).gt.1.0e-2_rp) then ! calcuate actuator disk force
                                                 ! in terms of veave
                 call nsi_ctcalc_vmed(xforc_material_nsi(1, pmate), pmate ) 
              else
                 call nsi_ctcalc(xforc_material_nsi(1, pmate), pmate)
              end if
             
            
           end if

        end do
     end if     
     if( INOTSLAVE ) then ! only master writes  
        do pmate = 1,kmate
           if( lforc_material_nsi(pmate) == 1 ) then
              xforc_material_nsi( 9:11,pmate) = xvelo(1:3,pmate)
              xforc_material_nsi(19:21,pmate) = xxvel(1:3, pmate) ! averaged velocity
              ! calculates Ct and power
              if (veave_nsi(1,2).gt.1.0e-2_rp) then! calcuate actuator disk force
                                                ! in terms of veave
                 call nsi_ctcalc_vmed(xforc_material_nsi(1, pmate), pmate)
              else
                 call nsi_ctcalc(xforc_material_nsi(1, pmate), pmate)
              end if
              vehub(1:3) = xvelo(1:3,pmate)
              veave(1:3) = xxvel(1:3, pmate)
              n_disk(1:3) = xforc_material_nsi( 3:5, pmate)
              xcowt(1:3) = xforc_material_nsi(6:8,pmate)
              Uinft = xforc_material_nsi(18, pmate)
              Cp = xforc_material_nsi(17, pmate)
              Ct = xforc_material_nsi(16, pmate)
              A_disk = xforc_material_nsi(15, pmate)/xforc_material_nsi(1, pmate)
              ! Uhub
              uhubn = -dot_product(vehub(1:3),n_disk(1:3))
              ! averaged velocity
              uaven = -dot_product(veave(1:3),n_disk(1:3))
              ! power =0.5*rho*cp*U^3*A (kW)
              power = 0.5_rp*xdensi(pmate)*A_disk*Cp*&
                   Uinft*Uinft*Uinft/1000.0_rp
              ! yaw disalignment cos(yaw) =  u·n/|u|
              yaw = 180.0_rp/pi*acos(uhubn/sqrt(dot_product(vehub(1:3),vehub(1:3))))
              !
              call vecpro(vehub,n_disk,veave,ndime)
              if (veave(3).gt.0) yaw = -yaw ! change sign
             
              write(lun_windt_nsi + pmate -2 ,101) ittim, cutim, Ct, Cp, Uinft, &
                   uhubn, uaven, power, yaw, xcowt(1), xcowt(2)
           end if
        end do
     end if
      
     call memory_deallo(mem_modul(1:2,modul),'XVELO','nsi_cfdwak',xvelo)
     call memory_deallo(mem_modul(1:2,modul),'XXVEL','nsi_cfdwak',xxvel)
     call memory_deallo(mem_modul(1:2,modul),'XDENSI','nsi_cfdwak',xdensi)
  end select
  
101 format(2x,i4,3x,f11.2,1x,6(f13.4,1x),f6.1,2(f10.1,1x) )
         
end subroutine nsi_cfdwak
   
