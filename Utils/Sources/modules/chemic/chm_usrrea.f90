subroutine chm_usrrea()
  !------------------------------------------------------------------------
  !****f* partis/chm_usrrea
  ! NAME 
  !    chm_usrrea
  ! DESCRIPTION
  !    Defines the chemical reactions
  ! USES
  ! USED BY
  !    chm_reaphy
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_chemic
  use def_domain
  use mod_memchk
  implicit none
  integer(ip) :: ireac,iodes,ispec,nn,nodei,nodev
  integer(4)  :: istat

  select case( wprob_chm )

  case( 'SIMPL' )

     !-------------------------------------------------------------------
     !
     ! Simple test case
     !
     !-------------------------------------------------------------------

     if( nclas_chm /= 1 ) call runend('SIMPLE REACTION: WRONG NUMBER OF CLASSES')
     if( nodes_chm /= 2 ) call runend('SIMPLE REACTION: WRONG NUMBER OF ODES')
     do ireac = 1,2
        allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
     end do
     !
     ! I  + I  <=> I2
     !
     lreac_chm( 1 ) %l ( 1 ) =  1  
     lreac_chm( 1 ) %l ( 2 ) =  1
     lreac_chm( 1 ) %l ( 3 ) = -2
     radiu_chm( 1 )          =  1.0_rp  
     !
     ! I  + I2 <=> I3
     !
     lreac_chm( 2 ) %l ( 1 ) =  1 
     lreac_chm( 2 ) %l ( 2 ) =  2
     lreac_chm( 2 ) %l ( 3 ) = -3
     radiu_chm( 2 )          =  1.0_rp  

     if( 2 /= nreac_chm ) call runend('WRONG NUMBER OF REACTIONS')
     
  case( 'COIN1' )

     !-------------------------------------------------------------------
     !
     ! Christophe Ortiz: COIN1D
     !
     !-------------------------------------------------------------------

     !
     ! 1. I + I <=> I2
     !
     ireac = 1 
     allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
     lreac_chm( ireac ) %l ( 1 ) =  1
     lreac_chm( ireac ) %l ( 2 ) =  1
     lreac_chm( ireac ) %l ( 3 ) = -2
     radiu_chm( ireac )          =  2.0_rp * radbi_chm 
     !
     ! 2. I + In <=> In+1
     !     
     do iodes = 1,nodes_chm-1
        ispec = iodes + nclas_chm
        ireac = ireac + 1
        allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
        lreac_chm( ireac ) %l ( 1 ) =  1
        lreac_chm( ireac ) %l ( 2 ) =  ispec
        lreac_chm( ireac ) %l ( 3 ) = -( ispec + 1 )
        radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(ispec)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp)

     end do

     if( ireac /= nreac_chm ) call runend('WRONG NUMBER OF REACTIONS')

  case( 'INVN1' )

     !-------------------------------------------------------------------
     !
     ! Christophe Ortiz: INVN1D
     !
     !-------------------------------------------------------------------

     !
     ! Species are organized as follows:
     ! 
     ! 1    2    3    4    5 ... nn+2  nn+3            2*nn+2
     ! I    V   I2   I3   I4 ... Inn     V1   V2   V3 ... Vnn
     ! ------   ---------------------------------------------
     ! PDE's    ODE's
     !
     nn = nodes_chm / 2
     ireac = 0
     if( 2*nn /= nodes_chm ) call runend('CHM_USRREA: WRONG NUMBER OF ODES')
     if( nclas_chm /= 2 )    call runend('CHM_USRREA: WRONG NUMBER OF CLASSES')
     !
     ! 1. I + I <=> I2
     !
     ireac = ireac + 1
     allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
     lreac_chm( ireac ) %l ( 1 ) =  1
     lreac_chm( ireac ) %l ( 2 ) =  1
     lreac_chm( ireac ) %l ( 3 ) = -3
     radiu_chm( ireac )          =  2.0_rp * radbi_chm 
     !
     ! 2. V + V <=> V2
     !
     ireac = ireac + 1
     allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
     lreac_chm( ireac ) %l ( 1 ) =  2
     lreac_chm( ireac ) %l ( 2 ) =  2
     lreac_chm( ireac ) %l ( 3 ) = -( nclas_chm + nn + 1 )
     radiu_chm( ireac )          =  2.0_rp * radbi_chm 
     !
     ! 3. I + In <=> In+1, 2 < n < n_ode-1
     !
     do iodes = 1,nn-1
        ispec = iodes + nclas_chm
        ireac = ireac + 1
        if( ireac > nreac_chm ) call runend('3. WRONG REACTION '//intost(iodes))
        allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
        lreac_chm( ireac ) %l ( 1 ) =  1
        lreac_chm( ireac ) %l ( 2 ) =  ispec
        lreac_chm( ireac ) %l ( 3 ) = -( ispec + 1 )
        radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(iodes+1)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp)
     end do
     !
     ! 4. V + In => In-1, 3 < n < n_ode
     !
     do iodes = 2,nn
        ispec = iodes + nclas_chm
        ireac = ireac + 1
        if( ireac > nreac_chm ) call runend('4. WRONG REACTION '//intost(iodes))
         allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
        lreac_chm( ireac ) %l ( 1 ) =  2
        lreac_chm( ireac ) %l ( 2 ) =  ispec
        lreac_chm( ireac ) %l ( 3 ) = -( ispec - 1 )
        radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(iodes+1)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp)
     end do
     !
     ! 5. V + I2 => I
     !
     iodes = 1
     ispec = iodes + nclas_chm
     ireac = ireac + 1
     allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
     lreac_chm( ireac ) %l ( 1 ) =  2
     lreac_chm( ireac ) %l ( 2 ) =  ispec
     lreac_chm( ireac ) %l ( 3 ) = -1
     radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(iodes+1)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp)
     !
     ! 6. V + Vn <=> Vn+1, 2 < n < n_ode-1
     !
     do iodes = 1,nn-1
        ispec = iodes + nclas_chm + nn
        ireac = ireac + 1
        if( ireac > nreac_chm ) call runend('6. WRONG REACTION '//intost(iodes))
        allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
        lreac_chm( ireac ) %l ( 1 ) =  2
        lreac_chm( ireac ) %l ( 2 ) =  ispec
        lreac_chm( ireac ) %l ( 3 ) = -( ispec + 1 )
        radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(iodes+1)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp)
     end do
     !
     ! 7. I + Vn => Vn-1, 3 < n < n_ode
     !
     do iodes = 2,nn
        ispec = iodes + nclas_chm + nn
        ireac = ireac + 1
        if( ireac > nreac_chm ) call runend('7. WRONG REACTION '//intost(iodes))
        allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
        lreac_chm( ireac ) %l ( 1 ) =  1
        lreac_chm( ireac ) %l ( 2 ) =  ispec
        lreac_chm( ireac ) %l ( 3 ) = -( ispec - 1 )
        radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(iodes+1)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp)
     end do
     !
     ! 8. I + V2 => V
     !
     iodes = 1
     ispec = iodes + nclas_chm + nn
     ireac = ireac + 1
     if( ireac > nreac_chm ) call runend('8. WRONG REACTION '//intost(iodes))
     allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
     lreac_chm( ireac ) %l ( 1 ) =  1
     lreac_chm( ireac ) %l ( 2 ) =  ispec
     lreac_chm( ireac ) %l ( 3 ) = -2
     radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(iodes+1)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp)
     !
     ! 9. I + V <=> 0
     !     
     ireac = ireac + 1
     if( ireac > nreac_chm ) call runend('9. WRONG REACTION '//intost(iodes))
     allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
     lreac_chm( ireac ) %l ( 1 ) =  1
     lreac_chm( ireac ) %l ( 2 ) =  2
     lreac_chm( ireac ) %l ( 3 ) =  0
     radiu_chm( ireac )          =  9.4379e-4_rp

     if( ireac /= nreac_chm ) call runend('WRONG NUMBER OF REACTIONS')

  case( 'INTE1' )

     !-------------------------------------------------------------------
     !
     ! Christophe Ortiz: INTE1D
     !
     !-------------------------------------------------------------------

     !
     ! Species are organized as follows:
     ! 
     ! 1    2    3    4    5 ... nn+2  nn+3            2*nn+2
     ! I    V   I2   I3   I4 ... Inn     V1   V2   V3 ... Vnn
     ! -----------   ----------------------------------------
     ! PDE's         ODE's
     !
     nodei = 298
     nodev = 299
     ireac = 0
     if( nodei + nodev /= nodes_chm ) call runend('CHM_USRREA: WRONG NUMBER OF ODES')
     if( nclas_chm /= 3 )             call runend('CHM_USRREA: WRONG NUMBER OF CLASSES')
     !
     ! 1. I + I <=> I2
     !
     ireac = ireac + 1
     allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
     lreac_chm( ireac ) %l ( 1 ) =  1
     lreac_chm( ireac ) %l ( 2 ) =  1
     lreac_chm( ireac ) %l ( 3 ) = -3
     radiu_chm( ireac )          =  2.0_rp * radbi_chm 
     !
     ! 2. V + V <=> V2
     !
     ireac = ireac + 1
     allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
     lreac_chm( ireac ) %l ( 1 ) =  2
     lreac_chm( ireac ) %l ( 2 ) =  2
     lreac_chm( ireac ) %l ( 3 ) = -( nclas_chm + nodei + 1 )
     radiu_chm( ireac )          =  2.0_rp * radbi_chm 
     !
     ! 3a. I + In <=> In+1, 2 < n < n_ode-1
     !
     ! I + I3 <=> I4
     ! I + I4 <=> I5 ...
     !
     do iodes = 1,nodei-1
        ispec = iodes + nclas_chm
        ireac = ireac + 1
        if( ireac > nreac_chm ) call runend('3. WRONG REACTION '//intost(iodes))
        allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
        lreac_chm( ireac ) %l ( 1 ) =  1
        lreac_chm( ireac ) %l ( 2 ) =  ispec
        lreac_chm( ireac ) %l ( 3 ) = -( ispec + 1 )
        radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(iodes+2)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp)
     end do
     !
     ! 3b. I2 + In => In+2, 3 < n < n_ode-2 ( I + I2 already taken into account)
     !
     ! I2 + I3 => I5
     ! I2 + I4 => I6 ...
     !
     do iodes = 1,nodei-2
        ispec = iodes + nclas_chm
        ireac = ireac + 1
        if( ireac > nreac_chm ) call runend('3. WRONG REACTION '//intost(iodes))
        allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
        lreac_chm( ireac ) %l ( 1 ) =  3
        lreac_chm( ireac ) %l ( 2 ) =  ispec
        lreac_chm( ireac ) %l ( 3 ) = -( ispec + 2 )
        radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(iodes+2)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp) &
             &                                            + (3.0_rp*real(2_ip)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp)
     end do
     !
     ! 3c. I2 + I2 => I4
     !
     ireac = ireac + 1
     allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
     lreac_chm( ireac ) %l ( 1 ) =  3
     lreac_chm( ireac ) %l ( 2 ) =  3
     lreac_chm( ireac ) %l ( 3 ) = -5
     radiu_chm( ireac )          =  2.0_rp * ( radbi_chm + (3.0_rp*real(2_ip)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp) )
     !
     ! 3d. I + I2 <=> I3
     !
     ireac = ireac + 1
     allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
     lreac_chm( ireac ) %l ( 1 ) =  1
     lreac_chm( ireac ) %l ( 2 ) =  3
     lreac_chm( ireac ) %l ( 3 ) = -4
     radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(2_ip)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp) 
     !
     ! 4. V + In => In-1, 2 < n < n_ode
     !
     ! V + I3 => I2
     ! V + I4 => I3 ...
     !
     do iodes = 1,nodei
        ispec = iodes + nclas_chm
        ireac = ireac + 1
        if( ireac > nreac_chm ) call runend('4. WRONG REACTION '//intost(iodes))
        allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
        lreac_chm( ireac ) %l ( 1 ) =  2
        lreac_chm( ireac ) %l ( 2 ) =  ispec
        lreac_chm( ireac ) %l ( 3 ) = -( ispec - 1 )
        radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(iodes+2)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp)
     end do
     !
     ! 5. V + I2 => I
     !
     ireac = ireac + 1
     allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
     lreac_chm( ireac ) %l ( 1 ) =  2
     lreac_chm( ireac ) %l ( 2 ) =  3
     lreac_chm( ireac ) %l ( 3 ) = -1
     radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(2_ip)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp)
     !
     ! 6. V + Vn <=> Vn+1, 1 < n < n_ode-1
     !
     ! V + V2 <=> V3
     ! V + V3 <=> V4 ...
     !
     do iodes = 1,nodev-1
        ispec = iodes + nclas_chm + nodei
        ireac = ireac + 1
        if( ireac > nreac_chm ) call runend('6. WRONG REACTION '//intost(iodes))
        allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
        lreac_chm( ireac ) %l ( 1 ) =  2
        lreac_chm( ireac ) %l ( 2 ) =  ispec
        lreac_chm( ireac ) %l ( 3 ) = -( ispec + 1 )
        radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(iodes+1)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp)
     end do
     !
     ! 7. I + Vn => Vn-1, 2 < n < n_ode
     !
     ! I + V3 => V2
     ! I + V4 => V3 ...
     !
     do iodes = 2,nodev
        ispec = iodes + nclas_chm + nodei
        ireac = ireac + 1
        if( ireac > nreac_chm ) call runend('7. WRONG REACTION '//intost(iodes))
        allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
        lreac_chm( ireac ) %l ( 1 ) =  1
        lreac_chm( ireac ) %l ( 2 ) =  ispec
        lreac_chm( ireac ) %l ( 3 ) = -( ispec - 1 )
        radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(iodes+1)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp)
     end do
     !
     ! 8. I + V2 => V
     !
     iodes = 1
     ispec = iodes + nclas_chm + nodei
     ireac = ireac + 1
     if( ireac > nreac_chm ) call runend('8. WRONG REACTION '//intost(iodes))
     allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
     lreac_chm( ireac ) %l ( 1 ) =  1
     lreac_chm( ireac ) %l ( 2 ) =  ispec
     lreac_chm( ireac ) %l ( 3 ) = -2
     radiu_chm( ireac )          =  2.0_rp * radbi_chm + (3.0_rp*real(2_ip)/(4.0_rp*pi*denma_chm))**(1.0_rp/3.0_rp)
     !
     ! 9. I + V <=> 0
     !     
     ireac = ireac + 1
     if( ireac > nreac_chm ) call runend('9. WRONG REACTION '//intost(iodes))
     allocate(lreac_chm( ireac ) %l ( 3 ),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LREAC_CHM','chm_usrrea',lreac_chm( ireac )%l)
     lreac_chm( ireac ) %l ( 1 ) =  1
     lreac_chm( ireac ) %l ( 2 ) =  2
     lreac_chm( ireac ) %l ( 3 ) =  0
     radiu_chm( ireac )          =  9.4379e-4_rp

     if( ireac /= nreac_chm ) call runend('WRONG NUMBER OF REACTIONS')
     
  end select

end subroutine chm_usrrea
