
subroutine hlm_dcost_sites(delta,ishot)

  use def_parame
  use def_master
  use def_domain
  use def_helmoz

  implicit none

  real(rp),intent(in) :: delta
  integer(ip),intent(in) :: ishot
  complex(rp)         :: elmat(4*mnode,4*mnode),elrhs(4*mnode)   !Element matrix, element RHS
  complex(rp)         :: pvecpo(ndime,mnode)                     !Primary magnetic vector potential
  complex(rp)         :: pscapo(mnode)                           !Primary electric scalar potential
  integer(ip)         :: ielem,igaus,idime,inode,ii              !Indices and dimensions
  integer(ip)         :: pelty,pmate,pnode,ipoin,poin,bpoin
  integer(ip)         :: pgaus,plapl,porde,ptopo
  real(rp)            :: elcod(ndime,mnode)
  real(rp)            :: gpvol(mgaus)                            !w * |J|(gaus point)
  complex(rp)         :: gprea(ncond_hlm,mgaus),dgprea(ncond_hlm,mgaus)  !Reaction terms
  complex(rp)         :: gprhs(4,mnode)                          !f (all terms)
  real(rp)            :: gpcar(ndime,mnode,mgaus)                !dNk/dxj ...
  real(rp)            :: gphes(ntens,mnode,mgaus)                !dNk/dxidx ...

  integer(ip)         :: iobs
  real(rp),dimension(3,nsite_hlm) :: vobs
  complex(rp),dimension(nsite_hlm) :: dobs
  complex(rp) :: tmpobs1, tmpobs2,tmpobs3,tmpobs4
  integer(ip) :: iobsok,ipoinok,iclosnode
  real(rp) :: rclosdistance, rnodedistance, shot 

  !call hlm_loadobs_K2(nsite_hlm,shot)

  shot = xoffsv_hlm(ishot) 

  do ii=1,size(dcostx)
     dcostx(ii)=cmplx(0.0_rp,0.0_rp,kind=rp)
  end do

  elements: do ielem = 1,nelem

     !Element dimensions
     pelty = ltype(ielem)       !Element type	
     pnode = nnode(pelty)       !Number of element nodes 
     pgaus = ngaus(pelty)       !Number of element Gauss points
     plapl = 0_ip               !Existence of Laplasian in formulation	
     porde = lorde(pelty)       !Element order
     ptopo = ltopo(pelty)       !Element topology
     
     !Check material
     pmate = 1_ip
     if ( nmate > 1_ip ) then
        pmate = lmate(ielem)
     end if
     
     iobsok=0_ip

     !Gather
     do inode = 1,pnode
        ipoin = lnods(inode,ielem)
        do iobs=1,nsite_hlm
           iclosnode =clsite1_hlm(iobs) 
           rclosdistance =  (xoffs_hlm - clsite_hlm(1,iobs))**2 + (yoffs_hlm - clsite_hlm(2,iobs))**2 + (zoffs_hlm - clsite_hlm(3,iobs))**2 
           rnodedistance =  1.0!((xoffs_hlm - coord(1, ipoin))**2 + (yoffs_hlm - coord(2, ipoin))**2 + (zoffs_hlm - coord(3, ipoin))**2 )**2


           !if(shot==-2000_rp)then
           !if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
           !)then
           !   rclosdistance = 0.0_rp
           !else
           !   if( (clsite_hlm(1,iobs)>shot-3000_rp .and. clsite_hlm(1,iobs)<shot+1000_rp) &!.or. &
           !   )then
           !      rclosdistance = 0.0_rp
           !   else
           !   rclosdistance = 1.0_rp
           !   !rclosdistance = sqrt(rclosdistance)**3  
           !   end if
           !end if
           !end if

           !if(shot==-1500_rp)then
           !if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
           !)then
           !   rclosdistance = 0.0_rp
           !else
           !   if( (clsite_hlm(1,iobs)>shot-3000_rp .and. clsite_hlm(1,iobs)<shot+1000_rp) &!.or. &
           !   )then
           !      rclosdistance = 0.0_rp
           !   else
           !   rclosdistance = 1.0_rp
           !   end if
           !end if
           !end if

           !if(shot==-1000_rp)then
           !if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
           !)then
           !   rclosdistance = 0.0_rp
           !else
           !   if( (clsite_hlm(1,iobs)>shot-3000_rp .and. clsite_hlm(1,iobs)<shot+1000_rp) &!.or. &
           !   )then
           !      rclosdistance = 0.0_rp
           !   else
           !   rclosdistance = 1.0_rp
           !   end if
           !end if
           !end if


           !if(shot==-500_rp)then
           !if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
           !)then
           !   rclosdistance = 0.0_rp
           !else
           !   if( (clsite_hlm(1,iobs)>shot-3000_rp .and. clsite_hlm(1,iobs)<shot+1000_rp) &!.or. &
           !   )then
           !      rclosdistance = 0.0_rp
           !   else
           !   rclosdistance = 1.0_rp
           !   end if
           !end if
           !end if

           !if(shot==500_rp)then
           !if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
           !)then
           !   rclosdistance = 0.0_rp
           !else
           !   if( (clsite_hlm(1,iobs)>shot-1000_rp .and. clsite_hlm(1,iobs)<shot+3000_rp) &!.or. &
           !   )then
           !      rclosdistance = 0.0_rp
           !   else
           !   rclosdistance = 1.0_rp
           !   end if
           !end if
           !end if

           !if(shot==1000_rp)then
           !if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
           !)then
           !   rclosdistance = 0.0_rp
           !else
           !   if( (clsite_hlm(1,iobs)>shot-1000_rp .and. clsite_hlm(1,iobs)<shot+3000_rp) &!.or. &
           !   )then
           !      rclosdistance = 0.0_rp
           !   else
           !   rclosdistance = 1.0_rp
           !   end if
           !end if
           !end if

           !if(shot==1500_rp)then
           !if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
           !)then
           !   rclosdistance = 0.0_rp
           !else
           !   if( (clsite_hlm(1,iobs)>shot-1000_rp .and. clsite_hlm(1,iobs)<shot+3000_rp) &!.or. &
           !   )then
           !      rclosdistance = 0.0_rp
           !   else
           !   rclosdistance = 1.0_rp
           !   end if
           !end if
           !end if

           !if(shot==2000_rp)then
           !if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
           !)then
           !   rclosdistance = 0.0_rp
           !else
           !   if( (clsite_hlm(1,iobs)>shot-1000_rp .and. clsite_hlm(1,iobs)<shot+3000_rp) &!.or. &
           !   )then
           !      rclosdistance = 0.0_rp
           !   else
           !   rclosdistance = 1.0_rp
           !   !rclosdistance = sqrt(rclosdistance)**3  
           !   end if
           !end if
           !end if

           rclosdistance = real(incidence_obs(ishot,iobs))

           if( lninv_loc(ipoin) ==  iclosnode ) then       

              !tmpobs1 = cmplx(0.0_rp,-1.0_rp*delta*rnodedistance*rclosdistance*(1.0_rp/real(countobs(iobs)))*(sum_obs(iobs) - aimag(smgvpX_obs(iobs)))) ! funciona ok
              tmpobs1 = cmplx(0.0_rp,-1.0_rp*delta*rnodedistance*rclosdistance*(1.0_rp/real(countobs(iobs)))*(sum_obs(iobs) - aimag(smgvpX_obs(ishot,iobs)))) 
              !tmpobs1 = cmplx(0.0_rp,(1.0_rp)*1.0_rp*delta*rnodedistance*rclosdistance *(aimag(smgvp_hlm(1,ipoin))/real(countobs(iobs)) - aimag(smgvpX_obs(iobs)))) 
              !tmpobs2 = (1.0_rp)*1.0_rp*delta*rnodedistance*rclosdistance *cmplx(0.0_rp,aimag(smgvp_hlm(2,ipoin) - smgvpY_obs(iobs))) 
              !tmpobs3 = (1.0_rp)*1.0_rp*delta*rnodedistance*rclosdistance *cmplx(0.0_rp,aimag(smgvp_hlm(3,ipoin) - smgvpZ_obs(iobs))) 

              !tmpobs1 =  0.0_rp 
              tmpobs2 =  cmplx(0.0_rp,0.0_rp) 
              tmpobs3 =  cmplx(0.0_rp,0.0_rp) 
              tmpobs4 =  cmplx(0.0_rp,0.0_rp)                                            

              poin  = 4_ip * (ipoin - 1_ip)
              if(countobs(iobs)>0_ip)then
                 dcostx(poin+1) = dcostx(poin+1) + tmpobs1 
                 dcostx(poin+2) = dcostx(poin+2) + tmpobs2 
                 dcostx(poin+3) = dcostx(poin+3) + tmpobs3 
                 dcostx(poin+4) = dcostx(poin+4) + tmpobs4 
 
                 modulo_dcost(lninv_loc(ipoin)) = modulo_dcost(lninv_loc(ipoin)) + tmpobs1*conjg(tmpobs1) + tmpobs2*conjg(tmpobs2)+ tmpobs3*conjg(tmpobs3)+tmpobs4*conjg(tmpobs4) 
              end if
           end if
        end do
     end do
  end do elements

end subroutine hlm_dcost_sites

