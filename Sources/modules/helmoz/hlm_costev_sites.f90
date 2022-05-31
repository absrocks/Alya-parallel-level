subroutine hlm_costev_sites(delta,ishot)
  !------------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_costev.f90
  ! NAME 
  !    hlm_costev
  ! DESCRIPTION
  !      This routine
  !      1. Computes element matrix and element RHS for each element in a mesh;
  !      2. Assembles elemet equations into the system equations.
  ! USES
  !    hlm_assemb
  ! USED BY
  !    hlm_matrix
  !------------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_helmoz
  implicit none
  
  real(rp), intent(in) :: delta
  integer(ip), intent(in) :: ishot

  integer(ip)         :: ielem,igaus,idime,inode,ii              !Indices and dimensions
  integer(ip)         :: iboun,inodb              !Indices and dimensions
  integer(ip)         :: pelty,pmate,pnode,ipoin,poin
  integer(ip)         :: pgaus,plapl,porde,ptopo
  integer(ip)         :: pblty,pmatb,pnodb,pgaub

  real(rp)            :: accum, shot
  integer(ip)         :: iobs
  real(rp),dimension(3,nsite_hlm) :: vobs
  complex(rp),dimension(nsite_hlm) :: dobs
  real(rp),dimension(nsite_hlm) :: sumobs 
  integer(ip)              :: iclos, iclosnode,nbnodes
  real(rp) :: normcoordsite,rclosdistance, rnodedistance ,raux

  if (IMASTER) then
     nbnodes = 0_ip       !Master will allocate minimum memory for working arrays since it does not perform any calculations
  else
     nbnodes = npoin
  endif

  costf_shot=0.0_rp
  sum_obs(:)=0.0_rp



  if( INOTMASTER ) then

     shot = xoffsv_hlm(ishot) 

     !call hlm_loadobs_K2(nsite_hlm,shot)
     do iobs=1,nsite_hlm 
        sumobs(iobs)=0.0_rp
        sum_obs(iobs)=0.0_rp
     end do

     elements: do ielem = 1,nelem
        !Element dimensions
        pelty = ltype(ielem)       !Element type	
        pnode = nnode(pelty)       !Number of element nodes 

        !Gather
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
                
           do iobs=1,nsite_hlm
              iclosnode =clsite1_hlm(iobs) 
              if( lninv_loc(ipoin) ==  iclosnode ) then       
                 !write(*,*) ' selsp_obs(',iobs,')=', selsp_hlm(ipoin)   
                 !write(*,*) 'smgvpX_obs(',iobs,')=',smgvp_hlm(1,ipoin) 
                 !write(*,*) 'smgvpY_obs(',iobs,')=',smgvp_hlm(2,ipoin) 
                 !write(*,*) 'smgvpZ_obs(',iobs,')=',smgvp_hlm(3,ipoin)
                 !write(*,*) int(shot),'p(',iobs,')=', selsp_hlm(ipoin)   
                 !write(*,*) int(shot),'x(',iobs,')=',smgvp_hlm(1,ipoin) 
                 !write(*,*) int(shot),'y(',iobs,')=',smgvp_hlm(2,ipoin) 
                 !write(*,*) int(shot),'z(',iobs,')=',smgvp_hlm(3,ipoin)

                 !sum_obs(iobs) = aimag(smgvp_hlm(1,ipoin)) 
                 !sum_obs(iobs) = sum_obs(iobs) +aimag(smgvp_hlm(1,ipoin)) 


                 sum_obs(iobs) = sum_obs(iobs) +aimag(smgvp_hlm(1,ipoin))*( 1.0  / real(countobs(iobs)) ) 
                 !sum_obs(iobs) = sum_obs(iobs) +aimag(smgvp_hlm(1,ipoin)) 
              end if
           end do
        end do
     end do elements

     costf_shot=0.0_rp
     
  end if


  call pararr('SUM',0_ip,nsite_hlm,sum_obs)


  if(INOTMASTER) then
     do iobs=1,nsite_hlm
        iclosnode =clsite1_hlm(iobs) 
        rclosdistance =  (xoffs_hlm - clsite_hlm(1,iobs))**2 + (yoffs_hlm - clsite_hlm(2,iobs))**2 + (zoffs_hlm - clsite_hlm(3,iobs))**2 
        rnodedistance =  1.0!((xoffs_hlm - coord(1, ipoin))**2 + (yoffs_hlm - coord(2, ipoin))**2 + (zoffs_hlm - coord(3, ipoin))**2 )**2

        !if(shot==-2000_rp)then
        !   if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
        !   )then
        !      rclosdistance = 0.0_rp
        !   else
        !      if( (clsite_hlm(1,iobs)>shot-3000_rp .and. clsite_hlm(1,iobs)<shot+1000_rp) &!.or. &
        !      )then
        !         rclosdistance = 0.0_rp
        !      else
        !      rclosdistance = 1.0_rp
        !      !rclosdistance = sqrt(rclosdistance)**3  
        !      end if
        !   end if
        !end if

        !if(shot==-1500_rp)then
        !   if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
        !   )then
        !      rclosdistance = 0.0_rp
        !   else
        !      if( (clsite_hlm(1,iobs)>shot-3000_rp .and. clsite_hlm(1,iobs)<shot+1000_rp) &!.or. &
        !      )then
        !         rclosdistance = 0.0_rp
        !      else
        !      rclosdistance = 1.0_rp
        !      !rclosdistance = sqrt(rclosdistance)**3  
        !      end if
        !   end if
        !end if

        !if(shot==-1000_rp)then
        !   if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
        !   )then
        !      rclosdistance = 0.0_rp
        !   else
        !      if( (clsite_hlm(1,iobs)>shot-3000_rp .and. clsite_hlm(1,iobs)<shot+1000_rp) &!.or. &
        !      )then
        !         rclosdistance = 0.0_rp
        !      else
        !      rclosdistance = 1.0_rp
        !      !rclosdistance = sqrt(rclosdistance)**3  
        !      end if
        !   end if
        !end if

        !if(shot==-500_rp)then
        !   if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
        !   )then
        !      rclosdistance = 0.0_rp
        !   else
        !      if( (clsite_hlm(1,iobs)>shot-3000_rp .and. clsite_hlm(1,iobs)<shot+1000_rp) &!.or. &
        !      )then
        !         rclosdistance = 0.0_rp
        !      else
        !      rclosdistance = 1.0_rp
        !      !rclosdistance = sqrt(rclosdistance)**3  
        !      end if
        !   end if
        !end if

        !if(shot==500_rp)then
        !   if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
        !   )then
        !      rclosdistance = 0.0_rp
        !   else
        !      if( (clsite_hlm(1,iobs)>shot-1000_rp .and. clsite_hlm(1,iobs)<shot+3000_rp) &!.or. &
        !      )then
        !         rclosdistance = 0.0_rp
        !      else
        !      rclosdistance = 1.0_rp
        !      !rclosdistance = sqrt(rclosdistance)**3  
        !      end if
        !   end if
        !end if

        !if(shot==1000_rp)then
        !   if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
        !   )then
        !      rclosdistance = 0.0_rp
        !   else
        !      if( (clsite_hlm(1,iobs)>shot-1000_rp .and. clsite_hlm(1,iobs)<shot+3000_rp) &!.or. &
        !      )then
        !         rclosdistance = 0.0_rp
        !      else
        !      rclosdistance = 1.0_rp
        !      !rclosdistance = sqrt(rclosdistance)**3  
        !      end if
        !   end if
        !end if

        !if(shot==1500_rp)then
        !   if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
        !   )then
        !      rclosdistance = 0.0_rp
        !   else
        !      if( (clsite_hlm(1,iobs)>shot-1000_rp .and. clsite_hlm(1,iobs)<shot+3000_rp) &!.or. &
        !      )then
        !         rclosdistance = 0.0_rp
        !      else
        !      rclosdistance = 1.0_rp
        !      !rclosdistance = sqrt(rclosdistance)**3  
        !      end if
        !   end if
        !end if

        !if(shot==2000_rp)then
        !   if( (clsite_hlm(1,iobs)<shot-3000_rp .or. clsite_hlm(1,iobs)>shot+3000_rp) &!.or. &
        !   )then
        !      rclosdistance = 0.0_rp
        !   else
        !      if( (clsite_hlm(1,iobs)>shot-1000_rp .and. clsite_hlm(1,iobs)<shot+3000_rp) &!.or. &
        !      )then
        !         rclosdistance = 0.0_rp
        !      else
        !      rclosdistance = 1.0_rp
        !      !rclosdistance = sqrt(rclosdistance)**3  
        !      end if
        !   end if
        !end if

        rclosdistance = real(incidence_obs(ishot,iobs))
        !print *,kfl_paral,'rclosdistance',rclosdistance,incidence_obs(ishot,iobs)


        sumobs(iobs) = 0.0_rp &
             !+ delta * rnodedistance * rclosdistance * (real(smgvp_hlm(1,ipoin) - smgvpX_obs(iobs))**2 + aimag(smgvp_hlm(1,ipoin) - smgvpX_obs(iobs))**2) & 
             !+ delta * rnodedistance * rclosdistance * (real(smgvp_hlm(2,ipoin) - smgvpY_obs(iobs))**2 + aimag(smgvp_hlm(2,ipoin) - smgvpY_obs(iobs))**2) &
             !+ delta * rnodedistance * rclosdistance * (real(smgvp_hlm(3,ipoin) - smgvpZ_obs(iobs))**2 + aimag(smgvp_hlm(3,ipoin) - smgvpZ_obs(iobs))**2) & 
             !+ delta * rnodedistance * rclosdistance * (real(selsp_hlm(ipoin)   - selsp_obs(iobs) )**2 + aimag(selsp_hlm(ipoin)   - selsp_obs(iobs) )**2)  
             !+ delta * rnodedistance * rclosdistance * (1.0_rp*aimag(sum_obs(iobs) - smgvpX_obs(iobs))**2)  

             !+ delta * rnodedistance * rclosdistance * ( (sum_obs(iobs) - aimag(smgvpX_obs(iobs)))**2)  ! funciona ok


             + delta * rnodedistance * rclosdistance * ( (sum_obs(iobs) - aimag(smgvpX_obs(ishot,iobs)))**2)  
             !+ delta * rnodedistance * rclosdistance * (1.0_rp*aimag(smgvp_hlm(2,ipoin) - smgvpY_obs(iobs))**2) & 
             !+ delta * rnodedistance * rclosdistance * (1.0_rp*aimag(smgvp_hlm(3,ipoin) - smgvpZ_obs(iobs))**2) / abs(aimag(smgvpZ_obs(iobs))**2) 

     end do

     do iobs=1,nsite_hlm
        if(countobs(iobs)>0_ip)then
           costf_shot = costf_shot + sumobs(iobs)
        end if
     end do

  end if

  call pararr('MAX',0_ip,1,costf_shot)

end subroutine hlm_costev_sites

