subroutine hlm_solite()

  !------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_solite.f90
  ! NAME 
  !    hlm_solite
  ! DESCRIPTION
  !    This routine solves one iteration of the incompletely gauged coupled 
  !    vector-scalar potential formulation of Maxwell's equations.
  ! USES
  !    hlm_matrix
  !    solvex
  !    hlm_updunk
  ! USED BY
  !    hlm_doiter
  !------------------------------------------------------------------------

  use def_parame
  use def_master
  use def_solver
  use def_domain
  use def_helmoz
  use def_kermod, only       :  kfl_ndvars_opt,costf
  use mod_messages, only : livinf
  use mod_communications
  implicit none

  integer(ip) :: ishot,icomp,h1,ii,ierror,indvars,isInside,nbnodes 
  real(rp)    :: cpu_refe1,cpu_refe2,cpu_refe
  real(rp)    :: cpu_alsy1,cpu_alsy2,cpu_alsy
  real(rp)    :: time1,time2,time3,time4,time5,time6, time7, time8,timeAss,timeProd, timeMat
  real(rp)    :: weight_shot, weight, normgradinv,normunknx

  real(rp) :: delta, shot 
  complex(rp) :: rr


  if( kfl_edges_hlm == 1 ) then

     call hlm_edge_element_assembly()

  else

     costf=0.0_rp 
     weight=0.0_rp


     if(kfl_servi(ID_OPTSOL)==1) then
        diffj(:)=0.0_rp
     end if

     shots: do ishot=1,nshot_hlm

        amatx(:)=cmplx(0.0_rp,0.0_rp,kind=rp)
        unknx(:)=cmplx(0.0_rp,0.0_rp,kind=rp)
        rhsix(:)=cmplx(0.0_rp,0.0_rp,kind=rp)

        pmgvp_hlm(:,:)=cmplx(0.0_rp,0.0_rp,kind=rp) 
        pelsp_hlm(  :)=cmplx(0.0_rp,0.0_rp,kind=rp)      
        smgvp_hlm(:,:)=cmplx(0.0_rp,0.0_rp,kind=rp) 
        selsp_hlm(  :)=cmplx(0.0_rp,0.0_rp,kind=rp)      

        if(kfl_servi(ID_OPTSOL)==1) then
           damatx(:)=cmplx(0.0_rp,0.0_rp,kind=rp)
           aunknx(:)=cmplx(0.0_rp,0.0_rp,kind=rp)
           drhsix(:)=cmplx(0.0_rp,0.0_rp,kind=rp)
           dcostx(:)=cmplx(0.0_rp,0.0_rp,kind=rp)
        end if

        xoffs_hlm  = xoffsv_hlm(ishot)
        yoffs_hlm  = yoffsv_hlm(ishot)
        zoffs_hlm  = zoffsv_hlm(ishot)
        elcur_hlm  = elcurv_hlm(ishot)
        length_hlm = lengthv_hlm(ishot)

        delta = 1.0E-00_rp 
        shot=xoffs_hlm

        call cputim(cpu_refe1)

        !Update inner iteration counter
        itinn(modul) = itinn(modul) + 1
        ittot_hlm    = ittot_hlm + 1
        call livinf(15_ip,' ',modul)   ! SOLVE HELMOZ

        !Construct the linear system of algebraic equations - construct the system matrix and the system right-hand-side
        call hlm_matrix()
        !Obtain the initial guess for vector of unknowns to be determined with an iterative solver
        call hlm_updunk(1_ip)


        !Solve the linear system of algebraic equations 
        call cputim(cpu_alsy1)

        call solvex(rhsix,unknx,amatx,pmatx)
        call PAR_BARRIER() 
        call cputim(cpu_alsy2)

        !Update unknowns: smgvr_hlm (secondary magnetic vector potential) and selsp_hlm (secondary electric scalar potential)
        call hlm_updunk(2_ip)

        call livinf(16_ip,' ',itinn(modul))
        call cputim(cpu_refe2)
        cpu_modul(3,modul) =  cpu_modul(3,modul) + cpu_refe2 - cpu_refe1
        if (INOTSLAVE) then
           cpu_refe = cpu_refe2 - cpu_refe1 
           write(*,*)'Total time:',cpu_refe
        endif

        if(kfl_servi(ID_OPTSOL)==1) then
           costf_shot=0.0_rp
           modulo_dcost(:) = 0.0_rp

           call cputim(time3)
           if(kfl_first_opt==1)then
              call hlm_countobs_sites(delta,shot)
           end if
           call hlm_costev_sites(delta,ishot)
           call cputim(time4)
           write(*,*)'kfl_paral=',kfl_paral,' iteration costev time: ',time4-time3
           if(kfl_curlin_opt==1)then
              call cputim(time5)
              if( INOTMASTER ) then
                 call hlm_dcost_sites(delta,ishot)
              end if
              call cputim(time6)
              write(*,*)'kfl_paral=',kfl_paral,' iteration dcost time: ',time6-time5

           end if
           call pararr('MAX',0_ip,kfl_ndvars_opt,diffj_illum)
           call pararr('SUM',0_ip,kfl_ndvars_opt,modulo_dcost)

           print *, 'kfl_paral=',kfl_paral,', costf_shot=',costf_shot

           costfv_hlm(ishot) = costf_shot

           !if(kfl_curlin_opt==1 .and. 1_ip==0_ip) then
           if(kfl_curlin_opt==1) then
              ! gradient calculation

              ! calculate dj/dgamma(gamma_k)
              do indvars=1,kfl_ndvars_opt
                 diffj_shot(indvars)=0.0_rp
              end do

              weightv_hlm(ishot)=0.0_rp
              rr=cmplx(0.0_rp,0.0_rp,kind=rp)

              call hlm_adjvar()

              timeAss=0.0_rp
              timeProd=0.0_rp
              timeMat=0.0_rp
              if(IMASTER)then
                 open(unit=10, file="dcostnew.txt", action='WRITE')
              end if

              if(IMASTER) then
                 call cputim(time1)
                 write(10,*)'indvars=',0_ip,'/', kfl_ndvars_opt,time1
              end if
              call cputim(time3)
              do indvars=1,kfl_ndvars_opt        
                 isInside=diffj_isInside(indvars)
                 if(IMASTER .and. modulo(indvars,10000_ip)==0_ip) then
                    call cputim(time1)
                    write(10,*)'indvars=',indvars,'/', kfl_ndvars_opt,' modulo_dcost=',modulo_dcost(indvars)
                 end if
                 if(isInside>0_ip)then
                    drhsix(:)=cmplx(0.0_rp,0.0_rp,kind=rp)
                    damatx(:)=cmplx(0.0_rp,0.0_rp,kind=rp)
                    if(INOTMASTER) then   
                       call cputim(time5)     
                       call hlm_elmopediffAssem(indvars)
                       call hlm_dirbcsdiffAssem(indvars)
                       call cputim(time6)    
                       timeAss=timeAss + (time6-time5) 
                    end if

                    call cputim(time5)    

                    call bcsplx2( npoin, 4_ip, damatx, c_sol, r_sol, unknx, drhsix)
                    call cputim(time6)     
                    timeMat=timeMat + (time6-time5) 
                    rr=cmplx(0.0_rp,0.0_rp,kind=rp)
                    call cputim(time6)     
                    call proptx(npoin,4_ip,aunknx,drhsix,rr)
                    diffj_shot(indvars) = -2.0_rp * real(rr)  
                    call cputim(time7)     
                    timeProd=timeProd + (time7-time6)

                    diffj_shot(indvars)=  diffj_shot(indvars) * exp(design_vars(indvars)) * diffj_illum(indvars)! chain rule

                    weightv_hlm(ishot) = weightv_hlm(ishot) + diffj_shot(indvars)*diffj_shot(indvars) 

                 else
                    diffj_shot(indvars) = 0.0_rp
                 end if

              end do
              call cputim(time4)
              write(*,*)'kfl_paral=',kfl_paral,' iteration diffj_shot time: ',time4-time3
              write(*,*)'kfl_p=',kfl_paral,' it diffj_shot-timeAss time: ',timeAss
              write(*,*)'kfl_p=',kfl_paral,' it diffj_shot-timeMat time: ',timeMat
              write(*,*)'kfl_p=',kfl_paral,' it diffj_shot-timeProd time: ',timeProd

              if(IMASTER)then
                 close(10)
              end if

              diffj = diffj + diffj_shot

              weightv_hlm(ishot) = sqrt(weightv_hlm(ishot)) 
              weight_shot = weight_shot + weightv_hlm(ishot) 
              !weight_shot = 1.0_rp / (costf_shot)

              !if(ishot == nshot_hlm)then
              !   do ii=1,nshot_hlm
              !      weightv_hlm(ii) = weightv_hlm(ii) / weight_shot 
              !   end do
              !end if

              ! end gradient calculation

           end if
        end if

        !costf=costf + costf_shot
        !print *, 'kfl_paral=',kfl_paral,', w_costf_shot=',(weightv_hlm(ishot) *costf_shot)
        !costf=costf + 1.0E+10 * weight_shot * costf_shot
        !costf=costf + weightv_hlm(ishot) * costf_shot

        !weight_shot = sqrt(dot_product(diffj_shot,diffj_shot)) 

        !diffj = diffj + diffj_shot
        !diffj = diffj + weight_shot * diffj_shot
        !if(ishot==1 .or. ishot==3)then
        !   diffj = diffj + 0.25 * diffj_shot
        !else if(ishot==2)then
        !   diffj = diffj + 0.5 * diffj_shot
        !end if
        !diffj = diffj + (1.0_rp/(weight_shot + 1.0e-30)) * diffj_shot
        !weight = weight + 1.0E+10 * weight_shot  

     end do shots


     if(kfl_servi(ID_OPTSOL)==1) then
        do ishot=1,nshot_hlm
           !weightv_hlm(ishot) = weightv_hlm(ishot) / weight_shot 
           !if(IMASTER)then
           !   print *, 'kfl_paral=',kfl_paral,', weightv_hlm=',weightv_hlm(ishot) 
           !   print *, 'kfl_paral=',kfl_paral,', w_costf_shot=',(weightv_hlm(ishot) *costfv_hlm(ishot))
           !end if
           !costf=costf + weightv_hlm(ishot) * costfv_hlm(ishot)
           costf=costf + costfv_hlm(ishot)
        end do
     end if
  end if

end subroutine hlm_solite
