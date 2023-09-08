
subroutine opt_endopt
  !------------------------------------------------------------------------
  !****f* Optsol/opt_endopt
  ! NAME
  !    opt_endopt
  ! DESCRIPTION
  ! OUTPUT
  ! USED BY
  !    Optsol
  !***
  !------------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_optsol
  use      def_inpout
  implicit none

  integer(ip) :: indvars, just_restarted
  real(rp) :: c1,c2,armijo,curvature,half,halfinv,twice,tempval
  real(rp) :: normdiffj,diffjdescdir,diffjdescdir2

  if(kfl_goopt==1) then

     c1=0.0001_rp
     c2=0.9_rp

     half=1.0_rp/2.0_rp
     twice=2.0_rp

     halfinv=2.0_rp
     normdiffj=0.0_rp  
     diffjdescdir=0.0_rp
     diffjdescdir2=0.0_rp
  
     if(kfl_curlin_opt==1) then 
        normdiffj= dot_product(diffj,diffj)
     else
        normdiffj= dot_product(diffj_prev,diffj_prev)
     end if
     normdiffj=sqrt(normdiffj)

     if(kfl_curstp_opt>=kfl_maxstp_opt .or. normdiffj<kfl_conver_opt ) then
!-------------- end of algorithm ---------------------
        kfl_goopt=0
     else 
!-------------- running algorithm ---------------------
        if(kfl_curstp_opt==1 .and. kfl_curlin_opt==1) then 
!-------------- very first iteration ---------------------

           costf_prev=costf
           costf_tmp = 1.0E+30_rp

           kfl_curlin_opt=kfl_curlin_opt + 1
           kfl_gotim=1
           kfl_goblk=1
           kfl_gocou=1
           ittim=0

           stepj = kfl_steplength_opt
           stepfactor = twice

           design_vars_prev(:) = design_vars(:)
           design_vars(:) = design_vars_tmp(:)
           if(kfl_curlin_opt>kfl_maxlin_opt) then
              kfl_goopt=0
           end if

           just_restarted=0_ip

        else 
!-------   ------- linesearch: fail, optimal point not found ---------------------
              if(kfl_curlin_opt>=kfl_maxlin_opt) then
!-------   ------- linesearch: max linesearch steps reached, restart needed ---------------------
                 if(kfl_curres_opt<kfl_maxres_opt) then              
!-------   ------- restart: success, restarting a new linesearch from different starting point ---------------------
                    if(IMASTER)then
                       print *,'CONTROL: restart'
                    end if
                   
                    costf_tmp =1.0E+30_rp
 
                    stepj = 1.0_rp

                    kfl_gotim=1
                    kfl_goblk=1
                    kfl_gocou=1
                    ittim=0
                    stepj=stepj*real(kfl_curres_opt+1)*1.618033_rp

                    kfl_curlin_opt=2
                    kfl_curres_opt=kfl_curres_opt + 1

                    do indvars=1,kfl_ndvars_opt
                       design_vars(indvars) = design_vars_prev(indvars) +  1.0E+00 * stepj * descdir_prev(indvars)
                    end do

                    !just_restarted=1_ip

                 else
!-------   ------- restart: max restart steps reached, end of algorithm ---------------------
                    if(IMASTER)then
                       print *,'CONTROL: end of line-search steps, kfl_goopt=0'
                    end if

                    kfl_goopt=0

                 end if
              else

                 if(kfl_curlin_opt==1)then

                    if(IMASTER)then
                       print *,'CONTROL: kfl_curlin=1'
                    end if

                    kfl_curlin_opt=kfl_curlin_opt + 1
                    kfl_gotim=1
                    kfl_goblk=1
                    kfl_gocou=1
                    ittim=0

                    costf_tmp = 1.0E+30_rp

                    stepfactor= twice
                    stepj = kfl_steplength_opt

                    design_vars = design_vars_tmp

                    just_restarted=1_ip

                 elseif(kfl_curlin_opt==2 .and. costf>costf_prev)then ! we need to decrease, not increase the step, or vice-versa
                 !elseif(costf>costf_prev)then ! we need to decrease, not increase the step, or vice-versa

                    if(IMASTER)then
                       print *,'CONTROL: kfl_curlin_opt==2 .and. costf>costf_prev'
                    end if

                    kfl_curlin_opt=kfl_curlin_opt + 1
                    kfl_gotim=1
                    kfl_goblk=1
                    kfl_gocou=1
                    ittim=0

                    stepfactor = half
                    stepj = stepfactor * kfl_steplength_opt

                    costf_tmp = 1.0E+30_rp

                    do indvars=1,kfl_ndvars_opt
                       design_vars(indvars) = design_vars_prev(indvars) +  1.0E+00 * stepj * scalefactor * descdir_prev(indvars)
                    end do

                    stepj = stepfactor * stepj

                 elseif(costf<costf_tmp)then ! line-search decreasing value of costf

                    if(IMASTER)then
                       print *,'CONTROL: costf<costf_tmp',costf,costf_tmp,just_restarted
                    end if

                    kfl_curlin_opt=kfl_curlin_opt + 1
                    kfl_gotim=1
                    kfl_goblk=1
                    kfl_gocou=1
                    ittim=0

                    costf_tmp = costf

                    stepj = stepfactor * stepj

                    design_vars_prev = design_vars
                    design_vars = design_vars_tmp

                 else ! line-search is not decreasing value of costf

                    if(IMASTER)then
                       print *,'CONTROL: costf>=costf_tmp',costf,costf_tmp,just_restarted
                    end if

                    costf_prev=costf_tmp
                    costf_tmp =1.0E+30_rp
                    kfl_curstp_opt=kfl_curstp_opt + 1
                    kfl_curlin_opt=1
                    kfl_gotim=1
                    kfl_goblk=1
                    kfl_gocou=1
                    ittim=0
                    stepj_prev=stepj*half
                    stepj = kfl_steplength_opt

                    just_restarted=1_ip

                    design_vars = design_vars_prev

                 end if 
              end if
        end if
     end if
  end if
  if(kfl_first_opt==1) then
     kfl_first_opt=0
  end if
  call Parall(20_ip) !MPI_Barrier
end subroutine opt_endopt
 
