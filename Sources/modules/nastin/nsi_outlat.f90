subroutine nsi_outlat(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_outlat
  ! NAME 
  !    nsi_outlat
  ! DESCRIPTION
  !    This routine writes info on the incompressible Navier-Stokes 
  !    equations in latex format
  ! USED BY
  !    nsi_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use def_kermod, only       :  turmu_ker
  implicit none
  integer(ip), intent(in) :: itask
  character(300)          :: equat
  character(20)           :: lvisc
  integer(ip)             :: ierhs
  character(500)          :: cgnup

  if( kfl_latex==1 .and. INOTSLAVE ) then

     select case(itask)

     case(1_ip)
        !
        ! Heading
        !
        write(lun_latex,1) wback//'subsection{'//namod(modul)//&
             ': incompressible Navier-Stokes equations}'
        !
        ! Physical problem
        !
        write(lun_latex,1) wback//'subsubsection{Physical problem}'
        !
        ! Description
        !
        write(lun_latex,1) 'The equations solved are the'
        if(kfl_timei_nsi==1) then
           equat='transient'
        else
           equat='stationary'
        end if
        if(kfl_visco_nsi/=0) then
           if(kfl_advec_nsi==0) then
              equat=trim(equat)//' Stokes equations'
           else
              if(turmu_ker % kfl_exist /= 0_ip) then
                 equat=trim(equat)//' Navier-Stokes with turbulence model'
              else
                 equat=trim(equat)//' Navier-Stokes equations'
              end if
           end if
        else
           equat=trim(equat)//' Euler equations'     
        end if
        write(lun_latex,1) trim(equat)//'.'
        equat=''
        if(fvnoa_nsi>zensi) then
           equat='They are solve in an accelerated frame of reference'&
                //' with angular velocity '//wback//'mbox{'//wback//'boldmath $'//wback//'omega$}'
           if(fvnol_nsi>zensi) then
              equat=trim(equat)//' and linear acceleration '//wback//'mbox{'//wback//'boldmath $a$}.'
           else
              equat=trim(equat)//'.'
           end if
           write(lun_latex,1) trim(equat)
        else if(fvnol_nsi>zensi) then
           write(lun_latex,1) 'They are solve in an accelerated frame of reference'&
                //' with linear acceleration '//wback//'mbox{'//wback//'boldmath $a$}.'
        end if
        if(kfl_cotem_nsi/=0) then
           write(lun_latex,1) 'They are coupled with the temperature field through the'&
                //' Boussinesq approximation.'
        end if
        write(lun_latex,1) 'They read:'
        !
        ! LHS
        !
        write(lun_latex,1) wback//'begin{eqnarray*}'
        if(kfl_timei_nsi==1) write(lun_latex,1) wback//'rho '//wback&
             //'frac{'//wback//'partial '//wback//'mbox{'//wback&
             //'boldmath $u$}}{'//wback//'partial t}'
        if(kfl_advec_nsi>=1) then
           if(int(fcons_nsi)==0) then
              write(lun_latex,1) &
                   '+ '//wback//'rho ('//wback//'mbox{'//wback//'boldmath $u$} '&
                   //wback//'cdot '//wback//'nabla)'//wback//'mbox{'//wback//'boldmath $u$}'
           else if(abs(fcons_nsi-0.5_rp)<=zensi) then
              write(lun_latex,1) &
                   '+ '//wback//'rho [('//wback//'mbox{'//wback//'boldmath $u$} '&
                   //wback//'cdot '//wback//'nabla)'//wback//'mbox{'&
                   //wback//'boldmath $u$}'//wback//'frac{1}{2}'&
                   //wback//'mbox{'//wback//'boldmath $u$}('&
                   //wback//'nabla '//wback//'cdot '//wback//'mbox{'&
                   //wback//'boldmath $u$})]'
           else 
              write(lun_latex,1) &
                   '+ '//wback//'rho [('//wback//'mbox{'//wback//'boldmath $u$} '&
                   //wback//'cdot '//wback//'nabla)'//wback//'mbox{'//wback//'boldmath $u$}'//&
                   ''//wback//'mbox{'//wback//'boldmath $u$} ('//wback//'nabla '&
                   //wback//'cdot '//wback//'mbox{'//wback//'boldmath $u$})]'
           end if
        end if
        if(kfl_visco_nsi/=0) then
           lvisc=''//wback//'mu'
           if( turmu_ker % kfl_exist /= 0_ip ) lvisc='('//wback//'mu+'//wback//'mu_t)'
           if( int(fvins_nsi)==0_ip ) then
              write(lun_latex,1) '-'//wback//'nabla '//wback//'cdot['//trim(lvisc)&
                   //''//wback//'nabla('//wback//'mbox{'//wback//'boldmath $u$})] '
           else if(int(fvins_nsi)==1) then
              write(lun_latex,1) '-'//wback//'nabla '//wback//'cdot [2'//trim(lvisc)&
                   //''//wback//'mbox{'//wback//'boldmath $'//wback//'varepsilon$}'&
                   //'('//wback//'mbox{'//wback//'boldmath $u$)}] '
           else
              write(lun_latex,1) '-'//wback//'nabla '//wback//'cdot [2'//trim(lvisc)&
                   //'('//wback//'mbox{'//wback//'boldmath $'//wback//'varepsilon$}'&
                   //'({'//wback//'boldmath u)} -'//wback//'frac{1}{3}'&
                   //wback//'nabla '//wback//'cdot '//wback//'mbox{'//wback//'boldmath $u$})]'
           end if
        end if
        if(fvnoa_nsi>zensi) then
           write(lun_latex,1) '+ 2 '//wback//'rho '//wback//'mbox{'//wback//'boldmath $'&
                //wback//'omega$} '//wback//'times '//wback//'mbox{'//wback//'boldmath $u$}'
        end if
        write(lun_latex,1) '+ '//wback//'nabla p '
        !!!!if(abs(poros_nsi(1,1))>zensi) write(lun_latex,1) '+ '//wback//'sigma '//wback//'mbox{'//wback//'boldmath $u$}'
        !
        ! RHS
        !
        ierhs=0
        equat=''
        if(grnor_nsi>zensi)  then
           ierhs=1
           equat=trim(equat)//'+'//wback//'rho '//wback//'mbox{'//wback//'boldmath $g$}'
        end if
        if(kfl_cotem_nsi/=0) then
           ierhs=1
           equat=trim(equat)//'-'//wback//'rho '//wback//'mbox{'//wback//'boldmath $g_0$}'//wback//'beta_0 (T-T_0)'
        end if
        if(fvnoa_nsi>zensi) then
           ierhs=1
           equat=trim(equat)//'- '//wback//'rho '//wback//'mbox{'//wback//'boldmath $'&
                //wback//'omega$} '//wback//'times'&
                //' ('//wback//'mbox{'//wback//'boldmath $'//wback//'omega$} '&
                //wback//'times '//wback//'mbox{'//wback//'boldmath $x$})'
           equat=trim(equat)//'- '//wback//'rho '//wback//'frac{d '//wback//'mbox{'&
                //wback//'boldmath $'//wback//'omega$}}{dt}'&
                //' '//wback//'times '//wback//'mbox{'//wback//'boldmath $x$}'
        end if
        if(fvnol_nsi>zensi) then
           ierhs=1
           equat=trim(equat)//'- '//wback//'rho '//wback//'mbox{'//wback//'boldmath $a$}'
        end if
        if( turmu_ker % kfl_exist /= 0_ip .and. kfl_grtur_nsi/=0 ) then
           ierhs=1
           equat=trim(equat)//'- '//wback//'frac{2}{3} '//wback//'rho '//wback//'nabla k'
        end if
        if(ierhs==0) then
           write(lun_latex,1) '='//wback//'mbox{'//wback//'boldmath $0$}'
        else
           write(lun_latex,1) ''//wback//''//wback//' ='//trim(equat)
        end if
        write(lun_latex,1) ''//wback//''//wback//' '//wback//'qquad'
        write(lun_latex,1) ''//wback//'nabla '//wback//'cdot '//wback//'mbox{'//wback//'boldmath $u$} = 0'
        write(lun_latex,1) ''//wback//'end{eqnarray*}'
        write(lun_latex,1) 'The physical properties are'
        write(lun_latex,1) ''//wback//'begin{eqnarray*}'
        write(lun_latex,4) '    '//wback//'rho &=& ',1.0
        write(lun_latex,4) ' '//wback//''//wback//' '//wback//'mu  &=& ',1.0
        if(grnor_nsi>zensi)  then
           write(lun_latex,4) ' '//wback//''//wback//' '//wback//'mbox{'//wback//'boldmath $g$} &=& ',grnor_nsi
           write(lun_latex,4) '[',gravi_nsi(1)
           write(lun_latex,4) ',',gravi_nsi(2)
           if(ndime==3) write(lun_latex,4) ',',gravi_nsi(3)
           write(lun_latex,1) ']'
        end if
        if(kfl_cotem_nsi/=0) then
           write(lun_latex,4) ' '//wback//''//wback//' '//wback//'beta_0 &=& ',boube_nsi
           write(lun_latex,4) ' '//wback//''//wback//'  T_0    &=& ',boutr_nsi
           write(lun_latex,4) ' '//wback//''//wback//' '//wback//'mbox{'//wback//'boldmath $g_0$} &=& ',bougr_nsi
           write(lun_latex,4) '[',gravi_nsi(1)
           write(lun_latex,4) ',',gravi_nsi(2)
           if(ndime==3) write(lun_latex,4) ',',gravi_nsi(3)
           write(lun_latex,1) ']'
        end if
        write(lun_latex,1) ''//wback//'end{eqnarray*}'
        !
        ! Numerical treatment
        !
        write(lun_latex,1) ''//wback//'subsubsection{Numerical treatment}'
        write(lun_latex,1) 'The different ingredients of the numerical treatment are:'
        write(lun_latex,1) ''//wback//'begin{itemize}'
        write(lun_latex,1) '  '//wback//'item Spatial treatment:'
        write(lun_latex,1) '  '//wback//'begin{itemize}'
        write(lun_latex,1) '    '//wback//'item Stabilization: Algebaric SubGrid Scale (ASGS)'
        if(kfl_sgsco_nsi==0) then
           write(lun_latex,1) ' without convection tracking'
        else
           write(lun_latex,1) ' with convection tracking'
        end if
        if(kfl_sgsti_nsi==0) then
           write(lun_latex,1) ' and without time tracking'
        else
           write(lun_latex,1) ' and with time tracking'
        end if
        write(lun_latex,1) ' and where the stabilization parameters are'
        write(lun_latex,1)    '          '//wback//'begin{eqnarray*}'
        equat=''//wback//'tau_1 &=& '//wback//'frac{1}{{'//wback//'displaystyle c_1 '&
             //wback//'frac{'//trim(lvisc)//'}{h^2}}'
        if(kfl_advec_nsi/=0) equat=trim(equat)//'+ {'//wback//'displaystyle c_2 '&
             //wback//'rho '//wback//'frac{|'//wback//'mbox{'//wback//'boldmath $u$}|}{h_c}}'
        if(fvnol_nsi>zensi)  equat=trim(equat)//'+ c_3 '//wback//'rho |'//wback//'mbox{'//wback//'boldmath $'//wback//'omega$}| '
        !if(poros_nsi(1,1)>zensi)      equat=trim(equat)//'+ c_3 '//wback//'sigma '
        equat=trim(equat)//'} '//wback//''//wback//''
        equat=trim(equat)//''//wback//'tau_2 &=& '//wback//'frac{c_4}{c_1 '//wback//'mu +c_2 '&
             //wback//'rho h_c |'//wback//'mbox{'//wback//'boldmath $u$}|'
        !if(poros_nsi(1,1)>zensi) equat=trim(equat)//'+c_3 '//wback//'rho h^2 |'//wback//'mbox{'&
        !     //wback//'boldmath $'//wback//'omega$}| '
        if(kfl_penal_nsi/=0)   equat=trim(equat)//'+'//wback//'epsilon'
        equat=trim(equat)//'}'
        write(lun_latex,1) trim(equat)
        write(lun_latex,1) '          '//wback//'end{eqnarray*}'
        if(kfl_ellen_nsi==1) then
           write(lun_latex,1) '          with $h$ and $h_c$ being the maximum element lengths'
        else if(kfl_ellen_nsi==0) then
           write(lun_latex,1) '          with $h$ and $h_c$ being the minimum element lengths'
        else if(kfl_ellen_nsi==2) then
           write(lun_latex,1) '          with $h$ and $h_c$ being the average element lengths'
        else if(kfl_ellen_nsi==3) then
           write(lun_latex,1) '          with $h$ and $h_c$ being the minimum element length'
           write(lun_latex,1) '          and the element length in the flow direction, respectively,'
        end if
        write(lun_latex,1) '          with the following values for the constants'
        write(lun_latex,1) '          '//wback//'begin{eqnarray*}'
        write(lun_latex,4) '              c_1 &=& ',staco_nsi(1)
        write(lun_latex,4) '           '//wback//''//wback//' c_2 &=& ',staco_nsi(2)
        write(lun_latex,4) '           '//wback//''//wback//' c_3 &=& ',staco_nsi(3)
        write(lun_latex,4) '           '//wback//''//wback//' c_4 &=& ',staco_nsi(4)
        write(lun_latex,1) '          '//wback//'end{eqnarray*}'
        if(kfl_penal_nsi==1) then
           write(lun_latex,4) '    '//wback//'item Penalization: classical with $'//wback//'epsilon=$',penal_nsi
           write(lun_latex,1) '          The continuity equation is replaced by'
           write(lun_latex,1) '          '//wback//'begin{eqnarray*}'
           write(lun_latex,1) '          '//wback//'epsilon p^i + '//wback//'nabla '//wback//'cdot '&
                //wback//'mbox{'//wback//'boldmath $u$}^i = 0'
           write(lun_latex,1) '          '//wback//'end{eqnarray*}'     
        else if(kfl_penal_nsi==2) then
           write(lun_latex,4) '    '//wback//'item Penalization: iterative with $'//wback//'epsilon=$',penal_nsi
           write(lun_latex,1) '          The continuity equation is replaced by'
           write(lun_latex,1) '          '//wback//'begin{eqnarray*}'
           write(lun_latex,1) '          '//wback//'epsilon p^i + '//wback//'nabla '//wback//'cdot '&
                //wback//'mbox{'//wback//'boldmath $u$}^i'&
                //'= '//wback//'epsilon p^{i-1}'
           write(lun_latex,1) '          '//wback//'end{eqnarray*}'     
        end if
        write(lun_latex,1) '  '//wback//'end{itemize}' 
        if(kfl_timei_nsi==1) then
           write(lun_latex,1) '  '//wback//'item Temporal treatment:'
           write(lun_latex,1) '  '//wback//'begin{itemize}'
           write(lun_latex,6,advance='no') '    '//wback//'item Time integration strategy:'
           if(kfl_tisch_nsi==1) then
              write(lun_latex,1) ' Trapezoidal'
           else if(kfl_tisch_nsi==2) then
              write(lun_latex,1) ' Backward finite difference'
           end if
           write(lun_latex,7) '     '//wback//'item Order of integration: ',kfl_tiacc_nsi
           write(lun_latex,1) ''//wback//'end{itemize}'
        end if
        write(lun_latex,1) ''//wback//'end{itemize}'
        !
        ! Convergence
        !
        write(lun_gnupl,1) 'reset'
        write(lun_gnupl,1) 'set xlabel ''Number of iterations'''
        write(lun_gnupl,1) 'set ylabel ''Residual'''
        write(lun_gnupl,1) 'set log y'
        cgnup='plot '''//trim(fil_conve_nsi)//''' u 16 t ''Momentum   (tran)'' w l,&
             &      '''//trim(fil_conve_nsi)//''' u 17 t ''Continuity (tran)'' w l,&
             &      '''//trim(fil_conve_nsi)//''' u 23 t ''Momentum   (stat)'' w l,&
             &      '''//trim(fil_conve_nsi)//''' u 24 t ''Continuity (stat)'' w l'
        if(kfl_sgsco_nsi==1.or.kfl_sgsti_nsi==1) then
           cgnup=trim(cgnup)//', '''//trim(fil_conve_nsi)//''' u 7 t ''Velocity SGS'' w l'
        end if
        write(lun_gnupl,1) trim(cgnup)
        write(lun_gnupl,1) 
        write(lun_gnupl,1) 'set term post eps noenhanced color '&
             //'dashed defaultplex ''Helvetica'' 24'
        write(lun_gnupl,1) 'set output ''latex-nastin-cvg.ps'''
        write(lun_gnupl,1) 'replot'
        write(lun_gnupl,1) 'set term unknown'
        write(lun_gnupl,1) 'set output'
        write(lun_gnupl,1) 'replot'
        !
        ! Solver
        !
        write(lun_gnupl,1) 'reset'
        write(lun_gnupl,1) 'set xlabel ''Number of iterations'''
        write(lun_gnupl,1) 'set ylabel ''Number of iterations'''
        cgnup='plot '''//trim(solve(1)%fil_solve)//''' u 2 t ''Momentum  '' w l,&
             &      '''//trim(solve(2)%fil_solve)//''' u 2 t ''Continuity'' w l'
        write(lun_gnupl,1) trim(cgnup)
        write(lun_gnupl,1) 
        write(lun_gnupl,1) 'set term post eps noenhanced color '&
             //'dashed defaultplex ''Helvetica'' 24'
        write(lun_gnupl,1) 'set output ''latex-nastin-sol.ps'''
        write(lun_gnupl,1) 'replot'
        write(lun_gnupl,1) 'set term unknown'
        write(lun_gnupl,1) 'set output'
        write(lun_gnupl,1) 'replot'


     case(2_ip)

        write(lun_latex,1) 'NASTIN convergence is shown in Figure '//wback//'ref{fig:NASTIN-convergence}.'
        write(lun_latex,1) 'NASTIN solvers numbers of iterations is shown in Figure '//wback//'ref{fig:NASTIN-solver}.'
        write(lun_latex,1) ''//wback//'begin{figure}'
        write(lun_latex,1) ''//wback//'centerline{'//wback//'psfig{figure=latex-nastin-cvg.ps,width=0.95'&
             &               //wback//'textwidth}}'
        write(lun_latex,1) ''//wback//'caption{NASTIN convergence.}'
        write(lun_latex,1) ''//wback//'label{fig:NASTIN-convergence}'
        write(lun_latex,1) ''//wback//'end{figure}'

        write(lun_latex,1) ''//wback//'begin{figure}'
        write(lun_latex,1) ''//wback//'centerline{'//wback//'psfig{figure=latex-nastin-sol.ps,width=0.95'&
             &               //wback//'textwidth}}'
        write(lun_latex,1) ''//wback//'caption{NASTIN solvers numbers of iterations.}'
        write(lun_latex,1) ''//wback//'label{fig:NASTIN-solver}'
        write(lun_latex,1) ''//wback//'end{figure}'

     end select

  end if
  !
  ! Formats
  !
1 format(a)
2 format(e12.6)
3 format(i6)
4 format(a,e12.6)
5 format(a,e12.6,a)
6 format(a)
7 format(a,i1)

end subroutine nsi_outlat
      
