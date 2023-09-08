subroutine tem_outlat(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_outlat
  ! NAME 
  !    tem_outlat
  ! DESCRIPTION
  !    This routine writes info on the heat equation in latex format
  ! USED BY
  !    tem_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_temper
  implicit none
  integer(ip),   intent(in) :: itask
  character(300)            :: equat
  character(20)             :: lvisc
  integer(ip)               :: ierhs
  character(500)            :: cgnup

  if(kfl_latex==1.and.kfl_paral<=0) then

     select case(itask)

     case(1_ip)
        !
        ! Heading
        !
        write(lun_latex,1) '\subsection{'//namod(modul)//&
             ': Heat equations}'
        !
        ! Physical problem
        !
        write(lun_latex,1) '\subsubsection{Physical problem}'
        !
        ! Description
        !
        write(lun_latex,1) 'The equation solved is the'
        if(kfl_timei_tem==1) then
           equat='transient'
        else
           equat='stationary'
        end if
        equat=trim(equat)//' Heat equation.'
        write(lun_latex,1) 'It is:'
        !
        ! LHS
        !
        write(lun_latex,1) '\begin{eqnarray*}'
        if(kfl_timei_tem==1) then
           if(kfl_regim_tem==1) then
              write(lun_latex,1) '\rho C_v \frac{d T}{dt}'
           else
              write(lun_latex,1) '\rho C_p \frac{d T}{dt}'
           end if
        end if
        if(kfl_advec_tem>=1) then
           if(kfl_regim_tem==1) then
              write(lun_latex,1) ' + \rho C_v \mbox{\boldmath $u$} \cdot \nabla T'
           else
              write(lun_latex,1) ' + \rho C_p \mbox{\boldmath $u$} \cdot \nabla T'
           end if
        end if
        lvisc='k'
        if(kfl_cotur_tem/=0) lvisc='(k+k_t)'
        write(lun_latex,1) ' - \nabla \cdot ('//trim(lvisc)//'\nabla T ) '
        if(kfl_regim_tem==1) then
           write(lun_latex,1) ' + [rho*R*div(u)]*T'
        end if
        !
        ! RHS
        !
        ierhs=0
        equat=''

        if(kfl_sourc_tem/=0)  then
           ierhs=1
           equat=trim(equat)//'+Q'
        end if
        if(kfl_joule_tem/=0) then

        end if
        if(ierhs==0) then
           write(lun_latex,1) '=\mbox{\boldmath $0$}'
        else
           write(lun_latex,1) ' ='//trim(equat)
        end if
        write(lun_latex,1) '\end{eqnarray*}'
        write(lun_latex,1) 'The physical properties are'
        write(lun_latex,1) '\begin{eqnarray*}'
        write(lun_latex,4) '    \rho &=& ', 0.0_rp
        write(lun_latex,4) ' \\  C_p &=& ', 0.0_rp
        write(lun_latex,4) ' \\    k &=& ',  0.0_rp
        if(kfl_joule_tem/=0) write(lun_latex,4) ' \\   \mu &=& ', 0.0_rp
        write(lun_latex,1) '\end{eqnarray*}'
        !
        ! Numerical treatment
        !
        write(lun_latex,1) '\subsubsection{Numerical treatment}'
        write(lun_latex,1) 'The different ingredients of the numerical treatment are:'
        write(lun_latex,1) '\begin{itemize}'
        write(lun_latex,1) '  \item Spatial treatment:'
        write(lun_latex,1) '  \begin{itemize}'
        write(lun_latex,1) '    \item Stabilization: Algebaric Sub Grid Scale (ASGS)'
        if(kfl_sgsti_tem==1) write(lun_latex,1) 'with time tracking of the subscales'                
        write(lun_latex,1)    '          \begin{eqnarray*}'
        equat='\tau_1 = \frac{1}{{\displaystyle c_1 \frac{'//trim(lvisc)//'}{h^2}}'
        if(kfl_advec_tem/=0) equat=trim(equat)//'+ {\displaystyle c_2 \rho \frac{|\mbox{\boldmath $u$}|}{h}}'
        if(react_tem>zetem)      equat=trim(equat)//'+ c_3 \sigma '
        equat=trim(equat)//'}'
        write(lun_latex,1) trim(equat)
        write(lun_latex,1) '          \end{eqnarray*}'
        write(lun_latex,1) '          with the fololowing values for the constants'
        write(lun_latex,1) '          \begin{eqnarray*}'
        write(lun_latex,4) '              c_1 &=& ',staco_tem(1)
        write(lun_latex,4) '           \\ c_2 &=& ',staco_tem(2)
        write(lun_latex,1) '          \end{eqnarray*}'
        write(lun_latex,1) '  \end{itemize}' 
        if(kfl_timei_tem==1) then
           write(lun_latex,1) '  \item Temporal treatment:'
           write(lun_latex,1) '  \begin{itemize}'
           write(lun_latex,6, advance='no') '    \item Time integration strategy:'
           write(lun_latex,1) ' Trapezoidal'
           write(lun_latex,7) '     \item Order of integration: ',kfl_tiacc_tem
           write(lun_latex,1) '\end{itemize}'
        end if
        write(lun_latex,1) '\end{itemize}'
        !
        ! Convergence
        !
        write(lun_gnupl,1) 'reset'
        write(lun_gnupl,1) 'set xlabel ''Number of iterations'''
        write(lun_gnupl,1) 'set ylabel ''Residual'''
        write(lun_gnupl,1) 'set log y'
        cgnup='plot '''//trim(momod(modul)%fil_conve)//''' u 5 t ''Temperature'' w l'
        if(kfl_sgsti_tem==1) then
           cgnup=trim(cgnup)//', '''//trim(momod(modul)%fil_conve)//''' u 7 t ''Temperature SGS'' w l'
        end if
        write(lun_gnupl,1) trim(cgnup)
        write(lun_gnupl,1) 
        write(lun_gnupl,1) 'set term post eps noenhanced color '&
             //'dashed defaultplex ''Helvetica'' 24'
        write(lun_gnupl,1) 'set output ''latex-'//trim(namod(modul))//'-cvg.ps'''
        write(lun_gnupl,1) 'replot'
        write(lun_gnupl,1) 'set term unknown'
        write(lun_gnupl,1) 'set output'
        write(lun_gnupl,1) 'replot'

     case(2_ip)

        write(lun_latex,1) trim(namod(modul))//' convergence is shown in Figure \ref{fig:'&
             //trim(namod(modul))//'-convergence}.'
        write(lun_latex,1) '\begin{figure}'
        write(lun_latex,1) '\centerline{\psfig{figure=latex-'//trim(namod(modul))//'-cvg.ps,width=0.85\textwidth}}'
        write(lun_latex,1) '\caption{'//trim(namod(modul))//' convergence.}'
        write(lun_latex,1) '\label{fig:'//trim(namod(modul))//'-convergence}'
        write(lun_latex,1) '\end{figure}'

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

end subroutine tem_outlat

