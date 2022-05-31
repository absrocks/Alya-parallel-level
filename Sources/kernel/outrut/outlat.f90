subroutine outlat(itask)
  !-----------------------------------------------------------------------
  !****f* Master/outlat
  ! NAME
  !    outlat
  ! DESCRIPTION
  !    This routine writes some information about the run in a Latex file
  ! OUTPUT
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_elmtyp
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ielty,imodu,ibloc,iorde
  character(100)          :: cmodu,cbloc,ctimc
  character(6)            :: crule
  character(3)            :: cquad
  character(20)           :: cblo2
  character(500)          :: cgnup

  if( kfl_latex/=0 .and. INOTSLAVE .and. .not. READ_AND_RUN() ) then

     select case (itask)

     case(1_ip)
        !
        ! Header
        !
        write(lun_latex,90)& 
             ''//wback//'documentclass[10pt]{article}',&
             ''//wback//'usepackage[dvips]{epsfig}',&
             ''//wback//'begin{document}',&
             ''//wback//'title{Sum up of Alya run for problem: ',trim(title),'}',&
             ''//wback//'maketitle'
        write(lun_latex,1)&
             ''//wback//'section{Problem definition}'
        !
        ! General data
        !
        cmodu=''
        cbloc=''
        do imodu=1,mmodu
           if(kfl_modul(imodu)/=0) then
              cmodu=trim(cmodu)//trim(namod(imodu))//','
           end if
        end do
        cmodu=trim(cmodu(1:len(trim(cmodu))-1))
        do ibloc=1,nblok
           cbloc=trim(cbloc)//trim(intost(ibloc))//': '
           do iorde=1,mmodu
              imodu=lmord(iorde,ibloc)
              if(imodu/=0) then
                 if(kfl_modul(imodu)/=0) then
                    cbloc=trim(cbloc)//trim(namod(imodu))//','
                 end if
              end if
           end do
           cbloc=trim(cbloc(1:len(trim(cbloc))-1))
        end do
        if(kfl_timco==0) then
           ctimc='PRESCRIBED'
        else if(kfl_timco==1) then
           ctimc='MINIMUM OF CRITICAL TIME STEPS'
        else
           ctimc='LOCAL TIME STEP'
        end if

        write(lun_latex,1)&
             ''//wback//'subsection{Physical problem}'
        write(lun_latex,1) 'In this run, Alya solves the modules '//trim(cmodu)//','
        write(lun_latex,1) 'according to the following coupling strategy:'
        write(lun_latex,1) ''//wback//'begin{itemize}'
        write(lun_latex,9) '  '//wback//'item[] {'//wback//'bf do} ',mitim,' time iterations '
        write(lun_latex,1) '  '//wback//'begin{itemize}'
        do ibloc=1,nblok
           cblo2=intost(ibloc)
           write(lun_latex,7)  ''//wback//'item[] block '//trim(cblo2)//': {'//wback//'bf do} ',micou(ibloc),' coupling iterations'
           write(lun_latex,1)  ''//wback//'begin{itemize}'
           do iorde=1,mmodu 
              imodu=lmord(iorde,ibloc)
              if(imodu/=0) then
                 if(kfl_modul(imodu)/=0) then                 
                    write(lun_latex,1)  '  '//wback//'item[-] Solve module '//namod(imodu)
                 end if
              end if
           end do
           write(lun_latex,1)  ''//wback//'end{itemize}'
           write(lun_latex,1)  ''//wback//'item[] {'//wback//'bf enddo} coupling iteration'
        end do
        write(lun_latex,1) '  '//wback//'end{itemize}'  
        write(lun_latex,1) ''//wback//'item[] {'//wback//'bf enddo} time iteration'
        write(lun_latex,1) ''//wback//'end{itemize}' 
        write(lun_latex,1) 'The time interval is $['
        write(lun_latex,2) timei
        write(lun_latex,1) ','
        write(lun_latex,2) timef
        write(lun_latex,1) ']$.'
        if(kfl_timco==0) then
           write(lun_latex,1) 'The time step $'//wback//'delta t$ is constant and such that $'//wback//'delta='
           write(lun_latex,2) dtime
           write(lun_latex,1) '$.'
        else if(kfl_timco==1) then
           write(lun_latex,1) &
                'Let $'//wback//'delta t_{{'//wback//'rm cri},i}$ be the critical time steps of module $i$,'
           write(lun_latex,1) &
                'and let $'//wback//'alpha_{{'//wback//'rm saf},i}$ be the safety factor of module $i$.'
           write(lun_latex,1) &
                'The time step $'//wback//'delta t$ is taken as the minimum of the module critical time steps'
           write(lun_latex,1) &
                'multiplied by the corresponding safety factor such that'
           write(lun_latex,1) ''//wback//'begin{eqnarray*}'
           write(lun_latex,1) ''//wback//'delta t = '//wback//'min_{i} '//wback//'alpha_{{'&
                //wback//'rm saf},i} '//wback//'delta t_{{'//wback//'rm cri},i}'
           write(lun_latex,1) ''//wback//'quad '//wback//'mbox{for all modules $i$ solved}'
           write(lun_latex,1) ''//wback//'end{eqnarray*}'             
        else
           write(lun_latex,1) 'The time step is computed for each element of the domain and'
           write(lun_latex,1) 'for each module as the element critical time step multiplied'
           write(lun_latex,1) 'by the module safety factor.'
        end if
        !
        ! Domain data
        !
        write(lun_latex,1)&
             ''//wback//'subsection{Geometrical data}'
        write(lun_latex,1) 'The computational domain belong to'
        if(ndime==2) then
           write(lun_latex,1) '$R^2$.'
        else
           write(lun_latex,1) '$R^3$.'
        end if
        write(lun_latex, 1) 'Let npoin be the number of nodes, nelem the number'
        write(lun_latex, 1) 'of voulme elements and nboun the number of boundary elements.'
        write(lun_latex, 1) 'We have:'
        write(lun_latex, 1) ''//wback//'begin{eqnarray*}'
        write(lun_latex, 4) ' '//wback//'mbox{npoin} &=&',npoin,' '//wback//''//wback//''
        write(lun_latex, 4) ' '//wback//'mbox{nelem} &=&',nelem,' '//wback//''//wback//''
        write(lun_latex, 4) ' '//wback//'mbox{nboun} &=&',nboun,''
        write(lun_latex, 1) ''//wback//'end{eqnarray*}'
        write(lun_latex, 1) 'The domain data are'
        write(lun_latex, 1) ''//wback//'begin{itemize}'
        write(lun_latex,10) '   '//wback//'item Total volume= ',vodom
        write(lun_latex,10) '   '//wback//'item Averaged element volume= ',voave
        write(lun_latex,10) '   '//wback//'item Minimum element volume= ',vomin
        write(lun_latex,11) '         for element ',elmin
        write(lun_latex,10) '   '//wback//'item Maximum element volume= ',vomax
        write(lun_latex,11) '         for element ',elmax
        write(lun_latex, 1) ''//wback//'end{itemize}'
        write(lun_latex, 1) 'The elements that compose the volume and boundary meshes, as well'     
        write(lun_latex, 1) 'as their corresponding integration rules and numbers of Gauss points'
        write(lun_latex, 1) 'is given in the following table.'
        write(lun_latex, 1) ''//wback//'begin{center}'
        write(lun_latex, 1) ''//wback//'begin{math}'
        write(lun_latex, 1) ''//wback//'begin{array}{lllllll}'
        write(lun_latex, 1) ''//wback//'hline'
        write(lun_latex, 1) ''//wback//'textrm{Element} & '//wback//'textrm{Total '//wback//'#}  & '&
             //wback//'textrm{Dimension} & '//wback//'textrm{Nodes '//wback//'#} & '
        write(lun_latex, 1) ''//wback//'textrm{Integr.} & '//wback//'textrm{Gauss}     & '&
             //wback//'textrm{Laplacian} '//wback//''//wback//''
        write(lun_latex, 1) ''//wback//'textrm{}        & '//wback//'textrm{}          & '&
             //wback//'textrm{}          & '//wback//'textrm{} & '
        write(lun_latex, 1) ''//wback//'textrm{rule}    & '//wback//'textrm{points '&
             //wback//'#} & '//wback//'textrm{} '//wback//''//wback//''
        write(lun_latex, 1) ''//wback//'hline'

        do ielty=1,nelty
           if(lexis(ielty)/=0) then
              if(llapl(ielty)==1) then
                 cquad='yes'
              else
                 cquad='no'
              end if
              if(lquad(ielty)==0) then
                 crule='open'
              else
                 crule='closed'
              end if
              write(lun_latex,5) &
                   ''//wback//'textrm{'//trim(cenam(ielty))//'}',&
                   lnuty(ielty),&
                   ldime(ielty),&
                   nnode(ielty),&
                   ''//wback//'textrm{'//trim(crule)//'}',&
                   ngaus(ielty),&
                   ''//wback//'textrm{'//trim(cquad)//'}',wback//wback
           end if
        end do
        write(lun_latex,1) ''//wback//'hline'
        write(lun_latex,1) ''//wback//'end{array}'
        write(lun_latex,1) ''//wback//'end{math}'
        write(lun_latex,1) ''//wback//'mbox{Table: Volume and boundary elements}'
        write(lun_latex,1) ''//wback//'end{center}'
        write(lun_latex,1)&
             ''//wback//'section{Module data}'
        call iofile_flush_unit(lun_latex)

        !
        ! Convergence
        ! 
        write(lun_gnupl,1) 'reset'
        write(lun_gnupl,1) 'set xlabel ''Number of iterations'''
        write(lun_gnupl,1) 'set ylabel ''Residual'''
        write(lun_gnupl,1) 'set log y'
        cgnup='plot 1 not, '
        do imodu=1,mmodu
           if(kfl_modul(imodu)/=0) then
              cgnup=trim(cgnup)&
                   //' '''//adjustl(trim(fil_conve))&
                   //''' u '//trim(intost(imodu+3_ip))//' title '''&
                   //adjustl(trim(namod(imodu)))//''' w lines,'
           end if
        end do
        cgnup=trim(cgnup(1:len(trim(cgnup))-1))
        write(lun_gnupl,1) trim(cgnup)
        write(lun_gnupl,1) 
        write(lun_gnupl,1) 'set term post eps noenhanced color '&
             //'dashed defaultplex ''Helvetica'' 24'       
        write(lun_gnupl,1) 'set output ''latex-cvg.ps'''
        write(lun_gnupl,1) 'replot'
        write(lun_gnupl,1) 'set term unknown'
        write(lun_gnupl,1) 'set output'
        write(lun_gnupl,1) 'replot'
        write(lun_gnupl,1) 

     case(2_ip)
        !
        ! Results
        !
        write(lun_latex,1) ''//wback//'section{Results}'
        write(lun_latex,1) 'The global convergence is shown in Figure '//wback//'ref{fig:global-convergence}.'
        write(lun_latex,1) ''//wback//'begin{figure}'
        write(lun_latex,1) ''//wback//'centerline{'//wback//'psfig{figure=latex-cvg.ps,width=0.95'&
             //wback//'textwidth}}'
        write(lun_latex,1) ''//wback//'caption{Global convergence.}'
        write(lun_latex,1) ''//wback//'label{fig:global-convergence}'
        write(lun_latex,1) ''//wback//'end{figure}'    

     case(3_ip)
        !
        ! Finish latex file
        !
        write(lun_latex,1) ''//wback//'end{document}'

     end select

  end if

1 format(a)
2 format(e12.6)
3 format(i6)
4 format(a,i9,a)
5 format(a,' & ',i9,' & ',i1,' & ',i2,' & ',a,' & ',i2,' & ',a,a)
7 format(a,i3,a)
9 format(a,i12,a)
10 format(a,e12.6)
11 format(a,i0)


90 format(a,/,a,/,a,/,a,a,a,/,a)


end subroutine outlat
