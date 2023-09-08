subroutine bouele(nelty,ngaus,lquad,ielty,ngaub,lquab,lexib)
  !-----------------------------------------------------------------------
  !****f* Domain/bouele
  ! NAME
  !    bouele
  ! DESCRIPTION
  !    This routine defines the boundary Gauss point from the volume
  !    Gauss point as well as the type of quadrature
  ! OUTPUT
  !    LEXIB ... If type exists
  !    NGAUB ... Boundary Gauss points (if IELTY > 0)
  !    LQUAB ... Boundary element quadrature (if IELTY > 0)
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  rp
  use def_elmtyp
  implicit none
  integer(ip), intent(in)  :: nelty,ngaus,lquad,ielty
  integer(ip), intent(out) :: ngaub(nelty),lquab(nelty),lexib(nelty)
  integer(ip)              :: jelty

  jelty = abs(ielty) 

  if( ielty < 0 ) then

     if( jelty < 10 ) then 

        lexib(POINT) = 1

     else if( jelty == TRI03 ) then 

        lexib(BAR02) = 1

     else if( jelty == TRI06 ) then   

        lexib(BAR03) = 1

     else if( jelty == QUA04 ) then   

        lexib(BAR02) = 1

     else if( jelty == QUA08 ) then   

        lexib(BAR03) = 1

     else if( jelty == QUA09 ) then   

        lexib(BAR03) = 1

     else if( jelty == QUA16 ) then   

        lexib(BAR04) = 1

     else if( jelty == TET04 ) then  

        lexib(TRI03) = 1

     else if( jelty == TET10 ) then 

        lexib(TRI06) = 1

     else if( jelty == PYR05) then

        lexib(TRI03) = 1
        lexib(QUA04) = 1

     else if( jelty == PYR14 ) then  

        lexib(TRI06) = 1
        lexib(QUA08) = 1

     else if( jelty == PEN06 ) then  

        lexib(TRI03) = 1
        lexib(QUA04) = 1

     else if( jelty == PEN15 ) then 

        lexib(TRI06) = 1
        lexib(QUA08) = 1
        call runend('BOUELE: PENTA_15 ELEMENT IS NOT READY')

     else if( jelty == PEN18 ) then 

        lexib(TRI06) = 1
        lexib(QUA09) = 1
        ! call runend('BOUELE: PENTA_18 ELEMENT IS NOT READY')

     else if( jelty == HEX08 ) then   

        lexib(QUA04) = 1

     else if( jelty == HEX20 ) then  

        lexib(QUA08) = 1

     else if( jelty == HEX27 ) then 

        lexib(QUA09) = 1

     else if( jelty == SHELL ) then   

        lexib(BAR02) = 1

     else if( jelty == BAR3D ) then   

        lexib(POINT) = 1

     end if

  else

     if( jelty < 10 ) then 

        lexib(POINT) = 1
        lquab(POINT) = lquad
        ngaub(POINT) = 1

     else if( jelty == TRI03 ) then 

        lexib(BAR02) = 1
        lquab(BAR02) = lquad

        if( lquad == 0 ) then               ! open rule
           if(ngaus<=1 ) then       
              if( ngaub(BAR02) == 0 ) ngaub(BAR02) = 1                 ! exact for P1 
           else if(ngaus<=4 ) then  
              if( ngaub(BAR02) == 0 ) ngaub(BAR02) = 2                 ! exact for P3
           else if(ngaus<=7 ) then  
              if( ngaub(BAR02) == 0 ) ngaub(BAR02) = 3                 ! exact for P5
           else if(ngaus<=13 ) then 
              if( ngaub(BAR02) == 0 ) ngaub(BAR02) = 3                 ! exact for P7
           end if
        else if( lquad == 1 ) then          ! closed rule
           if(ngaus<=4) then
              if( ngaub(BAR02) == 0 ) ngaub(BAR02) = 2
           else if(ngaus<=7) then
              if( ngaub(BAR02) == 0 ) ngaub(BAR02) = 3
           else if(ngaus<=10) then
              if( ngaub(BAR02) == 0 ) ngaub(BAR02) = 4
           end if
        end if

     else if( jelty == TRI06 ) then   

        lexib(BAR03) = 1
        lquab(BAR03) = lquad

        if( lquad == 0 ) then               ! open rule
           if(ngaus<=1 ) then       
              if( ngaub(BAR03) == 0 ) ngaub(BAR03) = 1                 ! exact for P1 
           else if(ngaus<=4 ) then  
              if( ngaub(BAR03) == 0 ) ngaub(BAR03) = 2                 ! exact for P3
           else if(ngaus<=7 ) then  
              if( ngaub(BAR03) == 0 ) ngaub(BAR03) = 3                 ! exact for P5
           else if(ngaus<=13 ) then 
              if( ngaub(BAR03) == 0 ) ngaub(BAR03) = 3                 ! exact for P7
           end if
        else if( lquad == 1 ) then          ! closed rule
           if(ngaus<=4) then
              if( ngaub(BAR03) == 0 ) ngaub(BAR03) = 2
           else if(ngaus<=7) then
              if( ngaub(BAR03) == 0 ) ngaub(BAR03) = 3
           else if(ngaus<=10) then
              if( ngaub(BAR03) == 0 ) ngaub(BAR03) = 4
           end if
        end if

     else if( jelty == QUA04 ) then   

        lexib(BAR02) = 1
        lquab(BAR02) = lquad
        if( ngaub(BAR02) == 0 ) ngaub(BAR02) = int(sqrt(real(ngaus,rp)))

     else if( jelty == QUA08 ) then   

        lexib(BAR03) = 1
        lquab(BAR03) = lquad
        if( ngaub(BAR03) == 0 ) ngaub(BAR03) = int(sqrt(real(ngaus,rp)))

     else if( jelty == QUA09 ) then   

        lexib(BAR03) = 1
        lquab(BAR03) = lquad
        if( ngaub(BAR03) == 0 ) ngaub(BAR03) = int(sqrt(real(ngaus,rp)))

     else if( jelty == QUA16 ) then   

        lexib(BAR04) = 1
        lquab(BAR04) = lquad
        if( ngaub(BAR04) == 0 ) ngaub(BAR04) = int(sqrt(real(ngaus,rp)))


     else if( jelty == TET04 ) then  

        lexib(TRI03) = 1
        lquab(TRI03) = lquad

        if( lquad == 0 ) then                  ! open rule
           if(ngaus==1) then
              if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 1                     ! exact for P1
           else if(ngaus<=5 ) then   
              if( ngaub(TRI03) == 0 ) ngaub(TRI03) = ngaus-1               ! exact for P2/P3
           else if(ngaus==11) then
              if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 6                     ! exact for P4
           else if(ngaus==14) then
              if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 7                     ! exact for P5
           end if
        else if( lquad == 1 ) then             ! closed rule
           if(ngaus<=5) then
              if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 3
           else if(ngaus<=11) then
              if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 6
           else if(ngaus==15) then
              if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 7
           else if(ngaus==20) then
              if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 10
           end if
        end if

     else if( jelty == TET10 ) then 

        lexib(TRI06) = 1
        lquab(TRI06) = lquad

        if( lquad == 0 ) then                  ! open rule
           if(ngaus==1) then
              if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 1                     ! exact for P1
           else if(ngaus<=5 ) then   
              if( ngaub(TRI06) == 0 ) ngaub(TRI06) = ngaus-1               ! exact for P2/P3
           else if(ngaus==8 ) then
              if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 6                     ! exact for P4
           else if(ngaus==11) then
              if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 6                     ! exact for P4
           else if(ngaus==14) then
              if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 7                     ! exact for P5
           else if(ngaus==15) then
              if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 7                     ! exact for P5
           else if(ngaus==45) then
              if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 7                     ! exact for P5
           else if(ngaus==65) then
              if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 7                     ! exact for P5
           end if
        else if( lquad == 1 ) then             ! closed rule
           if(ngaus<=5) then
              if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 3
           else if(ngaus<=11) then
              if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 6
           else if(ngaus==15) then
              if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 7
           else if(ngaus==20) then
              if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 10
           end if
        end if

     else if( jelty == PYR05) then

        lexib(TRI03) = 1
        lexib(QUA04) = 1
        if(ngaus==1) then
           if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 1
           if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 1
        else if(ngaus==5.or.ngaus==6.or.ngaus==8.or.ngaus==9) then
           if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 3
           if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 4        
        else if(ngaus==13.or.ngaus==18) then
           if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 6
           if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 9                  
        else
           call runend('BOUELE: WRONG NUMBER OF GAUSS POINTS FOR PYRA_5 ')
        end if

     else if( jelty == PYR14 ) then  

        lexib(TRI06) = 1
        lexib(QUA08) = 1
        call runend('BOUELE: PYRA_14 ELEMENT IS NOT READY')

     else if( jelty == PEN06 ) then  
        !
        ! PEN06
        !
        lexib(TRI03) = 1
        lexib(QUA04) = 1
        lquab(TRI03) = lquad
        lquab(QUA04) = lquad

        if(ngaus==1) then
           if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 1
           if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 1
        else if(ngaus==6) then
           if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 3
           if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 4
        else if(ngaus==8) then
           if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 4
           if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 4
        else if(ngaus==11) then
           if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 6
           if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 9
        else if(ngaus==16) then
           if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 7
           if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 9
        else if(ngaus==24) then
           if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 13
           if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 16
        else if(ngaus==29) then
           if( ngaub(TRI03) == 0 ) ngaub(TRI03) = 13
           if( ngaub(QUA04) == 0 ) ngaub(QUA04) = 16
        else
           call runend('PEN06: NOT CODED')
        end if

     else if( jelty == PEN15 ) then
        !
        ! PEN15
        !
        lexib(TRI06 ) = 1
        lexib(QUA08) = 1

        if(ngaus==1) then
           if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 1
           if( ngaub(QUA08) == 0 ) ngaub(QUA08) = 1
        else if(ngaus==6) then
           if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 3
           if( ngaub(QUA08) == 0 ) ngaub(QUA08) = 4
        end if

        call runend('BOUELE: PENTA_15 ELEMENT IS NOT READY')

     else if( jelty == PEN18 ) then
        !
        ! PEN18
        !
        lexib(TRI06) = 1
        lexib(QUA09) = 1

        if(ngaus==1) then
           if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 1
           if( ngaub(QUA09) == 0 ) ngaub(QUA09) = 1
        else if(ngaus==6) then
           if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 3
           if( ngaub(QUA09) == 0 ) ngaub(QUA09) = 4
        else if(ngaus==18) then
           if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 6
           if( ngaub(QUA09) == 0 ) ngaub(QUA09) = 9
        else if(ngaus==21) then
           if( ngaub(TRI06) == 0 ) ngaub(TRI06) = 7
           if( ngaub(QUA09) == 0 ) ngaub(QUA09) = 9
        end if

        ! call runend('BOUELE: PENTA_18 ELEMENT IS NOT READY')

     else if( jelty == HEX08 ) then   
        !
        ! HEX08
        !
        lexib(QUA04) = 1
        lquab(QUA04) = lquad
        if( ngaub(QUA04) == 0 ) ngaub(QUA04) = nint(real(ngaus,rp)**(2.0_rp/3.0_rp))

     else if( jelty == HEX20 ) then  

        lexib(QUA08) = 1
        lquab(QUA08) = lquad
        if(ngaus==20) then
           if( ngaub(QUA08) == 0 ) ngaub(QUA08) = 9
        else
           if( ngaub(QUA08) == 0 ) ngaub(QUA08) = nint(real(ngaus,rp)**(2.0_rp/3.0_rp))
        end if

     else if( jelty == HEX27 ) then 

        lexib(QUA09) = 1
        lquab(QUA09) = lquad
        if( ngaub(QUA09) == 0 ) ngaub(QUA09) = nint(real(ngaus,rp)**(2.0_rp/3.0_rp))

     else if( jelty == SHELL ) then   

        lexib(BAR02) = 1
        lquab(BAR02) = lquad
        if( ngaub(BAR02) == 0 ) ngaub(BAR02) = 1

     else if( jelty == BAR3D ) then   

        lexib(POINT) = 1
        lquab(POINT) = lquad
        ngaub(POINT) = 1

     end if

  end if

end subroutine bouele

