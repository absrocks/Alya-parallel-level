complex(rp) function ARMONICO(X,Y,Z,L,M)
  use def_kintyp, only    :  ip,rp
  use def_parame, only    :  PI
  implicit none
  complex(rp)             :: D
  real(rp)                :: TITA, PHI, ANCTE, ANLEG, DENOM, rta,funleg
  INTEGER(ip), intent(in) :: L,M
  real(rp),    intent(in) :: X,Y,Z
  integer(ip)             :: NFAC,M_M



  M_M=(-1)*M
  DENOM = SQRT(X**2 + Y**2 + Z**2)
  IF( DENOM /= 0.0_rp ) THEN
     TITA = Z/DENOM  ! COS(TITA)    
  ELSE
     TITA = 0.0_rp
  ENDIF

  IF( X /= 0.0_rp ) THEN 
     RTA = ATAN(Y/X)
     IF( X > 0.0_rp .AND. Y >= 0.0_rp) PHI = RTA
     IF( X < 0.0_rp .AND. Y >= 0.0_rp) PHI = PI        + RTA
     IF( X < 0.0_rp .AND. Y <  0.0_rp) PHI = PI        + RTA
     IF( X > 0.0_rp .AND. Y <  0.0_rp) PHI = 2.0_rp*PI + RTA
  ELSE
     IF( Y >= 0.0_rp ) PHI = PI*0.5_rp
     IF( Y <  0.0_rp ) PHI = 3.0_rp*PI*0.5_rp
  ENDIF

  ANCTE = (2.0_rp*real(L)+1.0_rp)/(4.0_rp*PI) 
  ANCTE = SQRT( ANCTE * real(NFAC(L-M)) / real(NFAC(L+M)) )

  ANLEG = 1.0_rp
  IF( M  <  0 ) THEN
     ANLEG = (-1.0_rp)**(M_M) * real(NFAC(L-M_M)) * 1.0_rp / real(NFAC(L+M_M))
  ENDIF

  D = (0.0_rp,1.0_rp)

  ARMONICO = ANCTE*ANLEG*FUNLEG(TITA,L,ABS(M))*EXP(D*real(M)*PHI) 

end function armonico


real(rp) FUNCTION FUNYLM(X,Y,Z,L)
  use def_master
  implicit none
  real(rp), intent(in) :: X,Y,Z
  real(rp)             :: TITA,DENOM
  integer(ip)          :: L
  integer(ip)          :: M,kk,NCTE,NCTEl,NLEG,NFAC,FUNLEG

  FUNYLM=0.0

  DENOM= SQRT(X**2 + Y**2 + Z**2)
  IF(DENOM.NE.0) THEN
     TITA = Z/DENOM   ! COS(TITA)
  ELSE
     TITA = 0.0_rp
  ENDIF

  NCTEL =  (2*L+1)/(16*ATAN(1.0)) 
  DO KK=1,2*L+1
     M=-L + KK - 1  
     NCTE = NCTEL * NFAC(L-M)/NFAC(L+M)
     NLEG = 1
     IF(M < 0) NLEG = (-1)**M *  NFAC(L-M)/NFAC(L+M)
     FUNYLM = FUNYLM + NCTE * (NLEG*FUNLEG(TITA,L,ABS(M)))**2  
  ENDDO


end function FUNYLM


real(rp) FUNCTION FUNLEG(X,N,M)
  !C
  !C       ========================================================
  !C       Purpose: Compute associated Legendre functions Pmn(x)
  !C                and Pmn'(x) for a given order
  !C       Input :  x --- Argument of Pmn(x)
  !C                m --- Order of Pmn(x),  m = 0,1,2,...,n
  !C                n --- Degree of Pmn(x), n = 0,1,2,...,N
  !C       Output:  PM(n) --- Pmn(x)
  !C                PD(n) --- Pmn'(x)
  !C       ========================================================

  use def_master
  implicit none
  real(rp),   intent(in) :: x
  integer(ip),intent(in) :: N,M
  real(rp)               :: PMK,pm0,pm1,pm2,x0        
  integer(ip)            :: k  

  FUNLEG = 0.0_rp

  IF ( ABS(X) == 1.0_rp ) THEN
     IF ( M == 0 ) THEN
        FUNLEG = 1.0_rp
        IF ( X < 0.0_rp ) THEN
           FUNLEG = (-1.0_rp)**N*FUNLEG
        ENDIF
     ENDIF
     RETURN
  ENDIF
  X0=DABS(1.0D0-X*X)
  PM0=1.0D0
  PMK=PM0
  DO K=1,M
     PMK=(2.0D0*K-1.0D0)*SQRT(X0)*PM0
     PM0=PMK
  ENDDO
  PM1=(2.0D0*M+1.0D0)*X*PM0

  IF(N.EQ.M) THEN 
     FUNLEG= PMK
  ELSEIF(N.EQ.M+1) THEN
     FUNLEG=PM1
  ELSE

     DO K=M+2,N
        PM2=((2.0D0*K-1.0D0)*X*PM1-(K+M-1.0D0)*PMK)/(K-M)
        FUNLEG=PM2
        PMK=PM1
        PM1=PM2
     ENDDO

  ENDIF

end function FUNLEG


INTEGER(ip) FUNCTION NFAC(N)
  use def_master
  implicit none
  integer(ip),intent(in) :: N
  integer(ip) :: kk

  NFAC = 1

  DO KK=1,N
     NFAC=NFAC*KK
  ENDDO

end function nfac


complex(rp) FUNCTION FUNTEOR(NCP,LL,MM,ZZ,X,Y,Z,X0,Y0,Z0)
  use def_master
  implicit none
  integer(ip),intent(in) :: NCP,LL,MM
  real(rp), intent(in)   :: ZZ,X,Y,Z,X0,Y0,Z0
  INTEGER(ip)            :: N1,N2,NFAC
  real(rp)               :: AN1,AN2,ACTE,ACTE2,RAD,HG,ALFA,FHGM

  RAD = SQRT((X-X0)**2 + (Y-Y0)**2 + (Z-Z0)**2)


  ALFA = 2.0_rp*ZZ/real(NCP)

  ACTE= ALFA**3 * real(NFAC(NCP-LL-1))/(2.0_rp*real(NCP)*real(NFAC(NCP+LL)))
  ACTE= SQRT(ACTE)

  ACTE2 = real(NFAC(NCP))/real(NFAC(NCP-LL-1))

  RAD = ALFA*RAD

  N1  = NCP-LL-1
  N2  = 2*LL+1
  AN1 = real(N1)
  AN2 = real(N2)

  HG =  real(NFAC(N1+N2))/real(NFAC(N2)) * FHGM(-AN1,AN2+1.0_rp,RAD)/real(NFAC(NCP))


  FUNTEOR = ACTE * EXP(-RAD*0.5_rp)*(RAD)**LL* HG  * ACTE2 

end function funteor



complex(rp) FUNCTION FUNTEOR1d(NCP,LL,MM,ZZ,X,X0)
  use def_master
  implicit none
  integer(ip),intent(in) ::NCP,LL,MM
  real(rp), intent(in) :: ZZ,X,X0
  INTEGER(ip)  :: N1,N2,NFAC
  real(rp)     :: AN1,AN2,ACTE,ACTE2,RAD,HG,ALFA,FHGM

  RAD = SQRT((X-X0)**2 )


  ALFA = 2*ZZ/NCP

  ACTE= ALFA**3 * NFAC(NCP-LL-1)/(2*NCP*NFAC(NCP+LL))
  ACTE= SQRT(ACTE)

  ACTE2 = NFAC(NCP)/NFAC(NCP-LL-1)

  RAD=ALFA*RAD

  N1=NCP-LL-1
  N2=2*LL+1
  AN1=N1
  AN2=N2

  HG =  NFAC(N1+N2)/NFAC(N2) * FHGM(-AN1,AN2+1,RAD)/NFAC(NCP)


  FUNTEOR1D = ACTE* EXP(-RAD*0.5)*(RAD)**LL* HG  * ACTE2 

end function funteor1d


real(rp) FUNCTION FHGM(A,B,X)
  use def_parame, only :  PI
  use def_master
  implicit none
  real(rp)             :: A,B,X
  INTEGER(ip)          :: N1,N2,m,nl,la,j,k,N,i
  real(rp)             :: A0,a1,x0,r,RG,xg,Y0,Y1,ta,tb,TBA,sum1,sum2,r1,r2,FHGM1,FHGM2
  real(rp)             :: ii,jj

  A0 = A
  A1 = A
  X0 = X
  FHGM=0.0_rp
  IF (B.EQ.0.0_rp.OR.B.EQ.-ABS(INT(B))) THEN
     FHGM=1.0D+300
  ELSE IF (A.EQ.0.0_rp.OR.X.EQ.0.0_rp) THEN
     FHGM=1.0_rp
  ELSE IF (A.EQ.-1.0_rp) THEN
     FHGM=1.0_rp-X/B
  ELSE IF (A.EQ.B) THEN
     FHGM=EXP(X)
  ELSE IF (A-B.EQ.1.0_rp) THEN
     FHGM=(1.0_rp+X/B)*EXP(X)
  ELSE IF (A.EQ.1.0_rp.AND.B.EQ.2.0_rp) THEN
     FHGM=(DEXP(X)-1.0_rp)/X
  ELSE IF (A.EQ.INT(A).AND.A < 0.0_rp) THEN
     M=INT(-A)
     R=1.0_rp
     FHGM=1.0_rp
     DO K=1,M
        R=R*(A+K-1.0_rp)/K/(B+K-1.0_rp)*X
        FHGM=FHGM+R
     enddo
  ENDIF
  IF (FHGM.NE.0.0_rp) RETURN
  IF (X < 0.0_rp) THEN
     A=B-A
     A0=A
     X=ABS(X)
  ENDIF
  IF (A < 2.0_rp) NL=0
  IF (A >= 2.0_rp) THEN
     NL=1
     LA=INT(A)
     A=A-LA-1.0_rp
  ENDIF
  DO N=0,NL
     IF (A0 >= 2.0_rp) A=A+1.0_rp
     IF (X.LE.30.0_rp+ABS(B).OR.A < 0.0_rp) THEN
        FHGM=1.0_rp
        RG=1.0_rp
        DO J=1,500
           JJ = REAL(J)
           RG=RG*(A+JJ-1.0_rp)/(JJ*(B+JJ-1.0_rp))*X
           FHGM=FHGM+RG
           IF (ABS(RG/FHGM) < 1.0D-15) GO TO 25
        enddo
     ELSE
        CALL GAMMA(A,TA)
        CALL GAMMA(B,TB)
        XG=B-A
        CALL GAMMA(XG,TBA)
        SUM1=1.0_rp
        SUM2=1.0_rp
        R1=1.0_rp
        R2=1.0_rp
        DO I=1,8
           II = REAL(I)
           R1=-R1*(A + II - 1.0_rp)*(A-B+II)/(X*II)
           R2=-R2*(B - A + II -1.0_rp)*(A-II)/(X*II)
           SUM1=SUM1+R1
           SUM2=SUM2+R2
        enddo
        FHGM1=TB/TBA*X**(-A)*COS(PI*A)*SUM1
        FHGM2=TB/TA*EXP(X)*X**(A-B)*SUM2
        FHGM=FHGM1+FHGM2
     ENDIF
25   IF (N.EQ.0) Y0=FHGM
     IF (N.EQ.1) Y1=FHGM
  enddo

  IF (A0 >= 2.0_rp) THEN
     DO  I=1,LA-1
        FHGM=((2.0_rp*A-B+X)*Y1+(B-A)*Y0)/A
        Y0=Y1
        Y1=FHGM
        A=A+1.0_rp
     enddo
  ENDIF
  IF (X0 < 0.0_rp) FHGM=FHGM*DEXP(X0)
  A=A1
  X=X0

end function fhgm


SUBROUTINE GAMMA(X,GA)
  use def_master
  implicit none
  real(rp), intent(in) :: X
  real(rp)             :: G(26)
  real(rp)          :: PI,z,r,GA, GR
  integer(ip)       :: k,m1,m

  DATA G/1.0_rp,0.5772156649015329_rp, -0.6558780715202538_rp, -0.420026350340952D-1,           &
       & 0.1665386113822915_rp,-.421977345555443D-1,-.96219715278770D-2, .72189432466630D-2,  &
       & -.11651675918591D-2, -.2152416741149D-3, .1280502823882D-3, -.201348547807D-4,      &
       & -.12504934821D-5, .11330272320D-5,-.2056338417D-6, .61160950D-8,                    &
       & .50020075D-8, -.11812746D-8,.1043427D-9, .77823D-11,-.36968D-11, .51D-12,           &
       & -.206D-13, -.54D-14, .14D-14, .1D-15/         


  PI=3.141592653589793_rp
  IF (X.EQ.INT(X)) THEN
     IF (X > 0.0_rp) THEN
        GA=1.0_rp
        M1=X-1
        DO  K=2,M1
           GA=GA*K
        enddo
     ELSE
        GA=1.0D+300
     ENDIF
  ELSE
     IF (DABS(X) > 1.0_rp) THEN
        Z=DABS(X)
        M=INT(Z)
        R=1.0_rp
        DO K=1,M
           R=R*(Z-K)
        enddo
        Z=Z-M
     ELSE
        Z=X
     ENDIF
     GR=G(26)
     DO K=25,1,-1
        GR=GR*Z+G(K)
     enddo
     GA=1.0_rp/(GR*Z)
     IF (DABS(X) > 1.0_rp) THEN
        GA=GA*R
        IF (X < 0.0_rp) GA=-PI/(X*GA*DSIN(PI*X))
     ENDIF
  ENDIF

end subroutine GAMMA
