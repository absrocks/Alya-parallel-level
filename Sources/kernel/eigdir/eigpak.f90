subroutine QZHES(N,A,B,matz,Z)

  use      def_solver

  implicit none

  integer(ip)             :: n
  real(rp)                :: a(n,n),b(n,n),z(n,n)
  integer(ip)             :: i,j,l,nm1,l1,nk1,nm2,k,lb
  real(rp)                :: s,r,rho,t,v1,v2,u2,u1
  logical                 :: matz


!C     :::::::::: INITIALIZE Z ::::::::::                                !QZH00520
      IF (.NOT. MATZ) GO TO 10                                          !!QZH00530

      DO 3 I = 1, n                                                     !!QZH00550

         DO 2 J = 1, n                                                  !!QZH00570
            Z(I,J) = 0.0D0                                              !!QZH00580
    2    CONTINUE                                                       !!QZH00590

         Z(I,I) = 1.0D0                                                 !!QZH00610
    3 CONTINUE                                                          !!QZH00620
!C     :::::::::: REDUCE B TO UPPER TRIANGULAR FORM ::::::::::           !QZH00630
   10 IF (N .LE. 1) GO TO 170                                           !!QZH00640
      NM1 = n - 1                                                       !!QZH00650

      DO 100 L = 1, NM1                                                 !!QZH00670
         L1 = L + 1                                                     !!QZH00680
         S = 0.0D0                                                      !!QZH00690

         DO 20 I = L1, N                                                !!QZH00710
            S = S + ABS(B(I,L))                                        !!QZH00720
   20    CONTINUE                                                       !!QZH00730

         IF (S .EQ. 0.0D0) GO TO 100                                    !!QZH00750
         S = S + ABS(B(L,L))                                           !!QZH00760
         R = 0.0D0                                                      !!QZH00770

         DO 25 I = L, N                                                 !!QZH00790
            B(I,L) = B(I,L) / S                                         !!QZH00800
            R = R + B(I,L)**2                                           !!QZH00810
   25    CONTINUE                                                       !!QZH00820

         R = SIGN(SQRT(R),B(L,L))                                     !QZH00840
         B(L,L) = B(L,L) + R                                            !QZH00850
         RHO = R * B(L,L)                                               !QZH00860

         DO 50 J = L1, N                                                !QZH00880
            T = 0.0D0                                                   !QZH00890

            DO 30 I = L, N                                              !QZH00910
               T = T + B(I,L) * B(I,J)                                  !QZH00920
   30       CONTINUE                                                    !QZH00930

            T = -T / RHO                                                !QZH00950

            DO 40 I = L, N                                              !QZH00970
               B(I,J) = B(I,J) + T * B(I,L)                             !QZH00980
   40       CONTINUE                                                    !QZH00990

   50    CONTINUE                                                       !QZH01010

         DO 80 J = 1, N                                                 !QZH01030
            T = 0.0D0                                                   !QZH01040

            DO 60 I = L, N                                              !QZH01060
               T = T + B(I,L) * A(I,J)                                  !QZH01070
   60       CONTINUE                                                    !QZH01080

            T = -T / RHO                                                !QZH01100

            DO 70 I = L, N                                              !QZH01120
               A(I,J) = A(I,J) + T * B(I,L)                             !QZH01130
   70       CONTINUE                                                    !QZH01140

   80    CONTINUE                                                       !QZH01160

         B(L,L) = -S * R                                                !QZH01180

         DO 90 I = L1, N                                                !QZH01200
            B(I,L) = 0.0D0                                              !QZH01210
   90    CONTINUE                                                       !QZH01220

  100 CONTINUE                                                          !QZH01240
!C     :::::::::: REDUCE A TO UPPER HESSENBERG FORM, WHILE               !QZH01250
!C                KEEPING B TRIANGULAR ::::::::::                        !QZH01260
      IF (N .EQ. 2) GO TO 170                                           !QZH01270
      NM2 = N - 2                                                       !QZH01280

      DO 160 K = 1, NM2                                                 !QZH01300
         NK1 = NM1 - K                                                  !QZH01310
!C     :::::::::: FOR L=N-1 STEP -1 UNTIL K+1 DO -- ::::::::::           !QZH01320
         DO 150 LB = 1, NK1                                             !QZH01330
            L = N - LB                                                  !QZH01340
            L1 = L + 1                                                  !QZH01350
!C     :::::::::: ZERO A(L+1,K) ::::::::::                               !QZH01360
            S = ABS(A(L,K)) + ABS(A(L1,K))                            !QZH01370
            IF (S .EQ. 0.0D0) GO TO 150                                 !QZH01380
            U1 = A(L,K) / S                                             !QZH01390
            U2 = A(L1,K) / S                                            !QZH01400
            R = SIGN(SQRT(U1*U1+U2*U2),U1)                            !QZH01410
            V1 =  -(U1 + R) / R                                         !QZH01420
            V2 = -U2 / R                                                !QZH01430
            U2 = V2 / V1                                                !QZH01440

            DO 110 J = K, N                                             !QZH01460
               T = A(L,J) + U2 * A(L1,J)                                !QZH01470
               A(L,J) = A(L,J) + T * V1                                 !QZH01480
               A(L1,J) = A(L1,J) + T * V2                               !QZH01490
  110       CONTINUE                                                    !QZH01500

            A(L1,K) = 0.0D0                                             !QZH01520

            DO 120 J = L, N                                             !QZH01540
               T = B(L,J) + U2 * B(L1,J)                                !QZH01550
               B(L,J) = B(L,J) + T * V1                                 !QZH01560
               B(L1,J) = B(L1,J) + T * V2                               !QZH01570
  120       CONTINUE                                                    !QZH01580
!C     :::::::::: ZERO B(L+1,L) ::::::::::                               !QZH01590
            S = ABS(B(L1,L1)) + ABS(B(L1,L))                          !QZH01600
            IF (S .EQ. 0.0D0) GO TO 150                                 !QZH01610
            U1 = B(L1,L1) / S                                           !QZH01620
            U2 = B(L1,L) / S                                            !QZH01630
            R = SIGN(SQRT(U1*U1+U2*U2),U1)                            !QZH01640
            V1 =  -(U1 + R) / R                                         !QZH01650
            V2 = -U2 / R                                                !QZH01660
            U2 = V2 / V1                                                !QZH01670

            DO 130 I = 1, L1                                            !QZH01690
               T = B(I,L1) + U2 * B(I,L)                                !QZH01700
               B(I,L1) = B(I,L1) + T * V1                               !QZH01710
               B(I,L) = B(I,L) + T * V2                                 !QZH01720
  130       CONTINUE                                                    !QZH01730

            B(L1,L) = 0.0D0                                             !QZH01750

            DO 140 I = 1, N                                             !QZH01770
               T = A(I,L1) + U2 * A(I,L)                                !QZH01780
               A(I,L1) = A(I,L1) + T * V1                               !QZH01790
               A(I,L) = A(I,L) + T * V2                                 !QZH01800
  140       CONTINUE                                                    !QZH01810

            IF (.NOT. MATZ) GO TO 150                                   !QZH01830

            DO 145 I = 1, N                                             !QZH01850
               T = Z(I,L1) + U2 * Z(I,L)                                !QZH01860
               Z(I,L1) = Z(I,L1) + T * V1                               !QZH01870
               Z(I,L) = Z(I,L) + T * V2                                 !QZH01880
  145       CONTINUE                                                    !QZH01890

  150    CONTINUE                                                       !QZH01910

  160 CONTINUE                                                          !QZH01930

  170 RETURN                                                            !QZH01950
!C     :::::::::: LAST CARD OF QZHES ::::::::::                          !QZH01960
      END
      

subroutine QZIT(N,A,B,eps1,matz,Z,ierr,niter)

  use      def_solver

  implicit none

  integer(ip)             :: n
  real(rp)     :: a(n,n),b(n,n),z(n,n)
  integer(ip)  :: I,J,K,L,K1,K2,LD,LL,L1,NA,NM,ISH,ITS,KM1,LM1,   & 
                  NENM2,IERR,LOR1,NENORN,NITER,nen
  real(rp)     :: R,S,T,A1,A2,A3,EP,SH,U1,U2,U3,V1,V2,V3,ANI,A11,A12,&        
                  A21,A22,A33,A34,A43,A44,BNI,B11,B12,B22,B33,B34,   &       
                  B44,EPSA,EPSB,EPS1,ANORM,BNORM
  logical      :: matz,notlas




!C        IERR IS SET TO                                                 !QZI00700
!C          ZERO       FOR NORMAL RETURN,                                !QZI00710
!C          J          IF NEITHER A(J,J-1) NOR A(J-1,J-2) HAS BECOME     !QZI00720
!C                     ZERO AFTER 500 ITERATIONS.                         !QZI00730

      IERR = 0                                                          !QZI00800
!C     :::::::::: COMPUTE EPSA,EPSB ::::::::::                           !QZI00810
      ANORM = 0.0D0                                                     !QZI00820
      BNORM = 0.0D0                                                     !QZI00830

      DO 30 I = 1, N                                                    !QZI00850
         ANI = 0.0D0                                                    !QZI00860
         IF (I .NE. 1) ANI = ABS(A(I,I-1))                             !QZI00870
         BNI = 0.0D0                                                    !QZI00880

         DO 20 J = I, N                                                 !QZI00900
            ANI = ANI + ABS(A(I,J))                                    !QZI00910
            BNI = BNI + ABS(B(I,J))                                    !QZI00920
   20    CONTINUE                                                       !QZI00930

         IF (ANI .GT. ANORM) ANORM = ANI                                !QZI00950
         IF (BNI .GT. BNORM) BNORM = BNI                                !QZI00960
   30 CONTINUE                                                          !QZI00970

      IF (ANORM .EQ. 0.0D0) ANORM = 1.0D0                               !QZI00990
      IF (BNORM .EQ. 0.0D0) BNORM = 1.0D0                               !QZI01000
      EP = EPS1                                                         !QZI01010
      IF (EP .GT. 0.0D0) GO TO 50                                       !QZI01020
!C     :::::::::: COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO ::::::::::      !QZI01030
      EP = 1.0D0                                                        !QZI01040
   40 EP = EP / 2.0D0                                                   !QZI01050
      IF (1.0D0 + EP .GT. 1.0D0) GO TO 40                               !QZI01060
   50 EPSA = EP * ANORM                                                 !QZI01070
      EPSB = EP * BNORM                                                 !QZI01080
!C     :::::::::: REDUCE A TO QUASI-TRIANGULAR FORM, WHILE               !QZI01090
!C                KEEPING B TRIANGULAR ::::::::::                        !QZI01100
      LOR1 = 1                                                          !QZI01110
      NENORN = N                                                         !QZI01120
      NEN = N                                                            !QZI01130
!C     :::::::::: BEGIN QZ STEP ::::::::::                               !QZI01140
   60 IF (NEN .LE. 2) GO TO 1001                                         !QZI01150
      IF (.NOT. MATZ) NENORN = NEN                                        !QZI01160
      ITS = 0                                                           !QZI01170
      NA = NEN - 1                                                       !QZI01180
      NENM2 = NA - 1                                                     !QZI01190
   70 ISH = 2                                                           !QZI01200
!C     :::::::::: CHECK FOR CONVERGENCE OR REDUCIBILITY.                 !QZI01210
!C                FOR L=NEN STEP -1 UNTIL 1 DO -- ::::::::::              !QZI01220
      DO 80 LL = 1, NEN                                                  !QZI01230
         LM1 = NEN - LL                                                  !QZI01240
         L = LM1 + 1                                                    !QZI01250
         IF (L .EQ. 1) GO TO 95                                         !QZI01260
         IF (ABS(A(L,LM1)) .LE. EPSA) GO TO 90                         !QZI01270
   80 CONTINUE                                                          !QZI01280

   90 A(L,LM1) = 0.0D0                                                  !QZI01300
      IF (L .LT. NA) GO TO 95                                           !QZI01310
!C     :::::::::: 1-BY-1 OR 2-BY-2 BLOCK ISOLATED ::::::::::             !QZI01320
      NEN = LM1                                                          !QZI01330
      GO TO 60                                                          !QZI01340
!C     :::::::::: CHECK FOR SMALL TOP OF B ::::::::::                    !QZI01350
   95 LD = L                                                            !QZI01360
  100 L1 = L + 1                                                        !QZI01370
      B11 = B(L,L)                                                      !QZI01380
      IF (ABS(B11) .GT. EPSB) GO TO 120                                !QZI01390
      B(L,L) = 0.0D0                                                    !QZI01400
      S = ABS(A(L,L)) + ABS(A(L1,L))                                  !QZI01410
      U1 = A(L,L) / S                                                   !QZI01420
      U2 = A(L1,L) / S                                                  !QZI01430
      R = SIGN(SQRT(U1*U1+U2*U2),U1)                                  !QZI01440
      V1 = -(U1 + R) / R                                                !QZI01450
      V2 = -U2 / R                                                      !QZI01460
      U2 = V2 / V1                                                      !QZI01470

      DO 110 J = L, NENORN                                               !QZI01490
         T = A(L,J) + U2 * A(L1,J)                                      !QZI01500
         A(L,J) = A(L,J) + T * V1                                       !QZI01510
         A(L1,J) = A(L1,J) + T * V2                                     !QZI01520
         T = B(L,J) + U2 * B(L1,J)                                      !QZI01530
         B(L,J) = B(L,J) + T * V1                                       !QZI01540
         B(L1,J) = B(L1,J) + T * V2                                     !QZI01550
  110 CONTINUE                                                          !QZI01560

      IF (L .NE. 1) A(L,LM1) = -A(L,LM1)                                !QZI01580
      LM1 = L                                                           !QZI01590
      L = L1                                                            !QZI01600
      GO TO 90                                                          !QZI01610
  120 A11 = A(L,L) / B11                                                !QZI01620
      A21 = A(L1,L) / B11                                               !QZI01630
      IF (ISH .EQ. 1) GO TO 140                                         !QZI01640
!C     :::::::::: ITERATION STRATEGY ::::::::::                          !QZI01650
      IF (ITS .EQ. niter) GO TO 1000                                    !QZI01660
      IF (ITS .EQ. 10) GO TO 155                                        !QZI01670
!C     :::::::::: DETERMINE TYPE OF SHIFT ::::::::::                     !QZI01680
      B22 = B(L1,L1)                                                    !QZI01690
      IF (ABS(B22) .LT. EPSB) B22 = EPSB                               !QZI01700
      B33 = B(NA,NA)                                                    !QZI01710
      IF (ABS(B33) .LT. EPSB) B33 = EPSB                               !QZI01720
      B44 = B(NEN,NEN)                                                    !QZI01730
      IF (ABS(B44) .LT. EPSB) B44 = EPSB                               !QZI01740
      A33 = A(NA,NA) / B33                                              !QZI01750
      A34 = A(NA,NEN) / B44                                              !QZI01760
      A43 = A(NEN,NA) / B33                                              !QZI01770
      A44 = A(NEN,NEN) / B44                                              !QZI01780
      B34 = B(NA,NEN) / B44                                              !QZI01790
      T = 0.5D0 * (A43 * B34 - A33 - A44)                               !QZI01800
      R = T * T + A34 * A43 - A33 * A44                                 !QZI01810
      IF (R .LT. 0.0D0) GO TO 150                                       !QZI01820
!C     :::::::::: DETERMINE SINGLE SHIFT ZEROTH COLUMN OF A ::::::::::   !QZI01830
      ISH = 1                                                           !QZI01840
      R = SQRT(R)                                                      !QZI01850
      SH = -T + R                                                       !QZI01860
      S = -T - R                                                        !QZI01870
      IF (ABS(S-A44) .LT. ABS(SH-A44)) SH = S                         !QZI01880
!C     :::::::::: LOOK FOR TWO CONSECUTIVE SMALL                         !QZI01890
!C                SUB-DIAGONAL ELEMENTS OF A.                            !QZI01900
!C                FOR L=NEN-2 STEP -1 UNTIL LD DO -- ::::::::::           !QZI01910
      DO 130 LL = LD, NENM2                                              !QZI01920
         L = NENM2 + LD - LL                                             !QZI01930
         IF (L .EQ. LD) GO TO 140                                       !QZI01940
         LM1 = L - 1                                                    !QZI01950
         L1 = L + 1                                                     !QZI01960
         T = A(L,L)                                                     !QZI01970
         IF (ABS(B(L,L)) .GT. EPSB) T = T - SH * B(L,L)                !QZI01980
         IF (ABS(A(L,LM1)) .LE. ABS(T/A(L1,L)) * EPSA) GO TO 100      !QZI01990
  130 CONTINUE                                                          !QZI02000

  140 A1 = A11 - SH                                                     !QZI02020
      A2 = A21                                                          !QZI02030
      IF (L .NE. LD) A(L,LM1) = -A(L,LM1)                               !QZI02040
      GO TO 160                                                         !QZI02050
!C     :::::::::: DETERMINE DOUBLE SHIFT ZEROTH COLUMN OF A ::::::::::   !QZI02060
  150 A12 = A(L,L1) / B22                                               !QZI02070
      A22 = A(L1,L1) / B22                                              !QZI02080
      B12 = B(L,L1) / B22                                               !QZI02090
      A1 = ((A33 - A11) * (A44 - A11) - A34 * A43 + A43 * B34 * A11) &  !QZI02100
          / A21 + A12 - A11 * B12                                      !QZI02110
      A2 = (A22 - A11) - A21 * B12 - (A33 - A11) - (A44 - A11)       & !QZI02120
          + A43 * B34                                                  !QZI02130
      A3 = A(L1+1,L1) / B22                                             !QZI02140
      GO TO 160                                                         !QZI02150
!C     :::::::::: AD HOC SHIFT ::::::::::                                !QZI02160
  155 A1 = 0.0D0                                                        !QZI02170
      A2 = 1.0D0                                                        !QZI02180
      A3 = 1.1605D0                                                     !QZI02190
  160 ITS = ITS + 1                                                     !QZI02200
      IF (.NOT. MATZ) LOR1 = LD                                         !QZI02210
!C     :::::::::: MAIN LOOP ::::::::::                                   !QZI02220
      DO 260 K = L, NA                                                  !QZI02230
         NOTLAS = K .NE. NA .AND. ISH .EQ. 2                            !QZI02240
         K1 = K + 1                                                     !QZI02250
         K2 = K + 2                                                     !QZI02260
         KM1 = MAX(K-1_ip,L)                                              !QZI02270
         LL = MIN(NEN,K1+ISH)                                           !QZI02280
         IF (NOTLAS) GO TO 190                                          !QZI02290
!C     :::::::::: ZERO A(K+1,K-1) ::::::::::                             !QZI02300
         IF (K .EQ. L) GO TO 170                                        !QZI02310
         A1 = A(K,KM1)                                                  !QZI02320
         A2 = A(K1,KM1)                                                 !QZI02330
  170    S = ABS(A1) + ABS(A2)                                        !QZI02340
         IF (S .EQ. 0.0D0) GO TO 190                                    !QZI02350
         U1 = A1 / S                                                    !QZI02360
         U2 = A2 / S                                                    !QZI02370
         R = SIGN(SQRT(U1*U1+U2*U2),U1)                               !QZI02380
         V1 = -(U1 + R) / R                                             !QZI02390
         V2 = -U2 / R                                                   !QZI02400
         U2 = V2 / V1                                                   !QZI02410

         DO 180 J = KM1, NENORN                                          !QZI02430
            T = A(K,J) + U2 * A(K1,J)                                   !QZI02440
            A(K,J) = A(K,J) + T * V1                                    !QZI02450
            A(K1,J) = A(K1,J) + T * V2                                  !QZI02460
            T = B(K,J) + U2 * B(K1,J)                                   !QZI02470
            B(K,J) = B(K,J) + T * V1                                    !QZI02480
            B(K1,J) = B(K1,J) + T * V2                                  !QZI02490
  180    CONTINUE                                                       !QZI02500

         IF (K .NE. L) A(K1,KM1) = 0.0D0                                !QZI02520
         GO TO 240                                                      !QZI02530
!C     :::::::::: ZERO A(K+1,K-1) AND A(K+2,K-1) ::::::::::              !QZI02540
  190    IF (K .EQ. L) GO TO 200                                        !QZI02550
         A1 = A(K,KM1)                                                  !QZI02560
         A2 = A(K1,KM1)                                                 !QZI02570
         A3 = A(K2,KM1)                                                 !QZI02580
  200    S = ABS(A1) + ABS(A2) + ABS(A3)                             !QZI02590
         IF (S .EQ. 0.0D0) GO TO 220                                    !QZI02600
         U1 = A1 / S                                                    !QZI02610
         U2 = A2 / S                                                    !QZI02620
         U3 = A3 / S                                                    !QZI02630
         R = SIGN(SQRT(U1*U1+U2*U2+U3*U3),U1)                         !QZI02640
         V1 = -(U1 + R) / R                                             !QZI02650
         V2 = -U2 / R                                                   !QZI02660
         V3 = -U3 / R                                                   !QZI02670
         U2 = V2 / V1                                                   !QZI02680
         U3 = V3 / V1                                                   !QZI02690

         DO 210 J = KM1, NENORN                                          !QZI02710
            T = A(K,J) + U2 * A(K1,J) + U3 * A(K2,J)                    !QZI02720
            A(K,J) = A(K,J) + T * V1                                    !QZI02730
            A(K1,J) = A(K1,J) + T * V2                                  !QZI02740
            A(K2,J) = A(K2,J) + T * V3                                  !QZI02750
            T = B(K,J) + U2 * B(K1,J) + U3 * B(K2,J)                    !QZI02760
            B(K,J) = B(K,J) + T * V1                                    !QZI02770
            B(K1,J) = B(K1,J) + T * V2                                  !QZI02780
            B(K2,J) = B(K2,J) + T * V3                                  !QZI02790
  210    CONTINUE                                                       !QZI02800

         IF (K .EQ. L) GO TO 220                                        !QZI02820
         A(K1,KM1) = 0.0D0                                              !QZI02830
         A(K2,KM1) = 0.0D0                                              !QZI02840
!C     :::::::::: ZERO B(K+2,K+1) AND B(K+2,K) ::::::::::                !QZI02850
  220    S = ABS(B(K2,K2)) + ABS(B(K2,K1)) + ABS(B(K2,K))            !QZI02860
         IF (S .EQ. 0.0D0) GO TO 240                                    !QZI02870
         U1 = B(K2,K2) / S                                              !QZI02880
         U2 = B(K2,K1) / S                                              !QZI02890
         U3 = B(K2,K) / S                                               !QZI02900
         R = SIGN(SQRT(U1*U1+U2*U2+U3*U3),U1)                         !QZI02910
         V1 = -(U1 + R) / R                                             !QZI02920
         V2 = -U2 / R                                                   !QZI02930
         V3 = -U3 / R                                                   !QZI02940
         U2 = V2 / V1                                                   !QZI02950
         U3 = V3 / V1                                                   !QZI02960

         DO 230 I = LOR1, LL                                            !QZI02980
            T = A(I,K2) + U2 * A(I,K1) + U3 * A(I,K)                    !QZI02990
            A(I,K2) = A(I,K2) + T * V1                                  !QZI03000
            A(I,K1) = A(I,K1) + T * V2                                  !QZI03010
            A(I,K) = A(I,K) + T * V3                                    !QZI03020
            T = B(I,K2) + U2 * B(I,K1) + U3 * B(I,K)                    !QZI03030
            B(I,K2) = B(I,K2) + T * V1                                  !QZI03040
            B(I,K1) = B(I,K1) + T * V2                                  !QZI03050
            B(I,K) = B(I,K) + T * V3                                    !QZI03060
  230    CONTINUE                                                       !QZI03070

         B(K2,K) = 0.0D0                                                !QZI03090
         B(K2,K1) = 0.0D0                                               !QZI03100
         IF (.NOT. MATZ) GO TO 240                                      !QZI03110

         DO 235 I = 1, N                                                !QZI03130
            T = Z(I,K2) + U2 * Z(I,K1) + U3 * Z(I,K)                    !QZI03140
            Z(I,K2) = Z(I,K2) + T * V1                                  !QZI03150
            Z(I,K1) = Z(I,K1) + T * V2                                  !QZI03160
            Z(I,K) = Z(I,K) + T * V3                                    !QZI03170
  235    CONTINUE                                                       !QZI03180
!C     :::::::::: ZERO B(K+1,K) ::::::::::                               !QZI03190
  240    S = ABS(B(K1,K1)) + ABS(B(K1,K))                             !QZI03200
         IF (S .EQ. 0.0D0) GO TO 260                                    !QZI03210
         U1 = B(K1,K1) / S                                              !QZI03220
         U2 = B(K1,K) / S                                               !QZI03230
         R = SIGN(SQRT(U1*U1+U2*U2),U1)                               !QZI03240
         V1 = -(U1 + R) / R                                             !QZI03250
         V2 = -U2 / R                                                   !QZI03260
         U2 = V2 / V1                                                   !QZI03270

         DO 250 I = LOR1, LL                                            !QZI03290
            T = A(I,K1) + U2 * A(I,K)                                   !QZI03300
            A(I,K1) = A(I,K1) + T * V1                                  !QZI03310
            A(I,K) = A(I,K) + T * V2                                    !QZI03320
            T = B(I,K1) + U2 * B(I,K)                                   !QZI03330
            B(I,K1) = B(I,K1) + T * V1                                  !QZI03340
            B(I,K) = B(I,K) + T * V2                                    !QZI03350
  250    CONTINUE                                                       !QZI03360

         B(K1,K) = 0.0D0                                                !QZI03380
         IF (.NOT. MATZ) GO TO 260                                      !QZI03390

         DO 255 I = 1, N                                                !QZI03410
            T = Z(I,K1) + U2 * Z(I,K)                                   !QZI03420
            Z(I,K1) = Z(I,K1) + T * V1                                  !QZI03430
            Z(I,K) = Z(I,K) + T * V2                                    !QZI03440
  255    CONTINUE                                                       !QZI03450

  260 CONTINUE                                                          !QZI03470
!C     :::::::::: END QZ STEP ::::::::::                                 !QZI03480
      GO TO 70                                                          !QZI03490
!C     :::::::::: SET ERROR -- NEITHER BOTTOM SUBDIAGONAL ELEMENT        !QZI03500
!C                HAS BECOME NEGLIGIBLE AFTER 50 ITERATIONS ::::::::::   !QZI03510
 1000 IERR = NEN                                                         !QZI03520
!C     :::::::::: SAVE EPSB FOR USE BY QZVAL AND QZVEC ::::::::::        !QZI03530
 1001 IF (N .GT. 1) B(N,1) = EPSB                                       !QZI03540
      RETURN                                                            !QZI03550
!C     :::::::::: LAST CARD OF QZIT ::::::::::                           !QZI03560
      END                                                               !QZI03570
      

!	   QZVAL(NM,N,A,B,ALFR,ALFI,BETA,MATZ,Z,NODO2)            !QZI00030

subroutine QZVAL(N,A,B,ALFR,ALFI,BETA,MATZ,Z)
  use      def_solver

  implicit none

  integer(ip)             :: n
  real(rp)     :: a(n,n),b(n,n),z(n,n),alfr(n),alfi(n),beta(n)
  integer(ip)  :: I,J,NEN,NA,NN,ISW

  real(rp)     :: C,D,E,R,S,T,AN,A1,A2,BN,CQ,CZ,DI,DR,EI,TI,TR,U1,U2,&       
                  V1,V2,A1I,A11,A12,A2I,A21,A22,B11,B12,B22,SQI,SQR, &
                  SSI,SSR,SZI,SZR,A11I,A11R,A12I,A12R,A22I,A22R,EPSB

  logical      :: matz




      EPSB = B(N,1)                                                     !QZI00750
      ISW = 1                                                           !QZI00760
!     :::::::::: FIND EIGENVALUES OF QUASI-TRIANGULAR MATRICES.         !QZI00770
!C                FOR NEN=N STEP -1 UNTIL 1 DO -- ::::::::::              !QZI00780
      DO 510 NN = 1, N                                                  !QZI00790
         NEN = N + 1 - NN                                                !QZI00800
         NA = NEN - 1                                                    !QZI00810
         IF (ISW .EQ. 2) GO TO 505                                      !QZI00820
         IF (NEN .EQ. 1) GO TO 410                                       !QZI00830
         IF (A(NEN,NA) .NE. 0.0D0) GO TO 420                             !QZI00840
!C     :::::::::: 1-BY-1 BLOCK, ONE REAL ROOT ::::::::::                 !QZI00850
  410    ALFR(NEN) = A(NEN,NEN)                                            !QZI00860
         IF (B(NEN,NEN) .LT. 0.0D0) ALFR(NEN) = -ALFR(NEN)                  !QZI00870
         BETA(NEN) = ABS(B(NEN,NEN))                                      !QZI00880
         ALFI(NEN) = 0.0D0                                               !QZI00890
         GO TO 510                                                      !QZI00900
!C     :::::::::: 2-BY-2 BLOCK ::::::::::                                !QZI00910
  420    IF (ABS(B(NA,NA)) .LE. EPSB) GO TO 455                        !QZI00920
         IF (ABS(B(NEN,NEN)) .GT. EPSB) GO TO 430                        !QZI00930
         A1 = A(NEN,NEN)                                                  !QZI00940
         A2 = A(NEN,NA)                                                  !QZI00950
         BN = 0.0D0                                                     !QZI00960
         GO TO 435                                                      !QZI00970
  430    AN = ABS(A(NA,NA)) + ABS(A(NA,NEN)) + ABS(A(NEN,NA)) &         !QZI00980
           + ABS(A(NEN,NEN))                                            !QZI00990
         BN = ABS(B(NA,NA)) + ABS(B(NA,NEN)) + ABS(B(NEN,NEN))          !QZI01000
         A11 = A(NA,NA) / AN                                            !QZI01010
         A12 = A(NA,NEN) / AN                                            !QZI01020
         A21 = A(NEN,NA) / AN                                            !QZI01030
         A22 = A(NEN,NEN) / AN                                            !QZI01040
         B11 = B(NA,NA) / BN                                            !QZI01050
         B12 = B(NA,NEN) / BN                                            !QZI01060
         B22 = B(NEN,NEN) / BN                                            !QZI01070
         E = A11 / B11                                                  !QZI01080
         C = 0.5D0 * ((A22 - E * B22) / B22 - A21 * B12 / (B11 * B22))  !QZI01090
         D = C * C + A21 * (A12 - E * B12) / (B11 * B22)                !QZI01100
         IF (D .LT. 0.0D0) GO TO 480                                    !QZI01110
!C     :::::::::: TWO REAL ROOTS.                                        !QZI01120
!C                ZERO BOTH A(NEN,NA) AND B(NEN,NA) ::::::::::             !QZI01130
         E = E + C + SIGN(SQRT(D),C)                                  !QZI01140
         A11 = A11 - E * B11                                            !QZI01150
         A12 = A12 - E * B12                                            !QZI01160
         A22 = A22 - E * B22                                            !QZI01170
         IF (ABS(A11) + ABS(A12) .LT.    &                             !QZI01180
                       ABS(A21) + ABS(A22)) GO TO 432                           !QZI01190
         A1 = A12                                                       !QZI01200
         A2 = A11                                                       !QZI01210
         GO TO 435                                                      !QZI01220
  432    A1 = A22                                                       !QZI01230
         A2 = A21                                                       !QZI01240
!C     :::::::::: CHOOSE AND APPLY REAL Z ::::::::::                     !QZI01250
  435    S = ABS(A1) + ABS(A2)                                        !QZI01260
         U1 = A1 / S                                                    !QZI01270
         U2 = A2 / S                                                    !QZI01280
         R = SIGN(SQRT(U1*U1+U2*U2),U1)                               !QZI01290
         V1 = -(U1 + R) / R                                             !QZI01300
         V2 = -U2 / R                                                   !QZI01310
         U2 = V2 / V1                                                   !QZI01320

         DO 440 I = 1, NEN                                               !QZI01340
            T = A(I,NEN) + U2 * A(I,NA)                                  !QZI01350
            A(I,NEN) = A(I,NEN) + T * V1                                  !QZI01360
            A(I,NA) = A(I,NA) + T * V2                                  !QZI01370
            T = B(I,NEN) + U2 * B(I,NA)                                  !QZI01380
            B(I,NEN) = B(I,NEN) + T * V1                                  !QZI01390
            B(I,NA) = B(I,NA) + T * V2                                  !QZI01400
  440    CONTINUE                                                       !QZI01410
         
		 IF (.NOT. MATZ) GO TO 450                                      !QZI01430

         DO 445 I = 1, N                                                !QZI01450
            T = Z(I,NEN) + U2 * Z(I,NA)                                  !QZI01460
            Z(I,NEN) = Z(I,NEN) + T * V1                                  !QZI01470
            Z(I,NA) = Z(I,NA) + T * V2                                  !QZI01480
  445    CONTINUE                                                       !QZI01490

  450    IF (BN .EQ. 0.0D0) GO TO 475                                   !QZI01510
         IF (AN .LT. ABS(E) * BN) GO TO 455                            !QZI01520
         A1 = B(NA,NA)                                                  !QZI01530
         A2 = B(NEN,NA)                                                  !QZI01540
         GO TO 460                                                      !QZI01550
  455    A1 = A(NA,NA)                                                  !QZI01560
         A2 = A(NEN,NA)                                                  !QZI01570
!C     :::::::::: CHOOSE AND APPLY REAL Q ::::::::::                     !QZI01580
  460    S = ABS(A1) + ABS(A2)                                        !QZI01590
         IF (S .EQ. 0.0D0) GO TO 475                                    !QZI01600
         U1 = A1 / S                                                    !QZI01610
         U2 = A2 / S                                                    !QZI01620
         R = SIGN(SQRT(U1*U1+U2*U2),U1)                               !QZI01630
         V1 = -(U1 + R) / R                                             !QZI01640
         V2 = -U2 / R                                                   !QZI01650
         U2 = V2 / V1                                                   !QZI01660

         DO 470 J = NA, N                                               !QZI01680
            T = A(NA,J) + U2 * A(NEN,J)                                  !QZI01690
            A(NA,J) = A(NA,J) + T * V1                                  !QZI01700
            A(NEN,J) = A(NEN,J) + T * V2                                  !QZI01710
            T = B(NA,J) + U2 * B(NEN,J)                                  !QZI01720
            B(NA,J) = B(NA,J) + T * V1                                  !QZI01730
            B(NEN,J) = B(NEN,J) + T * V2                                  !QZI01740
  470    CONTINUE                                                       !QZI01750

  475    A(NEN,NA) = 0.0D0                                               !QZI01770
         B(NEN,NA) = 0.0D0                                               !QZI01780
         ALFR(NA) = A(NA,NA)                                            !QZI01790
         ALFR(NEN) = A(NEN,NEN)                                            !QZI01800
         IF (B(NA,NA) .LT. 0.0D0) ALFR(NA) = -ALFR(NA)                  !QZI01810
         IF (B(NEN,NEN) .LT. 0.0D0) ALFR(NEN) = -ALFR(NEN)                  !QZI01820
         BETA(NA) = ABS(B(NA,NA))                                      !QZI01830
         BETA(NEN) = ABS(B(NEN,NEN))                                      !QZI01840
         ALFI(NEN) = 0.0D0                                               !QZI01850
         ALFI(NA) = 0.0D0                                               !QZI01860
         GO TO 505                                                      !QZI01870
!C     :::::::::: TWO COMPLEX ROOTS ::::::::::                           !QZI01880
  480    E = E + C                                                      !QZI01890
         EI = SQRT(-D)                                                 !QZI01900
         A11R = A11 - E * B11                                           !QZI01910
         A11I = EI * B11                                                !QZI01920
         A12R = A12 - E * B12                                           !QZI01930
         A12I = EI * B12                                                !QZI01940
         A22R = A22 - E * B22                                           !QZI01950
         A22I = EI * B22                                                !QZI01960
         IF (ABS(A11R) + ABS(A11I) + ABS(A12R) + ABS(A12I) .LT.  &  !QZI01970
            ABS(A21) + ABS(A22R) + ABS(A22I)) GO TO 482             !QZI01980
         A1 = A12R                                                      !QZI01990
         A1I = A12I                                                     !QZI02000
         A2 = -A11R                                                     !QZI02010
         A2I = -A11I                                                    !QZI02020
         GO TO 485                                                      !QZI02030
  482    A1 = A22R                                                      !QZI02040
         A1I = A22I                                                     !QZI02050
         A2 = -A21                                                      !QZI02060
         A2I = 0.0D0                                                    !QZI02070
!C     :::::::::: CHOOSE COMPLEX Z ::::::::::                            !QZI02080
  485    CZ = SQRT(A1*A1+A1I*A1I)                                      !QZI02090
         IF (CZ .EQ. 0.0D0) GO TO 487                                   !QZI02100
         SZR = (A1 * A2 + A1I * A2I) / CZ                               !QZI02110
         SZI = (A1 * A2I - A1I * A2) / CZ                               !QZI02120
         R = SQRT(CZ*CZ+SZR*SZR+SZI*SZI)                               !QZI02130
         CZ = CZ / R                                                    !QZI02140
         SZR = SZR / R                                                  !QZI02150
         SZI = SZI / R                                                  !QZI02160
         GO TO 490                                                      !QZI02170
  487    SZR = 1.0D0                                                    !QZI02180
         SZI = 0.0D0                                                    !QZI02190
  490    IF (AN .LT. (ABS(E) + EI) * BN) GO TO 492                     !QZI02200
         A1 = CZ * B11 + SZR * B12                                      !QZI02210
         A1I = SZI * B12                                                !QZI02220
         A2 = SZR * B22                                                 !QZI02230
         A2I = SZI * B22                                                !QZI02240
         GO TO 495                                                      !QZI02250
  492    A1 = CZ * A11 + SZR * A12                                      !QZI02260
         A1I = SZI * A12                                                !QZI02270
         A2 = CZ * A21 + SZR * A22                                      !QZI02280
         A2I = SZI * A22                                                !QZI02290
!C     :::::::::: CHOOSE COMPLEX Q ::::::::::                            !QZI02300
  495    CQ = SQRT(A1*A1+A1I*A1I)                                      !QZI02310
         IF (CQ .EQ. 0.0D0) GO TO 497                                   !QZI02320
         SQR = (A1 * A2 + A1I * A2I) / CQ                               !QZI02330
         SQI = (A1 * A2I - A1I * A2) / CQ                               !QZI02340
         R = SQRT(CQ*CQ+SQR*SQR+SQI*SQI)                               !QZI02350
         CQ = CQ / R                                                    !QZI02360
         SQR = SQR / R                                                  !QZI02370
         SQI = SQI / R                                                  !QZI02380
         GO TO 500                                                      !QZI02390
  497    SQR = 1.0D0                                                    !QZI02400
         SQI = 0.0D0                                                    !QZI02410
!C     :::::::::: COMPUTE DIAGONAL ELEMENTS THAT WOULD RESULT            !QZI02420
!C                IF TRANSFORMATIONS WERE APPLIED ::::::::::             !QZI02430
  500    SSR = SQR * SZR + SQI * SZI                                    !QZI02440
         SSI = SQR * SZI - SQI * SZR                                    !QZI02450
         I = 1                                                          !QZI02460
         TR = CQ * CZ * A11 + CQ * SZR * A12 + SQR * CZ * A21    &       !QZI02470
           + SSR * A22                                                 !QZI02480
         TI = CQ * SZI * A12 - SQI * CZ * A21 + SSI * A22               !QZI02490
         DR = CQ * CZ * B11 + CQ * SZR * B12 + SSR * B22                !QZI02500
         DI = CQ * SZI * B12 + SSI * B22                                !QZI02510
         GO TO 503                                                      !QZI02520
  502    I = 2                                                          !QZI02530
         TR = SSR * A11 - SQR * CZ * A12 - CQ * SZR * A21        &       !QZI02540
           + CQ * CZ * A22                                             !QZI02550
         TI = -SSI * A11 - SQI * CZ * A12 + CQ * SZI * A21              !QZI02560
         DR = SSR * B11 - SQR * CZ * B12 + CQ * CZ * B22                !QZI02570
         DI = -SSI * B11 - SQI * CZ * B12                               !QZI02580
  503    T = TI * DR - TR * DI                                          !QZI02590
         J = NA                                                         !QZI02600
         IF (T .LT. 0.0D0) J = NEN                                       !QZI02610
         R = SQRT(DR*DR+DI*DI)                                         !QZI02620
         BETA(J) = BN * R                                               !QZI02630
         ALFR(J) = AN * (TR * DR + TI * DI) / R                         !QZI02640
         ALFI(J) = AN * T / R                                           !QZI02650
         IF (I .EQ. 1) GO TO 502                                        !QZI02660
  505    ISW = 3 - ISW                                                  !QZI02670
  510 CONTINUE                                                          !QZI02680

      RETURN                                                            !QZI02700
!C     :::::::::: LAST CARD OF QZVAL ::::::::::                          !QZI02710
      END                                                               !QZI02720
      
subroutine QZVEC(N,A,B,ALFR,ALFI,BETA,Z)                 !QZI00030
  use      def_solver

  implicit none

  integer(ip)             :: n
  real(rp)     :: a(n,n),b(n,n),z(n,n), alfr(n),alfi(n),beta(n)
  integer(ip)  :: I,J,K,M,NEN,II,JJ,NA,NM,NN,ISW,NENM2

  real(rp)     :: D,Q,R,S,T,W,X,Y,DI,DR,RA,RR,SA,TI,TR,T1,T2,W1,X1,ZZ,Z1, &  
                  ALFM,ALMI,ALMR,BETM,EPSB 


  

      EPSB = B(N,1)                                                     !QZI00720
      ISW = 1                                                           !QZI00730
!C     :::::::::: FOR NEN=N STEP -1 UNTIL 1 DO -- ::::::::::              !QZI00740
      DO 800 NN = 1, N                                                  !QZI00750
         NEN = N + 1 - NN                                                !QZI00760
         NA = NEN - 1                                                    !QZI00770
         IF (ISW .EQ. 2) GO TO 795                                      !QZI00780
         IF (ALFI(NEN) .NE. 0.0D0) GO TO 710                             !QZI00790
!C     :::::::::: REAL VECTOR ::::::::::                                 !QZI00800
         M = NEN                                                         !QZI00810
         B(NEN,NEN) = 1.0D0                                               !QZI00820
         IF (NA .EQ. 0) GO TO 800                                       !QZI00830
         ALFM = ALFR(M)                                                 !QZI00840
         BETM = BETA(M)                                                 !QZI00850
!C     :::::::::: FOR I=NEN-1 STEP -1 UNTIL 1 DO -- ::::::::::            !QZI00860
         DO 700 II = 1, NA                                              !QZI00870
            I = NEN - II                                                 !QZI00880
            W = BETM * A(I,I) - ALFM * B(I,I)                           !QZI00890
            R = 0.0D0                                                   !QZI00900

            DO  J = M, NEN
				R = R + (BETM * A(I,J) - ALFM * B(I,J)) * B(J,NEN)
            enddo

            IF (I .EQ. 1 .OR. ISW .EQ. 2) GO TO 630                     !QZI00950
            IF (BETM * A(I,I-1) .EQ. 0.0D0) GO TO 630                   !QZI00960
            ZZ = W                                                      !QZI00970
            S = R                                                       !QZI00980
            GO TO 690                                                   !QZI00990
  630       M = I                                                       !QZI01000
            IF (ISW .EQ. 2) GO TO 640                                   !QZI01010
!C     :::::::::: REAL 1-BY-1 BLOCK ::::::::::                           !QZI01020
            T = W                                                       !QZI01030
            IF (W .EQ. 0.0D0) T = EPSB                                  !QZI01040
            B(I,NEN) = -R / T                                            !QZI01050
            GO TO 700                                                   !QZI01060
!C     :::::::::: REAL 2-BY-2 BLOCK ::::::::::                           !QZI01070
  640       X = BETM * A(I,I+1) - ALFM * B(I,I+1)                       !QZI01080
            Y = BETM * A(I+1,I)                                         !QZI01090
            Q = W * ZZ - X * Y                                          !QZI01100
            T = (X * S - ZZ * R) / Q                                    !QZI01110
            B(I,NEN) = T                                                 !QZI01120
            IF (ABS(X) .LE. ABS(ZZ)) GO TO 650                        !QZI01130
            B(I+1,NEN) = (-R - W * T) / X                                !QZI01140
            GO TO 690                                                   !QZI01150
  650       B(I+1,NEN) = (-S - Y * T) / ZZ                               !QZI01160
  690       ISW = 3 - ISW                                               !QZI01170
  700    CONTINUE                                                       !QZI01180
!C     :::::::::: END REAL VECTOR ::::::::::                             !QZI01190
         GO TO 800                                                      !QZI01200
!C     :::::::::: COMPLEX VECTOR ::::::::::                              !QZI01210
  710    M = NA                                                         !QZI01220
         ALMR = ALFR(M)                                                 !QZI01230
         ALMI = ALFI(M)                                                 !QZI01240
         BETM = BETA(M)                                                 !QZI01250
!C     :::::::::: LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT         !QZI01260
!C                EIGENVECTOR MATRIX IS TRIANGULAR ::::::::::            !QZI01270
         Y = BETM * A(NEN,NA)                                            !QZI01280
         B(NA,NA) = -ALMI * B(NEN,NEN) / Y                                !QZI01290
         B(NA,NEN) = (ALMR * B(NEN,NEN) - BETM * A(NEN,NEN)) / Y             !QZI01300
         B(NEN,NA) = 0.0D0                                               !QZI01310
         B(NEN,NEN) = 1.0D0                                               !QZI01320
         NENM2 = NA - 1                                                  !QZI01330
         IF (NENM2 .EQ. 0) GO TO 795                                     !QZI01340
!C     :::::::::: FOR I=NEN-2 STEP -1 UNTIL 1 DO -- ::::::::::            !QZI01350
         DO 790 II = 1, NENM2                                            !QZI01360
            I = NA - II                                                 !QZI01370
            W = BETM * A(I,I) - ALMR * B(I,I)                           !QZI01380
            W1 = -ALMI * B(I,I)                                         !QZI01390
            RA = 0.0D0                                                  !QZI01400
            SA = 0.0D0                                                  !QZI01410

            DO 760 J = M, NEN                                            !QZI01430
               X = BETM * A(I,J) - ALMR * B(I,J)                        !QZI01440
               X1 = -ALMI * B(I,J)                                      !QZI01450
               RA = RA + X * B(J,NA) - X1 * B(J,NEN)                     !QZI01460
               SA = SA + X * B(J,NEN) + X1 * B(J,NA)                     !QZI01470
  760       CONTINUE                                                    !QZI01480
!C                                                                       !QZI01490
            IF (I .EQ. 1 .OR. ISW .EQ. 2) GO TO 770                     !QZI01500
            IF (BETM * A(I,I-1) .EQ. 0.0D0) GO TO 770                   !QZI01510
            ZZ = W                                                      !QZI01520
            Z1 = W1                                                     !QZI01530
            R = RA                                                      !QZI01540
            S = SA                                                      !QZI01550
            ISW = 2                                                     !QZI01560
            GO TO 790                                                   !QZI01570
  770       M = I                                                       !QZI01580
            IF (ISW .EQ. 2) GO TO 780                                   !QZI01590
!C     :::::::::: COMPLEX 1-BY-1 BLOCK ::::::::::                        !QZI01600
            TR = -RA                                                    !QZI01610
            TI = -SA                                                    !QZI01620
  773       DR = W                                                      !QZI01630
            DI = W1                                                     !QZI01640
!C     :::::::::: COMPLEX DIVIDE (T1,T2) = (TR,TI) / (DR,DI) ::::::::::  !QZI01650
  775       IF (ABS(DI) .GT. ABS(DR)) GO TO 777                       !QZI01660
            RR = DI / DR                                                !QZI01670
            D = DR + DI * RR                                            !QZI01680
            T1 = (TR + TI * RR) / D                                     !QZI01690
            T2 = (TI - TR * RR) / D                                     !QZI01700
            GO TO (787,782), ISW                                        !QZI01710
  777       RR = DR / DI                                                !QZI01720
            D = DR * RR + DI                                            !QZI01730
            T1 = (TR * RR + TI) / D                                     !QZI01740
            T2 = (TI * RR - TR) / D                                     !QZI01750
            GO TO (787,782), ISW                                        !QZI01760
!C     :::::::::: COMPLEX 2-BY-2 BLOCK ::::::::::                        !QZI01770
  780       X = BETM * A(I,I+1) - ALMR * B(I,I+1)                       !QZI01780
            X1 = -ALMI * B(I,I+1)                                       !QZI01790
            Y = BETM * A(I+1,I)                                         !QZI01800
            TR = Y * RA - W * R + W1 * S                                !QZI01810
            TI = Y * SA - W * S - W1 * R                                !QZI01820
            DR = W * ZZ - W1 * Z1 - X * Y                               !QZI01830
            DI = W * Z1 + W1 * ZZ - X1 * Y                              !QZI01840
            IF (DR .EQ. 0.0D0 .AND. DI .EQ. 0.0D0) DR = EPSB            !QZI01850
            GO TO 775                                                   !QZI01860
  782       B(I+1,NA) = T1                                              !QZI01870
            B(I+1,NEN) = T2                                              !QZI01880
            ISW = 1                                                     !QZI01890
            IF (ABS(Y) .GT. ABS(W) + ABS(W1)) GO TO 785              !QZI01900
            TR = -RA - X * B(I+1,NA) + X1 * B(I+1,NEN)                   !QZI01910
            TI = -SA - X * B(I+1,NEN) - X1 * B(I+1,NA)                   !QZI01920
            GO TO 773                                                   !QZI01930
  785       T1 = (-R - ZZ * B(I+1,NA) + Z1 * B(I+1,NEN)) / Y             !QZI01940
            T2 = (-S - ZZ * B(I+1,NEN) - Z1 * B(I+1,NA)) / Y             !QZI01950
  787       B(I,NA) = T1                                                !QZI01960
            B(I,NEN) = T2                                                !QZI01970
  790    CONTINUE                                                       !QZI01980
!C     :::::::::: END COMPLEX VECTOR ::::::::::                          !QZI01990
  795    ISW = 3 - ISW                                                  !QZI02000
  800 CONTINUE                                                          !QZI02010
!C     :::::::::: END BACK SUBSTITUTION.                                 !QZI02020
!C                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.               !QZI02030
!C                FOR J=N STEP -1 UNTIL 1 DO -- ::::::::::               !QZI02040
      DO JJ = 1, N                                                  !QZI02050
         J = N + 1 - JJ                                                 !QZI02060

         DO I = 1, N                                                !QZI02080
            ZZ = 0.0D0                                                  !QZI02090

            DO  K = 1, J  
               ZZ = ZZ + Z(I,K) * B(K,J) 
            enddo
			                        
            Z(I,J) = ZZ                                                 !QZI02140
         enddo
      enddo 
  880 CONTINUE                                                          !QZI02150
!C     :::::::::: NORMALIZE SO THAT MODULUS OF LARGEST                   !QZI02160
!C                COMPONENT OF EACH VECTOR IS 1.                         !QZI02170
!C                (ISW IS 1 INITIALLY FROM BEFORE) ::::::::::            !QZI02180
      DO 950 J = 1,N                                                    !QZI02190
         D = 0.0D0                                                      !QZI02200
         IF (ISW .EQ. 2) GO TO 920                                      !QZI02210
         IF (ALFI(J) .NE. 0.0D0) GO TO 945                              !QZI02220

         DO 890 I = 1, N                                                !QZI02240
            IF (ABS(Z(I,J)) .GT. D) D = ABS(Z(I,J))                   !QZI02250
  890    CONTINUE                                                       !QZI02260

         DO I = 1, N
		    Z(I,J) = Z(I,J) / D                                            
	     enddo

         GO TO 950                                                      !QZI02310

  920    DO 930 I = 1, N                                                !QZI02330
            R = ABS(Z(I,J-1)) + ABS(Z(I,J))                           !QZI02340
            IF (R .NE. 0.0D0) R = R * SQRT((Z(I,J-1)/R)**2   &          !QZI02350
                                          +(Z(I,J)/R)**2)              !QZI02360
            IF (R .GT. D) D = R                                         !QZI02370
  930    CONTINUE                                                       !QZI02380

         DO 940 I = 1, N                                                !QZI02400
            Z(I,J-1) = Z(I,J-1) / D                                     !QZI02410
            Z(I,J) = Z(I,J) / D                                         !QZI02420
  940    CONTINUE                                                       !QZI02430

  945    ISW = 3 - ISW                                                  !QZI02450
  950 CONTINUE                                                          !QZI02460

      RETURN                                                            !QZI02480
!C     :::::::::: LAST CARD OF QZVEC ::::::::::                          !QZI02490
      END                                                               !QZI02500
      
      
