subroutine ibm_model_wind_turbine(itask,iimbo)
  !-----------------------------------------------------------------------
  !****f* Immbou/ibm_model_wind_turbine
  ! NAME
  !    ibm_model_wind_turbine
  ! DESCRIPTION
  !    Read IB
  ! OUTPUT
  ! USED BY
  !    ibm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use def_immbou
  use def_inpout
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), intent(in) :: iimbo
  integer(ip)             :: iboib,ipara,idime,ipoib,ii
  real(rp)                :: r0,r1,A0,A1,ww,rotat(3,3)
  real(rp)                :: nx0,ny0,nz0,nx1,ny1,nz1
  real(rp)                :: d1,coso,sino,theta,xx(3)

  select case ( itask )

  case ( 1_ip )
     !
     ! Dimensions
     !
     imbou(iimbo) % npari     =   1
     imbou(iimbo) % nparr     =  14
     imbou(iimbo) % npoib     =  68
     imbou(iimbo) % nboib     = 132
     imbou(iimbo) % kfl_typeb =   0
     imbou(iimbo) % kfl_coupl =   1
     lexib(TRI03)             =   1

  case ( 2_ip )
     !
     ! Geometry
     !
     do iboib = 1,imbou(iimbo) % nboib
        imbou(iimbo) % ltyib(iboib) = TRI03
     end do

     imbou(iimbo) % cooib(1:3, 1) = (/     0.293892626   , -0.404508497    ,  1.0 /)
     imbou(iimbo) % cooib(1:3, 2) = (/     0.154508497   , -0.475528258    ,  1.0 /)
     imbou(iimbo) % cooib(1:3, 3) = (/             0.0   ,         -0.5    ,  1.0 /)
     imbou(iimbo) % cooib(1:3, 4) = (/    -0.246478345   , -0.246478345    ,  1.0 /)
     imbou(iimbo) % cooib(1:3, 5) = (/     0.158248839   , -0.310580834    ,  1.0 /)
     imbou(iimbo) % cooib(1:3, 6) = (/   -0.0545288337   , -0.344281506    ,  1.0 /)
     imbou(iimbo) % cooib(1:3, 7) = (/     0.404508497   , -0.293892626    ,  0.0 /)
     imbou(iimbo) % cooib(1:3, 8) = (/     0.293892626   , -0.404508497    ,  0.0 /)
     imbou(iimbo) % cooib(1:3, 9) = (/     0.154508497   , -0.475528258    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,10) = (/     0.475528258   , -0.154508497    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,11) = (/     0.404508497   , -0.293892626    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,12) = (/    -0.154508497   , -0.475528258    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,13) = (/    -0.154508497   , -0.475528258    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,14) = (/    -0.293892626   , -0.404508497    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,15) = (/    -0.404508497   , -0.293892626    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,16) = (/    -0.475528258   , -0.154508497    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,17) = (/             0.0   ,         -0.5    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,18) = (/    0.0828601678   , -0.162622236    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,19) = (/    -0.146549106   , -0.074670499    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,20) = (/     0.344281506   , 0.0545288337    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,21) = (/    -0.344281506   ,-0.0545288337    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,22) = (/      0.29660958   , -0.116397518    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,23) = (/    0.0545288337   , -0.344281506    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,24) = (/     0.246478345   , -0.246478345    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,25) = (/    -0.131218882   , -0.340864889    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,26) = (/     0.475528258   , -0.154508497    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,27) = (/     0.404508497   ,  0.293892626    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,28) = (/     0.475528258   ,  0.154508497    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,29) = (/             0.5   ,-3.061617e-17    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,30) = (/    -0.293892626   , -0.404508497    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,31) = (/    -0.404508497   , -0.293892626    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,32) = (/    -0.475528258   , -0.154508497    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,33) = (/            -0.5   , 3.061617e-17    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,34) = (/            -0.5   , 3.061617e-17    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,35) = (/    -0.475528258   ,  0.154508497    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,36) = (/    -0.116302029   ,  0.116302029    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,37) = (/     0.110839793   , 0.0522927344    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,38) = (/    -0.333144855   ,  0.136792525    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,39) = (/     0.246478345   ,  0.246478345    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,40) = (/   -0.0828601678   , -0.162622236    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,41) = (/     0.146549106   , -0.074670499    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,42) = (/     0.344281506   ,-0.0545288337    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,43) = (/    -0.310580834   , -0.158248839    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,44) = (/    -0.344281506   , 0.0545288337    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,45) = (/     0.404508497   ,  0.293892626    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,46) = (/     0.475528258   ,  0.154508497    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,47) = (/             0.5   ,-3.061617e-17    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,48) = (/     0.293892626   ,  0.404508497    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,49) = (/    -0.475528258   ,  0.154508497    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,50) = (/    -0.404508497   ,  0.293892626    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,51) = (/    -0.404508497   ,  0.293892626    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,52) = (/    0.0545288337   ,  0.344281506    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,53) = (/    -0.158248839   ,  0.310580834    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,54) = (/    0.0828601678   ,  0.162622236    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,55) = (/    -0.146549106   ,  0.074670499    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,56) = (/    -0.246478345   ,  0.246478345    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,57) = (/     0.310580834   ,  0.158248839    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,58) = (/     0.293892626   ,  0.404508497    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,59) = (/     0.154508497   ,  0.475528258    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,60) = (/    -0.293892626   ,  0.404508497    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,61) = (/    -0.293892626   ,  0.404508497    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,62) = (/     0.131218882   ,  0.340864889    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,63) = (/   -0.0545288337   ,  0.344281506    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,64) = (/     0.154508497   ,  0.475528258    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,65) = (/    -0.154508497   ,  0.475528258    ,  0.0 /)
     imbou(iimbo) % cooib(1:3,66) = (/    -0.154508497   ,  0.475528258    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,67) = (/             0.0   ,          0.5    ,  1.0 /)
     imbou(iimbo) % cooib(1:3,68) = (/             0.0   ,          0.5    ,  0.0 /)

     imbou(iimbo) % lnoib(1:3,    1) = (/      67   ,    66    ,   68 /)
     imbou(iimbo) % lnoib(1:3,    2) = (/      66   ,    65    ,   68 /)
     imbou(iimbo) % lnoib(1:3,    3) = (/      66   ,    61    ,   65 /)
     imbou(iimbo) % lnoib(1:3,    4) = (/      61   ,    60    ,   65 /)
     imbou(iimbo) % lnoib(1:3,    5) = (/      61   ,    51    ,   60 /)
     imbou(iimbo) % lnoib(1:3,    6) = (/      51   ,    50    ,   60 /)
     imbou(iimbo) % lnoib(1:3,    7) = (/      51   ,    35    ,   50 /)
     imbou(iimbo) % lnoib(1:3,    8) = (/      35   ,    49    ,   50 /)
     imbou(iimbo) % lnoib(1:3,    9) = (/      35   ,    34    ,   49 /)
     imbou(iimbo) % lnoib(1:3,   10) = (/      34   ,    33    ,   49 /)
     imbou(iimbo) % lnoib(1:3,   11) = (/      34   ,    16    ,   33 /)
     imbou(iimbo) % lnoib(1:3,   12) = (/      16   ,    32    ,   33 /)
     imbou(iimbo) % lnoib(1:3,   13) = (/      16   ,    15    ,   32 /)
     imbou(iimbo) % lnoib(1:3,   14) = (/      15   ,    31    ,   32 /)
     imbou(iimbo) % lnoib(1:3,   15) = (/      15   ,    14    ,   31 /)
     imbou(iimbo) % lnoib(1:3,   16) = (/      14   ,    30    ,   31 /)
     imbou(iimbo) % lnoib(1:3,   17) = (/      14   ,    13    ,   30 /)
     imbou(iimbo) % lnoib(1:3,   18) = (/      13   ,    12    ,   30 /)
     imbou(iimbo) % lnoib(1:3,   19) = (/      13   ,     3    ,   12 /)
     imbou(iimbo) % lnoib(1:3,   20) = (/       3   ,    17    ,   12 /)
     imbou(iimbo) % lnoib(1:3,   21) = (/       3   ,     2    ,   17 /)
     imbou(iimbo) % lnoib(1:3,   22) = (/       2   ,     9    ,   17 /)
     imbou(iimbo) % lnoib(1:3,   23) = (/       2   ,     1    ,    9 /)
     imbou(iimbo) % lnoib(1:3,   24) = (/       1   ,     8    ,    9 /)
     imbou(iimbo) % lnoib(1:3,   25) = (/       1   ,    11    ,    8 /)
     imbou(iimbo) % lnoib(1:3,   26) = (/      11   ,     7    ,    8 /)
     imbou(iimbo) % lnoib(1:3,   27) = (/      11   ,    10    ,    7 /)
     imbou(iimbo) % lnoib(1:3,   28) = (/      10   ,    26    ,    7 /)
     imbou(iimbo) % lnoib(1:3,   29) = (/      10   ,    29    ,   26 /)
     imbou(iimbo) % lnoib(1:3,   30) = (/      29   ,    47    ,   26 /)
     imbou(iimbo) % lnoib(1:3,   31) = (/      29   ,    28    ,   47 /)
     imbou(iimbo) % lnoib(1:3,   32) = (/      28   ,    46    ,   47 /)
     imbou(iimbo) % lnoib(1:3,   33) = (/      28   ,    27    ,   46 /)
     imbou(iimbo) % lnoib(1:3,   34) = (/      27   ,    45    ,   46 /)
     imbou(iimbo) % lnoib(1:3,   35) = (/      27   ,    48    ,   45 /)
     imbou(iimbo) % lnoib(1:3,   36) = (/      48   ,    58    ,   45 /)
     imbou(iimbo) % lnoib(1:3,   37) = (/      48   ,    59    ,   58 /)
     imbou(iimbo) % lnoib(1:3,   38) = (/      59   ,    64    ,   58 /)
     imbou(iimbo) % lnoib(1:3,   39) = (/      59   ,    67    ,   64 /)
     imbou(iimbo) % lnoib(1:3,   40) = (/      67   ,    68    ,   64 /)
     imbou(iimbo) % lnoib(1:3,   41) = (/      46   ,    45    ,   57 /)
     imbou(iimbo) % lnoib(1:3,   42) = (/      30   ,    12    ,   25 /)
     imbou(iimbo) % lnoib(1:3,   43) = (/      68   ,    65    ,   63 /)
     imbou(iimbo) % lnoib(1:3,   44) = (/       8   ,     7    ,   24 /)
     imbou(iimbo) % lnoib(1:3,   45) = (/      49   ,    33    ,   44 /)
     imbou(iimbo) % lnoib(1:3,   46) = (/      58   ,    64    ,   62 /)
     imbou(iimbo) % lnoib(1:3,   47) = (/      60   ,    50    ,   56 /)
     imbou(iimbo) % lnoib(1:3,   48) = (/      17   ,     9    ,   23 /)
     imbou(iimbo) % lnoib(1:3,   49) = (/      32   ,    31    ,   43 /)
     imbou(iimbo) % lnoib(1:3,   50) = (/      26   ,    47    ,   42 /)
     imbou(iimbo) % lnoib(1:3,   51) = (/      65   ,    60    ,   56 /)
     imbou(iimbo) % lnoib(1:3,   52) = (/      45   ,    58    ,   57 /)
     imbou(iimbo) % lnoib(1:3,   53) = (/      12   ,    17    ,   25 /)
     imbou(iimbo) % lnoib(1:3,   54) = (/       7   ,    26    ,   24 /)
     imbou(iimbo) % lnoib(1:3,   55) = (/      33   ,    32    ,   43 /)
     imbou(iimbo) % lnoib(1:3,   56) = (/       9   ,     8    ,   24 /)
     imbou(iimbo) % lnoib(1:3,   57) = (/      50   ,    49    ,   56 /)
     imbou(iimbo) % lnoib(1:3,   58) = (/      64   ,    68    ,   62 /)
     imbou(iimbo) % lnoib(1:3,   59) = (/      31   ,    30    ,   43 /)
     imbou(iimbo) % lnoib(1:3,   60) = (/      47   ,    46    ,   57 /)
     imbou(iimbo) % lnoib(1:3,   61) = (/      44   ,    33    ,   43 /)
     imbou(iimbo) % lnoib(1:3,   62) = (/      63   ,    65    ,   56 /)
     imbou(iimbo) % lnoib(1:3,   63) = (/      42   ,    47    ,   57 /)
     imbou(iimbo) % lnoib(1:3,   64) = (/      23   ,     9    ,   24 /)
     imbou(iimbo) % lnoib(1:3,   65) = (/      49   ,    44    ,   56 /)
     imbou(iimbo) % lnoib(1:3,   66) = (/      68   ,    63    ,   62 /)
     imbou(iimbo) % lnoib(1:3,   67) = (/      58   ,    62    ,   57 /)
     imbou(iimbo) % lnoib(1:3,   68) = (/      30   ,    25    ,   43 /)
     imbou(iimbo) % lnoib(1:3,   69) = (/      26    ,   42     ,  24 /)
     imbou(iimbo) % lnoib(1:3,   70) = (/      17    ,   23     ,  25 /)
     imbou(iimbo) % lnoib(1:3,   71) = (/      24    ,   42     ,  41 /)
     imbou(iimbo) % lnoib(1:3,   72) = (/      24    ,   41     ,  23 /)
     imbou(iimbo) % lnoib(1:3,   73) = (/      41    ,   42     ,  57 /)
     imbou(iimbo) % lnoib(1:3,   74) = (/      56    ,   44     ,  55 /)
     imbou(iimbo) % lnoib(1:3,   75) = (/      55    ,   44     ,  43 /)
     imbou(iimbo) % lnoib(1:3,   76) = (/      56    ,   55     ,  63 /)
     imbou(iimbo) % lnoib(1:3,   77) = (/      57    ,   62     ,  54 /)
     imbou(iimbo) % lnoib(1:3,   78) = (/      57    ,   54     ,  41 /)
     imbou(iimbo) % lnoib(1:3,   79) = (/      41    ,   54     ,  55 /)
     imbou(iimbo) % lnoib(1:3,   80) = (/      54    ,   62     ,  63 /)
     imbou(iimbo) % lnoib(1:3,   81) = (/      25    ,   23     ,  40 /)
     imbou(iimbo) % lnoib(1:3,   82) = (/      25    ,   40     ,  43 /)
     imbou(iimbo) % lnoib(1:3,   83) = (/      40    ,   23     ,  41 /)
     imbou(iimbo) % lnoib(1:3,   84) = (/      40    ,   41     ,  55 /)
     imbou(iimbo) % lnoib(1:3,   85) = (/      40    ,   55     ,  43 /)
     imbou(iimbo) % lnoib(1:3,   86) = (/      55    ,   54     ,  63 /)
     imbou(iimbo) % lnoib(1:3,   87) = (/      61    ,   66     ,  53 /)
     imbou(iimbo) % lnoib(1:3,   88) = (/      10    ,   11     ,  22 /)
     imbou(iimbo) % lnoib(1:3,   89) = (/       3    ,   13     ,   6 /)
     imbou(iimbo) % lnoib(1:3,   90) = (/      48    ,   27     ,  39 /)
     imbou(iimbo) % lnoib(1:3,   91) = (/      16    ,   34     ,  21 /)
     imbou(iimbo) % lnoib(1:3,   92) = (/       1    ,    2     ,   5 /)
     imbou(iimbo) % lnoib(1:3,   93) = (/      14    ,   15     ,   4 /)
     imbou(iimbo) % lnoib(1:3,   94) = (/      67    ,   59     ,  52 /)
     imbou(iimbo) % lnoib(1:3,   95) = (/      35    ,   51     ,  38 /)
     imbou(iimbo) % lnoib(1:3,   96) = (/      28    ,   29     ,  20 /)
     imbou(iimbo) % lnoib(1:3,   97) = (/      13    ,   14     ,   4 /)
     imbou(iimbo) % lnoib(1:3,   98) = (/      11    ,    1     ,   5 /)
     imbou(iimbo) % lnoib(1:3,   99) = (/      66    ,   67     ,  53 /)
     imbou(iimbo) % lnoib(1:3,  100) = (/      27    ,   28     ,  39 /)
     imbou(iimbo) % lnoib(1:3,  101) = (/      34    ,   35     ,  38 /)
     imbou(iimbo) % lnoib(1:3,  102) = (/      15    ,   16     ,   4 /)
     imbou(iimbo) % lnoib(1:3,  103) = (/      59    ,   48     ,  39 /)
     imbou(iimbo) % lnoib(1:3,  104) = (/       2    ,    3     ,   5 /)
     imbou(iimbo) % lnoib(1:3,  105) = (/      51    ,   61     ,  53 /)
     imbou(iimbo) % lnoib(1:3,  106) = (/      29    ,   10     ,  22 /)
     imbou(iimbo) % lnoib(1:3,  107) = (/      38    ,   51     ,  53 /)
     imbou(iimbo) % lnoib(1:3,  108) = (/      21    ,   34     ,  38 /)
     imbou(iimbo) % lnoib(1:3,  109) = (/       6    ,   13     ,   4 /)
     imbou(iimbo) % lnoib(1:3,  110) = (/      22    ,   11     ,   5 /)
     imbou(iimbo) % lnoib(1:3,  111) = (/      20    ,   29     ,  22 /)
     imbou(iimbo) % lnoib(1:3,  112) = (/      52    ,   59     ,  39 /)
     imbou(iimbo) % lnoib(1:3,  113) = (/      16    ,   21     ,   4 /)
     imbou(iimbo) % lnoib(1:3,  114) = (/       3    ,    6     ,   5 /)
     imbou(iimbo) % lnoib(1:3,  115) = (/      28    ,   20     ,  39 /)
     imbou(iimbo) % lnoib(1:3,  116) = (/      67    ,   52     ,  53 /)
     imbou(iimbo) % lnoib(1:3,  117) = (/      39    ,   20     ,  37 /)
     imbou(iimbo) % lnoib(1:3,  118) = (/      39    ,   37     ,  52 /)
     imbou(iimbo) % lnoib(1:3,  119) = (/      37    ,   20     ,  22 /)
     imbou(iimbo) % lnoib(1:3,  120) = (/      38    ,   53     ,  36 /)
     imbou(iimbo) % lnoib(1:3,  121) = (/      36    ,   53     ,  52 /)
     imbou(iimbo) % lnoib(1:3,  122) = (/      38    ,   36     ,  21 /)
     imbou(iimbo) % lnoib(1:3,  123) = (/       4    ,   21     ,  19 /)
     imbou(iimbo) % lnoib(1:3,  124) = (/      19    ,   21     ,  36 /)
     imbou(iimbo) % lnoib(1:3,  125) = (/      19    ,   36     ,  37 /)
     imbou(iimbo) % lnoib(1:3,  126) = (/       4    ,   19     ,   6 /)
     imbou(iimbo) % lnoib(1:3,  127) = (/      22    ,    5     ,  18 /)
     imbou(iimbo) % lnoib(1:3,  128) = (/      18    ,    5     ,   6 /)
     imbou(iimbo) % lnoib(1:3,  129) = (/      22    ,   18     ,  37 /)
     imbou(iimbo) % lnoib(1:3,  130) = (/      37    ,   18     ,  19 /)
     imbou(iimbo) % lnoib(1:3,  131) = (/      37    ,   36     ,  52 /)
     imbou(iimbo) % lnoib(1:3,  132) = (/      19    ,   18     ,   6 /)    

  case ( 3_ip )
     !
     ! Model parameters
     !
     !
     ! RPARA(    1) =  A ... Rotor area
     !      (    2) =  w ... Rotor width
     !      ( 3: 5) =  x ... Rotor center
     !      ( 6: 8) =  n ... Rotor direction
     !      (    9) =  u ... Upstream velocity (not read but computed)
     !      (10:12) = xu ... Upstream position to compute u (not read but computed) 
     ! IPARA(    1) = Ct ... Thrust model 
     !
     print*,'a'
     call ecoute('MODEL')
     do while ( words(1) /= 'ENDMO' )
        if( words(1) == 'AREA ' ) then
           imbou(iimbo) % rpara(1) = getrea('AREA ',1.0_rp,'#Rotor area [m2]')
        else if( words(1) == 'RADIU' ) then
           r0 = getrea('RADIU',1.0_rp,'#Rotor area [m2]')
           imbou(iimbo) % rpara(1) = pi*r0*r0
        else if( words(1) == 'WIDTH' ) then
           imbou(iimbo) % rpara(2) = getrea('WIDTH',1.0_rp,'#Rotor width [m]')
        else if( words(1) == 'POSIT' ) then
           do idime = 1,ndime
              imbou(iimbo) % rpara(2+idime) = param(idime)
           end do
        else if( words(1) == 'DIREC' ) then
           imbou(iimbo) % rpara(6) = param(1)
        else if( words(1) == 'UPSTR' ) then
           imbou(iimbo) % rpara(7) = param(1)
        !else if( words(1) == 'DIREC' ) then
        !   do ipara = 1,nnpar
        !      imbou(iimbo) % rpara(5+ipara) = param(ipara)
        !   end do
        end if
        call ecoute('MODEL')           
     end do
     !
     ! Scale geometry
     !
     r0 = 0.5_rp
     A0 = pi*r0*r0
     A1 = imbou(iimbo) % rpara(1)
     r1 = sqrt(A1/pi)     
     d1 = 2.0_rp*r1
     r1 = r1/r0
     ww = imbou(iimbo) % rpara(2)
     do ipoib = 1,imbou(iimbo) % npoib
        do idime = 1,ndime
           imbou(iimbo) % cooib(idime,ipoib) = r1 * imbou(iimbo) % cooib(idime,ipoib)
        end do
        imbou(iimbo) % cooib(3,ipoib) = imbou(iimbo) % cooib(3,ipoib) / r1 * ww
     end do
     !
     ! Rotate wrt y 90 degrees by default
     !
     rotat      =  0.0_rp
     coso       =  0.0_rp
     sino       =  1.0_rp
     rotat(2,2) =  1.0_rp
     rotat(1,1) =  coso
     rotat(1,3) =  sino
     rotat(3,1) = -sino
     rotat(3,3) =  coso
     do ipoib = 1,imbou(iimbo) % npoib
        xx(1) = imbou(iimbo) % cooib(1,ipoib)
        xx(2) = imbou(iimbo) % cooib(2,ipoib)
        xx(3) = imbou(iimbo) % cooib(3,ipoib)
        call mbvatb(imbou(iimbo) % cooib(1,ipoib),rotat,xx,ndime,ndime)
     end do
     !
     ! Rotate geometry wrt z
     !
     rotat      =  0.0_rp
     theta      =  0.5_rp*pi
     coso       =  cos(theta)
     sino       =  sin(theta)
     rotat(1,1) =  coso
     rotat(1,2) = -sino
     rotat(2,1) =  sino
     rotat(2,2) =  coso
     rotat(3,3) =  1.0_rp
     do ipoib = 1,imbou(iimbo) % npoib
        xx(1) = imbou(iimbo) % cooib(1,ipoib)
        xx(2) = imbou(iimbo) % cooib(2,ipoib)
        xx(3) = imbou(iimbo) % cooib(3,ipoib)
        call mbvatb(imbou(iimbo) % cooib(1,ipoib),rotat,xx,ndime,ndime)
     end do
     !
     ! Rotate geometry wrt z: 0 angle: rotor is facing north
     ! Clockwise rotation: positive angle, hour is increasing
     !
     !                 theta
      !                 -->
     ! | y             | /
     ! |               |/
     ! |           +-------+
     ! |           |       |
     ! |           +-------+
     ! | z
     ! o-------------------------------> x
     !
     rotat      =  0.0_rp
     theta      =  imbou(iimbo) % rpara(6)/180.0_rp*pi
     coso       =  cos(theta)
     sino       =  sin(theta)
     rotat(1,1) =  coso
     rotat(1,2) = -sino
     rotat(2,1) =  sino
     rotat(2,2) =  coso
     rotat(3,3) =  1.0_rp
     do ipoib = 1,imbou(iimbo) % npoib
        xx(1) = imbou(iimbo) % cooib(1,ipoib)
        xx(2) = imbou(iimbo) % cooib(2,ipoib)
        xx(3) = imbou(iimbo) % cooib(3,ipoib)
        call mbvatb(imbou(iimbo) % cooib(1,ipoib),rotat,xx,ndime,ndime)
     end do
     !
     ! Place geometry
     !
     do ipoib = 1,imbou(iimbo) % npoib
        do idime = 1,ndime
           imbou(iimbo) % cooib(idime,ipoib) = imbou(iimbo) % cooib(idime,ipoib) + imbou(iimbo) % rpara(2+idime)
        end do
     end do
     !
     ! Upstream point
     !
     xx  = 0.0_rp     
     nx1 = sin(theta+pi)
     ny1 = cos(theta+pi)
     do idime = 1,ndime
        xx(idime) = imbou(iimbo) % rpara(2+idime)
     end do
     xx(1) = xx(1) + imbou(iimbo) % rpara(7) * nx1
     xx(2) = xx(2) + imbou(iimbo) % rpara(7) * ny1
     !
     ! Postprocess
     !
     open(unit=90,file='caca.post.msh',status='unknown')
     write(90,*) 'MESH cavtri03_TRI03 dimension 3 Elemtype Triangle Nnode  3'
     write(90,*) 'coordinates'
     do ipoib = 1,imbou(iimbo) % npoib
        write(90,'(i4,3(1x,e12.6))') ipoib,(imbou(iimbo) % cooib(idime,ipoib),idime=1,ndime)
     end do
     write(90,'(i4,3(1x,e12.6))') imbou(iimbo) % npoib+1,xx(1),xx(2),xx(3)
     write(90,'(i4,3(1x,e12.6))') imbou(iimbo) % npoib+2,xx(1),xx(2),xx(3)
     write(90,'(i4,3(1x,e12.6))') imbou(iimbo) % npoib+3,xx(1),xx(2),xx(3)
     write(90,*) 'end coordinates'
     write(90,*) 'elements'
     do iboib = 1,imbou(iimbo) % nboib
        write(90,'(5(1x,i4))') iboib,(imbou(iimbo) % lnoib(idime,iboib),idime=1,3),1
     end do
     ii = imbou(iimbo) % npoib
     write(90,'(5(1x,i4))') imbou(iimbo) % nboib+1,ii+1,ii+2,ii+3,2
     write(90,*) 'end elements'
     close(90)

  end select

end subroutine ibm_model_wind_turbine
