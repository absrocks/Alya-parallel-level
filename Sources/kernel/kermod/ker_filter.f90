subroutine ker_filter()
  !------------------------------------------------------------------------
  !****f* master/ker_filter
  ! NAME 
  !    ker_filter
  ! DESCRIPTION
  !    Filter for postprocess
  ! USES
  ! USED BY
  !    nsi_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  implicit none
  integer(ip) :: ipoin

  select case ( kfl_filte )

  case (1_ip) 
     !
     ! Coordinates box nose + trachea(new! with bronche)
     !
     !
     do ipoin = 1,npoin
        if (((coord(1,ipoin)<0.116_rp .and. coord(1,ipoin)>0.0752_rp).and.   &
             (coord(2,ipoin)<0.105_rp .and. coord(2,ipoin)>0.0146_rp).and.  &
             (coord(3,ipoin)<0.0539_rp .and. coord(3,ipoin)>0.00114_rp))    &
             .or.                                                           &
             ((coord(1,ipoin)<0.595_rp .and. coord(1,ipoin)>-0.396_rp).and. &
             (coord(2,ipoin)<0.326_rp .and. coord(2,ipoin)>0.105_rp).and.    &
             (coord(3,ipoin)<0.041_rp .and. coord(3,ipoin)>-0.488_rp)))     &
             gefil(ipoin) = 1
     end do

  case (2_ip)
     !
     ! Coordinates box nose
     !
     !
     do ipoin = 1,npoin
        if (((coord(1,ipoin)<0.116_rp .and. coord(1,ipoin)>0.0752_rp).and.   &
             (coord(2,ipoin)<0.105_rp .and. coord(2,ipoin)>0.0146_rp).and.  &
             (coord(3,ipoin)<0.0539_rp .and. coord(3,ipoin)>0.00114_rp)))    &
             gefil(ipoin) = 1
     end do
     
  case (3_ip)
     !
     ! Plane
     !
     call platri()

  case (4_ip)
     !
     ! Plane + boundary
     !
     call platri()

     do ipoin = 1,npoin
        if(lpoty(ipoin)>0_ip) gefil(ipoin) = 1
     end do

  case (5_ip)
     !
     ! several plans for the nose
     !
     !call noseplan (0.99760002_rp, -0.10000000_rp, 6.98999986E-02_rp,-9.90286469E-02_rp)
     !call noseplan (0.99760002_rp, -0.10000000_rp, 6.98999986E-02_rp,-9.60358456E-02_rp)
     !call noseplan (0.99760002_rp, -0.10000000_rp, 6.98999986E-02_rp,-9.50382426E-02_rp)
     !call noseplan (0.99760002_rp, -0.10000000_rp, 6.98999986E-02_rp,-9.40406471E-02_rp)
     !call noseplan (0.99760002_rp, -0.10000000_rp, 6.98999986E-02_rp,-9.30430442E-02_rp)
     !call noseplan (0.99760002_rp, -0.10000000_rp, 6.98999986E-02_rp,-9.00502503E-02_rp)
     !call noseplan (0.99260002_rp, 5.00000007E-02_rp, -0.20000000_rp,-8.34573060E-02_rp)
     !call noseplan (0.99260002_rp, 5.00000007E-02_rp, -0.20000000_rp,-8.64351019E-02_rp)
     !call noseplan (0.99260002_rp, 5.00000007E-02_rp, -0.20000000_rp,-8.74277055E-02_rp)
     !call noseplan (0.99260002_rp, 5.00000007E-02_rp, -0.20000000_rp,-8.84203017E-02_rp)
     !call noseplan (0.99260002_rp, 5.00000007E-02_rp, -0.20000000_rp,-8.94128978E-02_rp)
     !call noseplan (0.99260002_rp, 5.00000007E-02_rp, -0.20000000_rp,-9.23907012E-02_rp)
     !call noseplan ( 0.0000000_rp, 0.34490001_rp, 0.93860000_rp,-2.46693268E-02_rp)
     !call noseplan ( 0.0000000_rp, 0.70709997_rp, -0.70709997_rp,-8.11750814E-02_rp)
     !call noseplan ( 0.0000000_rp, 0.81650001_rp, 0.57730001_rp,-4.31093611E-02_rp)
     !call noseplan ( 0.0000000_rp, 0.99500000_rp, 9.95000005E-02_rp,-5.38792498E-02_rp)
     !call noseplan ( 0.0000000_rp, 1.0000000_rp, 0.0000000_rp, -6.15999997E-02_rp)
     !call noseplan ( 0.0000000_rp, 1.0000000_rp, 0.0000000_rp, -6.66000023E-02_rp)
     !call noseplan ( 0.0000000_rp, 1.0000000_rp, 0.0000000_rp, -7.15999976E-02_rp)
     !call noseplan ( 0.0000000_rp, 1.0000000_rp, 0.0000000_rp, -9.48000029E-02_rp)
     !call noseplan (-5.00000024E-04_rp, 0.97979999_rp, -0.19990000_rp,-0.10485945_rp)
     !call noseplan ( 1.0_rp, 0.0_rp, 0.0_rp, -9.89999995E-02_rp)
     !call noseplan ( 0.0_rp, 0.0_rp, 1.0_rp, 9.00000036E-02_rp)
     !call noseplan ( 0.0_rp, 0.0_rp, 1.0_rp, 6.49999976E-02_rp)
     !call noseplan ( 0.0_rp, 0.0_rp, 1.0_rp, 5.56841157E-02_rp)
     !call noseplan ( 0.0_rp, 0.0_rp, 1.0_rp, 9.16194767E-02_rp)
     !
     !add some plan 
     !
     call noseplan ( 1.0_rp, 0.0_rp, 0.0_rp, -0.1_rp)
     call noseplan ( 0.97912198_rp, 4.93011065E-02_rp, -0.19720443_rp, -8.52620006E-02_rp)
     call noseplan ( 0.0_rp, 1.0_rp, 0.0_rp, -0.05_rp)
     call noseplan ( 0.0_rp, 1.0_rp, 0.0_rp, -0.08_rp)
     call noseplan ( 0.0_rp, 1.0_rp, 0.0_rp, -0.10_rp)
     call noseplan ( 0.0_rp, 1.0_rp, 0.0_rp, -0.11_rp)    
     call noseplan ( 0.0_rp, -0.9_rp, 0.4_rp, 9.84368473E-02_rp)
     call noseplan ( 0.0_rp, -0.8_rp, 0.5_rp, 9.43255052E-02_rp)
     call noseplan ( 0.0_rp, 0.0_rp, 1.0_rp, -4.00139857E-03_rp)
     call noseplan ( 0.0_rp, 0.0_rp, 1.0_rp, 1.0E-02_rp)
     call noseplan ( 0.0_rp, 0.0_rp, 1.0_rp, 3.0E-02_rp)
     call noseplan ( 0.0_rp, 0.0_rp, 1.0_rp, 5.0E-02_rp)
     call noseplan ( 0.0_rp, 0.0_rp, 1.0_rp, 7.0E-02_rp)
     call noseplan ( 0.0_rp, 0.0_rp, 1.0_rp, 9.0E-02_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.13872665_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.15_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.17_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.19_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.21_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.23_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.25_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.26_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.27_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.28_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.29_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.30_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.32_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.34_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.35_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.36_rp)
     call noseplan ( -3.86460312E-02_rp, -0.42482808_rp, 0.90444875_rp, 0.37_rp)
     call noseplanclean()

  case (6_ip)
     !
     ! Coordinates box trachea reduced for POD 
     !
     !
     do ipoin = 1,npoin
        if (((coord(1,ipoin)<0.113_rp .and. coord(1,ipoin)>0.0863_rp).and.   &
             (coord(2,ipoin)<0.149_rp .and. coord(2,ipoin)>0.122_rp).and.  &
             (coord(3,ipoin)<-0.0542_rp .and. coord(3,ipoin)>-0.109_rp)))    &
             gefil(ipoin) = 1
     end do
     

  case (7_ip)

  case default

     call runend('KER_FILTER: THIS FILTER WAS NOT CODED')

  end select

end subroutine ker_filter
