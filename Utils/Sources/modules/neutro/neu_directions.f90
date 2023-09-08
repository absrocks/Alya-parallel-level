!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_directions.f90
!> @date    31/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Directions
!> @details Directions
!> @} 
!-----------------------------------------------------------------------

subroutine neu_directions()

  use def_kintyp, only : ip,rp
  use def_master, only : INOTMASTER,momod,modul,INOTSLAVE
  use def_domain, only : ndime
  use mod_maths,  only : maths_angle_of_a_vector
  use def_neutro, only : kfl_icosa_neu,direc_neu
  use def_neutro, only : kfl_snord_neu,num_directions_neu
  implicit none
  integer(ip) :: idire
  real(rp)    :: angl1,angl2

  if( kfl_icosa_neu == 1 )  then 
     !
     ! ICOSAHEDRON DIRECTION
     !
     call neu_ordico()  
     
  else if( kfl_snord_neu /= 0 ) then
     !
     ! SN8 DIRECTIONS
     !
     call neu_snordi()
     
  end if
  !
  ! Write direction in output file
  !
  if( INOTSLAVE ) then
     do idire = 1,num_directions_neu 
        call maths_angle_of_a_vector(ndime,direc_neu(1:ndime,idire),angl1,angl2)
        write(momod(modul)%lun_outpu,*) 'Direction: ',idire,' angle= ',angl1
     end do
     flush(momod(modul)%lun_outpu)
  end if

end subroutine neu_directions

subroutine neu_snordi()  

  !-----------------------------------------------------------------------
  ! constructs SN ordinates
  !-----------------------------------------------------------------------

  use def_kintyp, only : ip,rp
  use def_domain, only : ndime
  use def_neutro, only : kfl_snord_neu
  use def_neutro, only : num_directions_neu
  use def_neutro, only : direc_neu
  use def_neutro, only : weigd_neu
  implicit none

  real(rp)    :: dirpr(6) 
  integer(ip) :: idire,idime

  select case ( kfl_snord_neu )

  case ( 4_ip )

     !-------------------------------------------------------------------
     !
     ! SN4
     !
     !-------------------------------------------------------------------
     !
     ! cosen directions
     !
     dirpr(1) = 0.2958759_rp
     dirpr(2) = 0.9082483_rp
     !
     !   first octant
     !
     direc_neu (1, 1) = dirpr(1)
     direc_neu (1, 2) = dirpr(1)
     direc_neu (1, 3) = dirpr(2)

     direc_neu (2, 1) = dirpr(1)
     direc_neu (2, 2) = dirpr(2)
     direc_neu (2, 3) = dirpr(1)

     direc_neu (3, 1) = dirpr(2)
     direc_neu (3, 2) = dirpr(1)
     direc_neu (3, 3) = dirpr(1)


     weigd_neu(1)=0.523598775598299_rp
     weigd_neu(2)=0.523598775598299_rp
     weigd_neu(3)=0.523598775598299_rp
     !
     ! octant II z>0, y>0, x<0 
     !
     do idire =4, 6
        direc_neu(1,idire)= - direc_neu(1,idire-3)
        do idime =2,3
           direc_neu(idime,idire)= direc_neu(idime,idire-3)
        end do
        weigd_neu(idire) = weigd_neu(idire-3)
     end do
     !
     ! octants III and IV z>0, y<0, x>0 and x<0
     !
     do idire =7, 12
        direc_neu(2,idire) = -direc_neu(2,idire-6)
        direc_neu(1,idire) =  direc_neu(1,idire-6)
        direc_neu(3,idire) =  direc_neu(3,idire-6)    
        weigd_neu(idire)   =  weigd_neu(idire-6)
     end do
     !
     ! octants with z<0
     !
     if( ndime == 3 ) then
        do idire =13, 24
           direc_neu(2,idire) =  direc_neu(2,idire-12)
           direc_neu(1,idire) =  direc_neu(1,idire-12)
           direc_neu(3,idire) = -direc_neu(3,idire-12)    
           weigd_neu(idire)   =  weigd_neu(idire-12)
        end do
     else if( ndime == 2 ) then  
        do idire=1, num_directions_neu
           weigd_neu(idire) = 2.0_rp * weigd_neu(idire)
        end do
     end if

  case ( 8_ip ) 

     !-------------------------------------------------------------------
     !
     ! S8-approximation
     !
     !-------------------------------------------------------------------
     !
     ! cosen directions
     !
     dirpr(1) = 0.1422555_rp
     dirpr(2) = 0.5773503_rp
     dirpr(3) = 0.8040087_rp
     dirpr(4) = 0.9795543_rp
     !
     ! first octant
     !
     do idire=1, 4
        direc_neu(1,idire) = dirpr(1)   
     end do
     do idire=5, 7
        direc_neu(1,idire) = dirpr(2)   
     end do
     do idire=8,9
        direc_neu(1,idire) = dirpr(3)   
     end do
     direc_neu(1,10) = dirpr(4)   

     do idire=1, 3
        direc_neu(2,idire)    = dirpr(idire)   
        direc_neu(2,idire +4) = dirpr(idire)   
        direc_neu(2,idire +7) = dirpr(idire)   
     end do

     direc_neu(2, 4) = dirpr(4)   
     direc_neu(2,10) = dirpr(1) 

     do idire=1, 4
        direc_neu(3,idire)    = dirpr(5-idire)
        direc_neu(3,idire +3) = dirpr(5-idire)
     end do
     direc_neu(3,8 )  = dirpr(2)
     direc_neu(3,9 )  = dirpr(1)
     direc_neu(3,10)  = dirpr(1)

     weigd_neu( 1) = 0.1712359_rp
     weigd_neu( 2) = 0.0992284_rp
     weigd_neu( 3) = 0.0992284_rp
     weigd_neu( 4) = 0.1712359_rp
     weigd_neu( 5) = 0.0992284_rp
     weigd_neu( 6) = 0.4617179_rp
     weigd_neu( 7) = 0.0992284_rp
     weigd_neu( 8) = 0.0992284_rp
     weigd_neu( 9) = 0.0992284_rp
     weigd_neu(10) = 0.1712359_rp
     !
     !octant II z>0, y>0, x<0 
     !
     do idire =11, 20
        direc_neu(1,idire)= - direc_neu(1,idire-10)
        do idime =2,3
           direc_neu(idime,idire)= direc_neu(idime,idire-10)
        end do
        weigd_neu(idire) = weigd_neu(idire-10)
     end do
     !
     !octants III and IV z>0, y<0, x>0 and x<0
     !
     do idire =21, 40
        direc_neu(2,idire)= - direc_neu(2,idire-20)
        direc_neu(1,idire)= direc_neu(1,idire-20)
        direc_neu(3,idire)= direc_neu(3,idire-20)    
        weigd_neu(idire) = weigd_neu(idire-20)
     end do
     !
     !octants with z<0
     !
     if( ndime == 3 ) then
        do idire = 41,80
           direc_neu(2,idire) =  direc_neu(2,idire-40)
           direc_neu(1,idire) =  direc_neu(1,idire-40)
           direc_neu(3,idire) = -direc_neu(3,idire-40)    
           weigd_neu(idire)       =  weigd_neu(idire-40)
        end do
     else if( ndime == 2 ) then  
        do idire = 1,num_directions_neu
           weigd_neu(idire) = 2*weigd_neu(idire)
        end do
     end if

  case ( 6_ip )

     !-------------------------------------------------------------------
     !  
     ! SN6: S6-approximation
     !
     !-------------------------------------------------------------------
     !
     ! cosen directions
     !
     dirpr(1) = 0.1838670_rp
     dirpr(2) = 0.6950514_rp
     dirpr(3) = 0.9656013_rp
     !
     ! first octant
     !
     do idire=1, 3
        direc_neu(1,idire) = dirpr(1)   
     end do
     do idire=4, 5
        direc_neu(1,idire) = dirpr(2)   
     end do

     direc_neu(1,6) = dirpr(3)   

     do idire=1, 3
        direc_neu(2,idire)    = dirpr(idire)   
        direc_neu(2,idire +3) = dirpr(idire)   

     end do

     direc_neu(2, 6) = dirpr(1)   

     do idire=1, 3
        direc_neu(3,idire)    = dirpr(4-idire)   
     end do

     direc_neu(3,4) = dirpr(2)
     direc_neu(3,5) = dirpr(1)
     direc_neu(3,6) = dirpr(1)

     weigd_neu(1) = 0.1609517_rp
     weigd_neu(2) = 0.3626469_rp
     weigd_neu(3) = 0.1609517_rp
     weigd_neu(4) = 0.3626469_rp
     weigd_neu(5) = 0.3626469_rp
     weigd_neu(6) = 0.1609517_rp
     !
     ! Octant II z>0, y>0, x<0 
     !
     do idire =7, 12
        direc_neu(1,idire)= - direc_neu(1,idire-6)
        do idime =2,3
           direc_neu(idime,idire)= direc_neu(idime,idire-6)
        end do
        weigd_neu(idire) = weigd_neu(idire-6)
     end do
     !
     ! Octants III and IV z>0, y<0, x>0 and x<0
     !
     do idire =13, 24
        direc_neu(2,idire)= - direc_neu(2,idire-12)
        direc_neu(1,idire)= direc_neu(1,idire-12)
        direc_neu(3,idire)= direc_neu(3,idire-12)    
        weigd_neu(idire) = weigd_neu(idire-12)
     end do
     !
     ! Octants with z<0
     !
     if( ndime == 3 ) then
        do idire =25, 48
           direc_neu(2,idire)= direc_neu(2,idire-24)
           direc_neu(1,idire)= direc_neu(1,idire-24)
           direc_neu(3,idire)= -direc_neu(3,idire-24)    
           weigd_neu(idire) = weigd_neu(idire-24)
        end do
     else if( ndime == 2 ) then  
        do idire=1, num_directions_neu
           weigd_neu(idire)=2*weigd_neu(idire)
        end do
     end if

  case ( 10_ip )   

     !-------------------------------------------------------------------
     !
     ! SN10, S10-approximation
     !
     !-------------------------------------------------------------------
     !
     !cosen directions
     !
     dirpr(1) = 0.13727193312_rp
     dirpr(2) = 0.50468891002_rp
     dirpr(3) = 0.70041288408_rp
     dirpr(4) = 0.85231773445_rp
     dirpr(5) = 0.98097544961_rp
     !
     !   first octant
     !
     do idire = 1,3
        do idime = 1,3
           direc_neu(idime,idire) = dirpr(1)   
        end do
        direc_neu(idire, idire)= dirpr(5)
     end do
     do idire = 4,6
        do idime = 1,3
           direc_neu(idime,idire) = dirpr(2)   
        end do
        direc_neu(idire-3, idire)= dirpr(3)
     end do
     do idire = 13, 15
        do idime = 1,3
           direc_neu(idime,idire) = dirpr(3)   
        end do
        direc_neu(idire-12, idire)= dirpr(1)
     end do
     do idire = 7,9
        direc_neu(idire-6, idire)   = dirpr(1)
        direc_neu(idire-6, idire+3) = dirpr(1)
     end do
     do idire = 7,8
        direc_neu(idire-5, idire)   = dirpr(2)
        direc_neu(idire-5, idire+3) = dirpr(4)
     end do
     direc_neu(1, 9)  = dirpr(2)
     direc_neu(1, 12) = dirpr(4)
     do idire = 8,9
        direc_neu(idire-7, idire)   = dirpr(4)
        direc_neu(idire-7, idire+3) = dirpr(2)
     end do
     direc_neu(3, 7)  = dirpr(4)
     direc_neu(3, 10) = dirpr(2)

     do idire =1,3
        weigd_neu(idire) = 0.0944411600_rp
     end do
     do idire =4,6
        weigd_neu(idire) = 0.1149971656_rp
     end do
     do idire =7,12
        weigd_neu(idire) = 0.1483951164_rp
     end do
     do idire =13,15
        weigd_neu(idire) = 0.0173702170_rp
     end do
     !
     ! octant II z>0, y>0, x<0 
     !
     do idire =16, 30
        direc_neu(1,idire)= - direc_neu(1,idire-15)
        do idime =2,3
           direc_neu(idime,idire)= direc_neu(idime,idire-15)
        end do
        weigd_neu(idire) = weigd_neu(idire-15)
     end do
     !
     ! octants III and IV z>0, y<0, x>0 and x<0
     !
     do idire =31, 60 
        direc_neu(2,idire) = -direc_neu(2,idire-30)
        direc_neu(1,idire) =  direc_neu(1,idire-30)
        direc_neu(3,idire) =  direc_neu(3,idire-30)    
        weigd_neu(idire)   =  weigd_neu(idire-30)
     end do
     !
     ! octants with z<0
     !
     if( ndime == 3 ) then
        do idire =61, 120
           direc_neu(2,idire) =  direc_neu(2,idire-60)
           direc_neu(1,idire) =  direc_neu(1,idire-60)
           direc_neu(3,idire) = -direc_neu(3,idire-60)    
           weigd_neu(idire)   =  weigd_neu(idire-60)
        end do
     else if( ndime == 2 ) then  
        do idire=1, num_directions_neu
           weigd_neu(idire) = 2.0_rp * weigd_neu(idire)
        end do
     end if

  case ( 12_ip )

     !-------------------------------------------------------------------
     !
     ! SN12: S12-approximation
     !
     !-------------------------------------------------------------------
     !
     ! cosen directions
     !
     dirpr(1) = 0.12028966644554817_rp
     dirpr(2) = 0.45363844804142206_rp
     dirpr(3) = 0.63016353371905381_rp
     dirpr(4) = 0.76708820673840594_rp
     dirpr(5) = 0.88303032485016963_rp
     dirpr(6) = 0.98542416871763827_rp
     !
     ! first octant
     !
     do idire=1, 3
        do idime =1,3
           direc_neu(idime,idire) = dirpr(1)   
        end do
        direc_neu(idire, idire)= dirpr(6)
     end do
     do idire=4, 6
        do idime =1,3
           direc_neu(idime,idire) = dirpr(2)   
        end do
        direc_neu(idire-3, idire)= dirpr(4)
     end do
     do idire=7, 9
        do idime =1,3
           direc_neu(idime,idire) = dirpr(3)   
        end do
        direc_neu(idire-6, idire)= dirpr(2)
     end do

     do idire=10, 12
        direc_neu(idire-9, idire)   = dirpr(1)
        direc_neu(idire-9, idire+3) = dirpr(1)
     end do
     do idire=10, 11
        direc_neu(idire-8, idire)   = dirpr(2)
        direc_neu(idire-8, idire+3) = dirpr(5)
     end do
     direc_neu(1, 12)  = dirpr(2)
     direc_neu(1, 15) = dirpr(5)
     do idire=11, 12
        direc_neu(idire-10, idire)   = dirpr(5)
        direc_neu(idire-10, idire+3) = dirpr(2)
     end do
     direc_neu(3, 10)  = dirpr(5)
     direc_neu(3, 13)  = dirpr(2)

     do idire=16, 18
        direc_neu(idire-15, idire)   = dirpr(1)
        direc_neu(idire-15, idire+3) = dirpr(1)
     end do
     do idire=16, 17
        direc_neu(idire-14, idire)   = dirpr(3)
        direc_neu(idire-14, idire+3) = dirpr(4)
     end do
     direc_neu(1, 18)  = dirpr(3)
     direc_neu(1, 21) = dirpr(4)
     do idire=17, 18
        direc_neu(idire-16, idire)   = dirpr(4)
        direc_neu(idire-16, idire+3) = dirpr(3)
     end do
     direc_neu(3, 16)  = dirpr(4)
     direc_neu(3, 19)  = dirpr(3) 


     do idire =1,3
        weigd_neu(idire) = 0.0801404674_rp
     end do
     do idire =4,6
        weigd_neu(idire) = 0.1357133976_rp
     end do
     do idire =7,9
        weigd_neu(idire) = 0.0239769589_rp
     end do
     do idire =10,15
        weigd_neu(idire) = 0.0974977375_rp
     end do
     do idire =16,21
        weigd_neu(idire) = 0.0443862383_rp
     end do
     !
     ! octant II z>0, y>0, x<0 
     !
     do idire =22, 42
        direc_neu(1,idire)= - direc_neu(1,idire-21)
        do idime =2,3
           direc_neu(idime,idire)= direc_neu(idime,idire-21)
        end do
        weigd_neu(idire) = weigd_neu(idire-21)
     end do
     !
     ! octants III and IV z>0, y<0, x>0 and x<0
     !
     do idire =43, 84 
        direc_neu(2,idire)= - direc_neu(2,idire-42)
        direc_neu(1,idire)= direc_neu(1,idire-42)
        direc_neu(3,idire)= direc_neu(3,idire-42)    
        weigd_neu(idire) = weigd_neu(idire-42)
     end do
     !
     ! octants with z<0
     !
     if( ndime == 3 ) then
        do idire = 85, 168
           direc_neu(2,idire)= direc_neu(2,idire-84)
           direc_neu(1,idire)= direc_neu(1,idire-84)
           direc_neu(3,idire)= -direc_neu(3,idire-84)    
           weigd_neu(idire) = weigd_neu(idire-84)
        end do
     else if( ndime == 2 ) then  
        do idire=1, num_directions_neu
           weigd_neu(idire)= 2.0_rp * weigd_neu(idire)
        end do
     end if

  end select


end subroutine neu_snordi

subroutine neu_ordico()
  !-----------------------------------------------------------------------
  !constructs ordinates for an icosahedron
  !defines direc_neu and weigd_neu
  !
  !-----------------------------------------------------------------------

  use def_kintyp, only : ip,rp
  use def_parame, only : pi
  use def_domain, only : ndime
  use def_neutro
  use def_master
  implicit none

  integer(ip) :: istat, idire, idime, ipoin, i, kfl_normx_rad
  real(rp)    :: vertp(3,12), auxve(12), vertr(3, 12),  auxdi(3,40)
  real(rp)    :: goldr, normv, rdumm

  kfl_normx_rad = 0
  !  
  ! construction of icosahedron
  !
  !
  ! ICOSAHEDRON VERTEXS  vertp(idime,ipoin)  ipoin =1, 12
  ! golden ratio
  !
  goldr = 0.5_rp*(1.0_rp+sqrt(5.0_rp))
  goldr = 1.0_rp/goldr
  auxve = 0.0_rp

  do ipoin=1,4
     vertp(1,ipoin)   = 0.0_rp
     vertp(3,ipoin+4) = 0.0_rp
     vertp(2,ipoin+8) = 0.0_rp
  end do
  vertp(2, 1) =  1.0_rp
  vertp(2, 2) =  1.0_rp
  vertp(2, 3) = -1.0_rp
  vertp(2, 4) = -1.0_rp

  vertp(3, 1) =  goldr
  vertp(3, 2) = -goldr
  vertp(3, 3) =  goldr
  vertp(3, 4) = -goldr

  vertp(1, 5) =  1.0_rp
  vertp(1, 6) =  1.0_rp
  vertp(1, 7) = -1.0_rp
  vertp(1, 8) = -1.0_rp

  vertp(2, 5) =  goldr
  vertp(2, 6) = -goldr
  vertp(2, 7) =  goldr
  vertp(2, 8) = -goldr

  vertp(1, 9) =  goldr
  vertp(1,10) =  goldr
  vertp(1,11) = -goldr
  vertp(1,12) = -goldr

  vertp(3, 9) =  1.0_rp
  vertp(3,10) = -1.0_rp
  vertp(3,11) =  1.0_rp
  vertp(3,12) = -1.0_rp  

  do idire =1, 12
     call vecuni(3_ip,vertp(1,idire),normv)  
  end do

  !--------------  
  ! ROTATION MATRIX: Rotation around x vector to have symmetry (casi)respect to xy plane
  !--------------
  !    |1.0           0.0                 0.0            |      
  !    |                                                 |      
  !    |0.0  gr/sqrt(gr*gr+1)    -1.0/sqrt(gr*gr+1)      |
  ! R= |                                                 |
  !    |0.0  1.0/sqrt(gr*gr+1)   gr/sqrt(gr*gr+1)        |
  !    |                                                 |
  !

  rdumm = 1.0_rp/sqrt(goldr*goldr+1.0_rp)

  vertr (1, :) = vertp(1,:)
  vertr (2, :) = rdumm*goldr*vertp(2,:) - rdumm*vertp(3,:)
  vertr (3, :) = rdumm*vertp(2,:) + goldr*rdumm*vertp(3,:)

  if( kfl_normx_rad == 1 ) then
     auxve (:)    = vertr (3,:)
     vertr (3, :) = vertr (2,:) 
     vertr (2, :) = vertr (1,:)
     vertr (1, :) = auxve (:)
  end if
  !
  !
  ! I HAVE 3 VERTEXS PER ICOSAHEDRON FACE
  !
  if( num_directions_neu == 20 ) then

     do idime = 1,3 !direction vectors pointing to the center of icosahedron faces

        direc_neu(idime, 1) = vertp(idime,9)+ vertp(idime,11)+ vertp(idime,1)
        direc_neu(idime, 2) = vertp(idime,9)+ vertp(idime,1)+ vertp(idime,5)
        direc_neu(idime, 3) = vertp(idime,1)+ vertp(idime,5)+ vertp(idime,2)
        direc_neu(idime, 4) = vertp(idime,1)+ vertp(idime,7)+ vertp(idime,2)
        direc_neu(idime, 5) = vertp(idime,1)+ vertp(idime,7)+ vertp(idime,11)
        direc_neu(idime, 6) = vertp(idime,9)+ vertp(idime,11)+ vertp(idime,3)
        direc_neu(idime, 7) = vertp(idime,9)+ vertp(idime,3)+ vertp(idime,6)
        direc_neu(idime, 8) = vertp(idime,3)+ vertp(idime,6)+ vertp(idime,4)
        direc_neu(idime, 9) = vertp(idime,3)+ vertp(idime,4)+ vertp(idime,8)
        direc_neu(idime, 10) = vertp(idime,3)+ vertp(idime,8)+ vertp(idime,11)
        direc_neu(idime, 11) = vertp(idime,10)+ vertp(idime,2)+ vertp(idime,5)
        direc_neu(idime, 12) = vertp(idime,10)+ vertp(idime,6)+ vertp(idime,5)
        direc_neu(idime, 13) = vertp(idime,10)+ vertp(idime,4)+ vertp(idime,6)
        direc_neu(idime, 14) = vertp(idime,10)+ vertp(idime,2)+ vertp(idime,12)
        direc_neu(idime, 15) = vertp(idime,10)+ vertp(idime,4)+ vertp(idime,12)
        direc_neu(idime, 16) = vertp(idime,9)+ vertp(idime,6)+ vertp(idime,5)
        direc_neu(idime, 17) = vertp(idime,11)+ vertp(idime,8)+ vertp(idime,7)
        direc_neu(idime, 18) = vertp(idime,8)+ vertp(idime,4)+ vertp(idime,12)
        direc_neu(idime, 19) = vertp(idime,12)+ vertp(idime,7)+ vertp(idime,2)
        direc_neu(idime, 20) = vertp(idime,12)+ vertp(idime,7)+ vertp(idime,8)     
     end do
     do idire =1, num_directions_neu
        call vecuni(3_ip,direc_neu(:,idire),normv) 
     end do

  else

     i =0
     call neu_icodir(vertp(1,9), vertp(1,11), vertp(1,1), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,9), vertp(1,1), vertp(1,5), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,1), vertp(1,5), vertp(1,2), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,1), vertp(1,7), vertp(1,2), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,1), vertp(1,7), vertp(1,11), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,9), vertp(1,11), vertp(1,3), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,9), vertp(1,3), vertp(1,6), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,3), vertp(1,6), vertp(1,4), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,3), vertp(1,4), vertp(1,8), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,3), vertp(1,8), vertp(1,11), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,10), vertp(1,2), vertp(1,5), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,10), vertp(1,6), vertp(1,5), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,10), vertp(1,4), vertp(1,6), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,10), vertp(1,2), vertp(1,12), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,10), vertp(1,4), vertp(1,12), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,9), vertp(1,6), vertp(1,5), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,11), vertp(1,8), vertp(1,7), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,8), vertp(1,4), vertp(1,12), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,12), vertp(1,7), vertp(1,2), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))
     i= i +1
     call neu_icodir(vertp(1,12), vertp(1,7), vertp(1,8), direc_neu(1,i+1), direc_neu(1,20+3*i+1), direc_neu(1,20+ 3*i +2), direc_neu(1,20 + 3*i +3))

  end if
  !
  ! STORES ONLY HALF DIRECTIONS FOR 2 DIMENSIONAL CASE (z>0) 
  !
  rdumm = 0.0_rp
  if( ndime == 2 .and. num_directions_neu /= 20 ) then !only directions with z component gt 0.0, save half directions 
     call runend('NEU_DIRECTIONS: NOT CODED')
  end if

  do idire = 1,num_directions_neu
     weigd_neu(idire) = 4-0_rp*pi/real(num_directions_neu,rp)
  end do

end subroutine neu_ordico

subroutine neu_icodir(vert1, vert2, vert3, direm, dire1, dire2, dire3)
  !-----------------------------------------------------------------------
  !subroutine that obtains 4 directions in an icosahedron face
  !vert1, vert2, vert3 are the 3 vertex coordinates
  !direm will be the central direction
  !dire1, dire2, dire3 will be the other 3        /\
  !                                              /--\
  !                                             / \/ \
  !                                            /------\ 
  !-----------------------------------------------------------------------
  use def_kintyp, only : ip,rp
  implicit none
  real(rp), intent(in)  :: vert1(3), vert2(3), vert3(3)
  real(rp), intent(out) :: direm(3), dire1(3), dire2(3), dire3(3)
  real(rp)              :: vert12(3), vert13(3), vert23(3)
  real(rp)              :: normv
  integer(ip)           :: idime
  !    do idime = 1, 3
  !     direm(idime) = vert1(idime)+ vert2(idime)+ vert3(idime)
  !     vert12(idime) = vert1(idime) + vert2(idime)
  !     vert13(idime) = vert1(idime) + vert3(idime)
  !     vert23(idime) = vert2(idime) + vert3(idime)
  !    end do
  do idime = 1, 3
     direm(idime)  = vert1(idime)+ vert2(idime)+ vert3(idime)
     vert12(idime) = vert1(idime) + vert2(idime)
     vert13(idime) = vert1(idime) + vert3(idime)
     vert23(idime) = vert2(idime) + vert3(idime)
  end do

  call vecuni(3_ip,direm, normv)  
  call vecuni(3_ip,vert12,normv)  
  call vecuni(3_ip,vert13,normv)  
  call vecuni(3_ip,vert23,normv) 

  do idime = 1,3
     dire1(idime) = vert1(idime) + vert12(idime)+ vert13(idime)
     dire2(idime) = vert2(idime) + vert12(idime)+ vert23(idime)
     dire3(idime) = vert3(idime) + vert13(idime)+ vert23(idime)
  end do

  call vecuni(3_ip,dire1,normv)  
  call vecuni(3_ip,dire2,normv)  
  call vecuni(3_ip,dire3,normv)  

end subroutine neu_icodir
