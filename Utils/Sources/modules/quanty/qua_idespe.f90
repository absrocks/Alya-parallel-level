subroutine qua_especiename(nelec,ntipo,wprob1,wprob2)
  use mod_outfor, only : outfor
  use def_master
  implicit none
  integer(ip), intent(in)  :: nelec,ntipo
  character(50), intent(out) :: wprob1,wprob2


  SELECT CASE (nelec)
  CASE (1)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'H.CPI'
  CASE (2)
     wprob1 = 'Helium'
     if(ntipo==0) wprob2 = 'HE.CPI'
  CASE (3)
     wprob1 = 'Lytium'
     if(ntipo==0)  wprob2 = 'LI.CPI'
  CASE (4)
     wprob1 = 'Berilium'
     if(ntipo==0) wprob2 = 'BE.CPI'
  CASE (5)
     wprob1 = 'Boro'
     if(ntipo==0) wprob2 = 'B.CPI'
  CASE (6)
     wprob1 = 'Carbono'
     if(ntipo==0) wprob2 = 'C.CPI'
  CASE (7)
     wprob1 = 'Nitrogeno'
     if(ntipo==0) wprob2 = 'N.CPI'
  CASE (8)
     wprob1 = 'Oxigeno'
     if(ntipo==0) wprob2 = 'O.CPI'
  CASE (9)
     wprob1 = 'Fluor'
     if(ntipo==0) wprob2 = 'F.CPI'
  CASE (10)
     wprob1 = 'Neon'
     if(ntipo==0) wprob2 = 'NE.CPI'
  CASE (11)
     wprob1 = 'Sodium'
     if(ntipo==0) wprob2 = 'NA.CPI'
  CASE (12)
     wprob1 = 'Magnesium'
     if(ntipo==0) wprob2 = 'MG.CPI'
  CASE (13)
     wprob1 = 'Aluminium'
     if(ntipo==0) wprob2 = 'AL.CPI'
  CASE (14)
     wprob1 = 'Silicium'
     if(ntipo==0) wprob2 = 'SI.CPI'
  CASE (15)
     wprob1 = 'fosforo'
     if(ntipo==0) wprob2 = 'P.CPI'
  CASE (16)
     wprob1 = 'Azufre'
     if(ntipo==0) wprob2 = 'S.CPI'
  CASE (17)
     wprob1 = 'Cloro'
     if(ntipo==0) wprob2 = 'CL.CPI'
  CASE (18)
     wprob1 = 'Argon'
     if(ntipo==0) wprob2 = 'AR.CPI'
  CASE (19)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'K.CPI'
  CASE (20)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'CA.CPI'
  CASE (21)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'SC.CPI'
  CASE (22)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'TI.CPI'
  CASE (23)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'V.CPI'
  CASE (24)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'CR.CPI'
  CASE (25)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'MN.CPI'
  CASE (26)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'FE.CPI'
  CASE (27)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'CO.CPI'
  CASE (28)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'NI.CPI'
  CASE (29)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'CU.CPI'
  CASE (30)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'ZN.CPI'
  CASE (31)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'GA.CPI'
  CASE (32)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'GE.CPI'
  CASE (33)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'AS.CPI'
  CASE (34)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'SE.CPI'
  CASE (35)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'BR.CPI'
  CASE (36)
     wprob1 = 'Hydrogen'
     if(ntipo==0) wprob2 = 'ZR.CPI'
  END SELECT

end subroutine qua_especiename


subroutine leeespecie(kato)
!leeespecie(wprob2,atnumber,ncp,ncl,ncm,nlmax,nlocal,valencia,nrad,&
!                      radio,ppseu,ppphi,nspin)

use def_master
use def_quanty
implicit none
integer(ip), intent(in)  :: kato
!integer(ip), intent(in)  :: atnumber,nspin
!character(50), intent(in) :: wprob2

!integer(ip), intent(out)  :: ncp,ncl,ncm,nlmax,nlocal,valencia,nrad
!real(rp) ,allocatable     :: radio(:),ppseu(:,:),ppphi(:,:)
real(rp) :: A1,A2,A3,A4,amesh
integer(ip) :: KK,II,M
character(40) :: fname


  call qua_idespe(especie_qua(kato)%atnumber,especie_qua(kato)%nspin,especie_qua(kato)%ncp,especie_qua(kato)%ncl,especie_qua(kato)%ncm)

  fname="data\\"//adjustl(especie_qua(kato)%wprob2)
  open(unit=111,file= fname,status = 'unknown')

      
  if(especie_qua(kato)%ncp.eq.-1.or.especie_qua(kato)%ncl.eq.-1) then
     call outfor(4_ip,lun_outpu_qua,'SPECIES with bad number NCP and/or NCL')
  endif

  read(111,*) especie_qua(kato)%valencia, especie_qua(kato)%nlmax, especie_qua(kato)%nlocal
  
  especie_qua(kato)%nlmax=especie_qua(kato)%nlmax-1

  read(111,*) A1,A2,A3,A4 
      
  do kk=1,9
     read(111,*) A1,A2,A3 
  enddo
    
  read(111,*) especie_qua(kato)%nrad,amesh 
      
  allocate(especie_qua(kato)%radio(especie_qua(kato)%nrad), &
  especie_qua(kato)%ppseu(especie_qua(kato)%nlmax+1,especie_qua(kato)%nrad),&
  especie_qua(kato)%ppphi(especie_qua(kato)%nlmax+1,especie_qua(kato)%nrad))
      
  do ii=1,especie_qua(kato)%nlmax+1
        
     do kk=1,especie_qua(kato)%nrad
           read(111,*) M,especie_qua(kato)%radio(kk),especie_qua(kato)%ppphi(ii,kk),especie_qua(kato)%ppseu(ii,kk)  
           especie_qua(kato)%ppphi(ii,kk)=especie_qua(kato)%ppphi(ii,kk)/especie_qua(kato)%radio(kk)
   enddo
        
     if(ii.lt.especie_qua(kato)%nlmax+1) read(111,*) especie_qua(kato)%nrad,amesh 
   
  enddo

  close(111)  

end subroutine leeespecie

subroutine qua_idespe(atnumber,nspin,ncp,ncl,ncm)
  use def_master
  implicit none
  integer(ip), intent(in)  :: atnumber,nspin
  integer(ip), intent(out) :: ncp,ncl,ncm

  ! ojo para numeros atomicos grances hay lio en esta subrutina
  SELECT CASE (atnumber)
  CASE (1,2)
     NCP=1
     NCL=0
     NCM=0
  CASE (3,4)
     NCP=2
     NCL=0
     NCM=0
  CASE (5:10)
     NCP=2
     NCL=1
     SELECT CASE (atnumber)
     CASE (5,6)
        NCM=0
     CASE (7,8)
        NCM=-1
     CASE (9,10)
        NCM=1
     END SELECT
  CASE (11,12)
     NCP=3
     NCL=0
     NCM=0
  CASE (13:18)
     NCP=3
     NCL=1
     SELECT CASE (atnumber)
     CASE (13,14)
        NCM=0
     CASE (15,16)
        NCM=-1
     CASE (17,18)
        NCM=1
     END SELECT
  CASE (19:28)
     NCP=3
     NCL=2
     SELECT CASE (atnumber)
     CASE (19:20)
        NCM=0
     CASE (21:22)
        NCM=-1
     CASE (23:24)
        NCM=-2
     CASE (25:26)
        NCM=1
     CASE (27:28)
        NCM=2
     END SELECT
  CASE (29:30)
     NCP=4
     NCL=0
     NCM=0
  CASE (31:36)
     NCP=4
     NCL=1
     SELECT CASE (atnumber)
     CASE (31:32)
        NCM=0
     CASE (33:34)
        NCM=-1
     CASE (35:36)
        NCM=1
     END SELECT
  CASE (37:46)
     NCP=4
     NCL=2
     SELECT CASE (atnumber)
     CASE (37:38)
        NCM=0
     CASE (39:40)
        NCM=-1
     CASE (41:42)
        NCM=-2
     CASE (43:44)
        NCM=1
     CASE (45:46)
        NCM=2
     END SELECT
  CASE (47:48)
     NCP=5
     NCL=0
     NCM=0
  CASE (49:54)
     NCP=5
     NCL=1
     SELECT CASE (atnumber)
     CASE (49:50)
        NCM=0
     CASE (51:52)
        NCM=-1
     CASE (53:54)
        NCM=1
     END SELECT
  CASE (55:68)
     NCP=4
     NCL=3
     SELECT CASE (atnumber)
     CASE (55:56)
        NCM=0
     CASE (57:58)
        NCM=-1
     CASE (59:60)
        NCM=-2
     CASE (61:62)
        NCM=-3
     CASE (63:64)
        NCM=1
     CASE (65:66)
        NCM=2
     CASE (67:68)
        NCM=3
     END SELECT
  CASE (69:78)
     NCP=5
     NCL=2
     SELECT CASE (atnumber)
     CASE (69:70)
        NCM=0
     CASE (71:72)
        NCM=-1
     CASE (73:74)
        NCM=-2
     CASE (75:76)
        NCM=1
     CASE (77:78)
        NCM=2
     END SELECT
  CASE (79:80)
     NCP=6
     NCL=0
     NCM=0
  CASE (81:86)
     NCP=6
     NCL=1
     SELECT CASE (atnumber)
     CASE (81:82)
        NCM=0
     CASE (83:84)
        NCM=-1
     CASE (85:86)
        NCM=1
     END SELECT
  CASE (87:100)
     NCP=5
     NCL=4
     SELECT CASE (atnumber)
     CASE (87:88)
        NCM=0
     CASE (89:90)
        NCM=-1
     CASE (91:92)
        NCM=-2
     CASE (93:94)
        NCM=-3
     CASE (95:96)
        NCM=-4
     CASE (97:98)
        NCM=1
     CASE (99:100)
        NCM=2
        !             CASE (102:103)
        !           NCM=-3
        !             CASE (104:101)
        !           NCM=-4
     END SELECT
  CASE (101:110)
     NCP=6
     NCL=2
     SELECT CASE (atnumber)
     CASE (101:102)
        NCM=0
     CASE (103:104)
        NCM=-1
     CASE (105:106)
        NCM=-2
     CASE (107:108)
        NCM=1
     CASE (109:110)
        NCM=2
     END SELECT
  CASE (111:112)
     NCP=7
     NCL=0
     NCM=0
  CASE (113:118)
     NCP=7
     NCL=1
     SELECT CASE (atnumber)
     CASE (113:114)
        NCM=0
     CASE (115:116)
        NCM=-1
     CASE (117:118)
        NCM=1
     END SELECT
  CASE DEFAULT
     NCP = -1
     NCL = -1
  END SELECT

end subroutine qua_idespe


