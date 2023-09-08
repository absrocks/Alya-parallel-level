subroutine dod_qualit(ipoin,jpoin,kpoin,poext,Q,sign)
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master
  implicit none
  integer(ip),intent(in)     :: ipoin,jpoin,kpoin,poext
  real(rp),   intent(out)    :: Q
  integer(ip),intent(out)    :: sign
  integer(ip)                :: ipex1_tmp,idime,jpoin_tmp,pnode
  integer(ip)                :: lnods_tmp(mnode)
  real(rp)                   :: hmax, alpha,volum,sur1,sur2,sur3,sur4,surf,volum2
  real(rp)                   :: side1,side2,side3,side4,side5,side6,critv
  real(rp)                   :: cora1(ndime),cora2(ndime),cora4(ndime),increm
  real(rp)                   :: corb1(ndime),corb2(ndime),corb3(ndime)

  increm = epsilon(1.0_rp)
  alpha = sqrt(6.0_rp)/12.0_rp
  sign = 1_ip
  lnods_tmp(1)=ipoin
  lnods_tmp(2)=jpoin
  lnods_tmp(3)=kpoin
  lnods_tmp(4)=poext

  do idime=1,ndime
     cora1(idime) = coord(idime,lnods_tmp(1))-coord(idime,lnods_tmp(2))  
     corb1(idime) = coord(idime,lnods_tmp(4))-coord(idime,lnods_tmp(2)) 

     cora2(idime) = coord(idime,lnods_tmp(4))-coord(idime,lnods_tmp(1)) !side1
     corb2(idime) = coord(idime,lnods_tmp(3))-coord(idime,lnods_tmp(1)) !side2

!!!cora3(idime) = coord(idime,lnods_tmp(1))-coord(idime,lnods_tmp(2))===cora1  
     corb3(idime) = coord(idime,lnods_tmp(3))-coord(idime,lnods_tmp(2)) !side5

     cora4(idime) = coord(idime,lnods_tmp(4))-coord(idime,lnods_tmp(3)) !side6
!!!corb4(idime) = coord(idime,lnods_tmp(4))-coord(idime,lnods_tmp(2)) !side4
  end do


  !
  ! hmax
  !
  side1=0.0_rp
  side2=0.0_rp
  side3=0.0_rp
  side4=0.0_rp
  side5=0.0_rp
  side6=0.0_rp

  do idime=1,ndime 
     side1 = side1 + cora1(idime)*cora1(idime) !in base
     side2 = side2 + corb1(idime)*corb1(idime)
     side3 = side3 + cora2(idime)*cora2(idime)
     side4 = side4 + corb2(idime)*corb2(idime) !in base
     side5 = side5 + corb3(idime)*corb3(idime) !in base
     side6 = side6 + cora4(idime)*cora4(idime)
  end do


  side1 = sqrt(side1)
  side2 = sqrt(side2)
  side3 = sqrt(side3)
  side4 = sqrt(side4)
  side5 = sqrt(side5)
  side6 = sqrt(side6)
  hmax = max(side1,side2,side3,side4,side5,side6)
  !
  ! total surface
  !
   sur1 = ( cora1(1)*corb1(2) - cora1(2)*corb1(1) )*( cora1(1)*corb1(2) - cora1(2)*corb1(1) ) +&
          ( cora1(3)*corb1(1) - cora1(1)*corb1(3) )*( cora1(3)*corb1(1) - cora1(1)*corb1(3) ) +&
          ( cora1(2)*corb1(3) - cora1(3)*corb1(2) )*( cora1(2)*corb1(3) - cora1(3)*corb1(2) )
   sur1=sqrt(sur1)*0.5_rp 


   sur2 = ( cora2(1)*corb2(2) - corb2(1)*cora2(2) )*( cora2(1)*corb2(2) - corb2(1)*cora2(2) ) +&
          ( cora2(3)*corb2(1) - cora2(1)*corb2(3) )*( cora2(3)*corb2(1) - cora2(1)*corb2(3) ) +&
          ( cora2(2)*corb2(3) - cora2(3)*corb2(2) )*( cora2(2)*corb2(3) - cora2(3)*corb2(2) )
   sur2=sqrt(sur2)*0.5_rp


   sur3 = ( cora1(1)*corb3(2) - corb3(1)*cora1(2) )*( cora1(1)*corb3(2) - corb3(1)*cora1(2) ) +&
          ( cora1(3)*corb3(1) - cora1(1)*corb3(3) )*( cora1(3)*corb3(1) - cora1(1)*corb3(3) ) +&
          ( cora1(2)*corb3(3) - cora1(3)*corb3(2) )*( cora1(2)*corb3(3) - cora1(3)*corb3(2) )
   sur3=sqrt(sur3)*0.5_rp

   sur4 = ( cora4(1)*corb1(2) - corb1(1)*cora4(2) )*( cora4(1)*corb1(2) - corb1(1)*cora4(2) ) +&
          ( cora4(3)*corb1(1) - cora4(1)*corb1(3) )*( cora4(3)*corb1(1) - cora4(1)*corb1(3) ) +&
          ( cora4(2)*corb1(3) - cora4(3)*corb1(2) )*( cora4(2)*corb1(3) - cora4(3)*corb1(2) )
   sur4=sqrt(sur4)*0.5_rp
   surf = sur1 + sur2 + sur3 + sur4

   !
   ! volume
   !
   

   volum = (  cora1(1)*cora2(2)*corb2(3) + cora1(2)*cora2(3)*corb2(1) + cora1(3)*cora2(1)*corb2(2) &
        &   - cora1(3)*cora2(2)*corb2(1) - cora1(2)*cora2(1)*corb2(3) - cora1(1)*cora2(3)*corb2(2) ) /6.0_rp
 
   if (volum <= 0.0_rp) sign=-1_ip

   !
   ! Quality
   !

     Q = alpha * hmax * surf /(3.0_rp*abs(volum)+increm )


end subroutine dod_qualit
