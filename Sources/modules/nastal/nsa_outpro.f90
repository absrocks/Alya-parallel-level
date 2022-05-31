subroutine nsa_outpro(itask)
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_outpro
  ! NAME 
  !    nsa_outpro
  ! DESCRIPTION
  !    This routine outputs 2d profile distributions
  ! USES
  ! USED BY
  !    nsa_posbri
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_nastal
  use      mod_iofile
  use      mod_memchk
  implicit none
  integer(ip),intent(in)   :: itask
  integer(ip)              :: ipoin,kpoin,ifibo,idofn,ipopr,idime,jfibo
  real(rp)                 :: velmo,xvelo,xcpre

  if (kfl_paral==0.and.kfl_postp_par==0) return
 
  if (itask==0 .or. itask==1) then
     ipopr=0
     do ipoin = 1,npoin
        linod_nsa(ipoin) = 0
        ifibo=0
        do idime= 1,ndime
           jfibo= kfl_fixno_nsa(idime,ipoin)
           if (jfibo>0) ifibo= ifibo+jfibo
        end do
        if (ifibo>0) then
           velmo=0.0_rp
           do idime=1,ndime
              xvelo= veloc(idime,ipoin,3)
              velmo= velmo + xvelo*xvelo
           end do
           if((kfl_fixno_nsa(1,ipoin)==2) .or. (velmo<zensa)) then
              ipopr=ipopr+1
              if (itask == 0) then
                 npopr_nsa = npopr_nsa + 1
              else 
                 linod_nsa(ipopr) = ipoin
              end if
           end if
        end if
     end do
     
  else if (itask==2) then
     
     write(lun_pro2d_nsa,*) '# COORDINATES '
     write(lun_pro2d_nsa,*) '# Points      = ',npopr_nsa
     do ipopr=1,npopr_nsa
        ipoin=linod_nsa(ipopr)
        write(lun_pro2d_nsa,100) ipopr,ipoin,(coord(idime,ipoin),idime=1,ndime)
     end do

     kfl_pro2d_nsa = 2
              
  else if (itask==3) then
     
     write(lun_pro2d_nsa,*) '# PRESSURE_COEFFICIENT '
     write(lun_pro2d_nsa,*) '# Iterations  = ',ittim,itcou,itinn(modul)
     write(lun_pro2d_nsa,*) '# Cutim       = ',cutim
     write(lun_pro2d_nsa,*) '# Points      = ',npopr_nsa
     do ipopr=1,npopr_nsa
        ipoin= linod_nsa(ipopr)
        xcpre= 2.0_rp*(press(ipoin,1)-press_nsa)/speed_nsa/speed_nsa
        unkno(ipopr)=xcpre
        write(lun_pro2d_nsa,100) ipopr,ipoin,unkno(ipopr)
     end do
              
  end if

  
100 format(2x,i9,2x,i9,3(2x,e12.6))

  
end subroutine nsa_outpro
