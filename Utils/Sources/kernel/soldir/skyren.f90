subroutine skyren(&
     nsist,ndofn,npoin,ntotv,veori,veopt,&
     lpntn,itask)
  !-----------------------------------------------------------------------
  !
  ! This routine redefines vectors according to the renumbering strategy
  !
  !-----------------------------------------------------------------------
  use def_kintyp
  implicit none
  integer(ip) :: nsist,ndofn,npoin,ntotv,itask
  integer(ip) :: lpntn(npoin)
  real(rp)    :: veori(ntotv,nsist), veopt(ntotv,nsist)
  integer(ip) :: isist,ipoin,itoto,itotn,idofn

  if(itask==1) then
     !
     ! Original --> optimized
     !
     do isist = 1,nsist
        do ipoin = 1,npoin
           itoto = (ipoin-1)*ndofn
           itotn = (lpntn(ipoin)-1)*ndofn
           do idofn = 1,ndofn
              itoto = itoto + 1
              itotn = itotn + 1
              veopt(itotn,isist) = veori(itoto,isist)
           end do
        end do
     end do

  else if(itask==2) then
     !
     ! Optimized --> original
     !
     do isist = 1,nsist
        do ipoin = 1,npoin
           itoto = (ipoin-1)*ndofn
           itotn = (lpntn(ipoin)-1)*ndofn
           do idofn = 1,ndofn
              itoto = itoto + 1
              itotn = itotn + 1
              veori(itoto,isist) = veopt(itotn,isist) 
           end do
        end do
     end do
  end if

end subroutine skyren
