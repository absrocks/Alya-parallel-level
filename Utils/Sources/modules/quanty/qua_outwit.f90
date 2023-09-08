subroutine qua_outwit()
  !------------------------------------------------------------------------
  !****f* Quanty/qua_outwit
  ! NAME 
  !    qua_outwit
  ! DESCRIPTION
  !    Output values of phi on witness points
  ! USES
  ! USED BY
  !    qua_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_quanty
  implicit none
  integer(ip), save :: ipass=0
  integer(ip)       :: iwitn,ielem,inode,pnode,pelty,ipoin
  real(rp)          :: wtemp

  if(kfl_paral==-1) then
     if(nwitn>0.and.maxval(npp_witne_qua)>0) then
        if(kfl_paral<=0.and.ipass==0) then
           ipass=1
           write(lun_witne_qua,1)
        end if
        if(kfl_paral<=0) write(lun_witne_qua,2) cutim
        !
        ! Results on witness points
        !
        if(kfl_paral/=0) then 
           do iwitn=1,nwitn
              ielem=lewit(iwitn)
              if(ielem>0) then
                 pelty=ltype(ielem)
                 pnode=nnode(pelty)
                 if(npp_witne_qua(1)==1) then
                    wtemp=0.0_rp
                    do inode=1,pnode
                       ipoin=lnods(inode,ielem)
                       wtemp=wtemp+shwit(inode,iwitn)*phion(ipoin,1)
                    end do
                    write(lun_witne_qua,2) wtemp
                 end if
              end if
           end do
        end if
        !
        ! Exchange master-slaves
        !
        if(kfl_paral>=0) then

        end if
     end if
  end if

1 format('$ Time desnsity ->',100(1x,e12.6))
2 format(1x,e12.6,$)

end subroutine qua_outwit

