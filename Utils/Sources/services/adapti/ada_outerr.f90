subroutine ada_outerr
!-----------------------------------------------------------------------
!****f* adapti/ada_outerr
! NAME 
!    ada_outerr
! DESCRIPTION
!    This routine writes the chosen error estimator in a gid post.res file
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_inpout
  use      def_adapti
  implicit none
  integer(ip)              :: &
       lutmp,ielem,iesva,jesva,idime,inode,pelty,pnode,ipoin,iunkn

  if (kfl_redgr_ada < 10) return
  
  if (ndime == 2) then

     do ipoin= 1,npoin
        do iesva= 1,ndime+2
           iunkn= (ipoin-1)*(ndime+2)+iesva 
           unkno(iunkn)= 0.0_rp
        end do
     end do

     do ielem= 1,nelem
        pelty=ltype(ielem)
        pnode=nnode(pelty)
        do iesva= 1,ndime+2
           jesva= (ielem-1)*(ndime+2)+iesva
           do inode= 1,pnode
              ipoin= lnods(inode,ielem)
              iunkn= (ipoin-1)*(ndime+2)+iesva 
              unkno(iunkn)= unkno(iunkn) + erres(jesva)
           end do
        end do
     end do

     write(lun_ouerr_ada,'(a)') 'GiD Post Results File 1.0'
     write(lun_ouerr_ada,'(a)') 

     write(lun_ouerr_ada,'(a)') 'Result ERR_CONTI ANALYSIS 1 Scalar OnNodes '
     write(lun_ouerr_ada,'(a)') 'Values'
     do ipoin= 1,npoin
        iunkn= (ipoin-1)*(ndime+2)+ndime+1 
        write(lun_ouerr_ada,*) ipoin, unkno(iunkn)      
     end do
     write(lun_ouerr_ada,'(a)') 'End Values'
     
     write(lun_ouerr_ada,'(a)') 'Result ERR_ENERG ANALYSIS 1 Scalar OnNodes '
     write(lun_ouerr_ada,'(a)') 'Values'
     do ipoin= 1,npoin
        iunkn= (ipoin-1)*(ndime+2)+ndime+2 
        write(lun_ouerr_ada,*) ipoin, unkno(iunkn)      
     end do
     write(lun_ouerr_ada,'(a)') 'End Values'
     
     write(lun_ouerr_ada,'(a)') 'Result ERR_MOMEN ANALYSIS 1 Vector OnNodes'
     write(lun_ouerr_ada,'(a)') 'ComponentNames MOMEN_X , MOMEN_Y , MOMEN_Z'    
     write(lun_ouerr_ada,'(a)') 'Values'
     do ipoin= 1,npoin
        iunkn= (ipoin-1)*(ndime+2) 
        write(lun_ouerr_ada,*) ipoin, (unkno(iunkn+idime),idime=1,ndime)      
     end do
     write(lun_ouerr_ada,'(a)') 'End Values'

  end if

end subroutine ada_outerr
