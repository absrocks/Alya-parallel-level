subroutine ibm_boboib(iimbo)
  !-----------------------------------------------------------------------
  !****f* ibm_boboib/ibm_boboib
  ! NAME
  !    ibm_boboib
  ! DESCRIPTION
  !    Determine the bounding box for each particle and for each 
  !    face in a particle
  ! USED BY
  !    ibm_iniunk
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,nnode
  use def_master, only       :  imbou
  implicit none
  integer(ip),    intent(in) :: iimbo
  integer(ip)                :: idime,iboib,ipoib,inoib    
  real(rp)                   :: dista,toler
  !
  ! Loop over faces of a particle
  !
  toler = 1.0e-3_rp

  do iboib = 1,imbou(iimbo) % nboib       
    
     do idime = 1,ndime

        imbou(iimbo) % fabox(1,idime,iboib) =  1.0e12_rp
        imbou(iimbo) % fabox(2,idime,iboib) = -1.0e12_rp

        do inoib = 1,nnode(imbou(iimbo) % ltyib(iboib))
           ipoib = imbou(iimbo) % lnoib(inoib,iboib)
           imbou(iimbo) % fabox(1,idime,iboib) = min( imbou(iimbo) % fabox(1,idime,iboib) , imbou(iimbo) % cooib(idime,ipoib) )
           imbou(iimbo) % fabox(2,idime,iboib) = max( imbou(iimbo) % fabox(2,idime,iboib) , imbou(iimbo) % cooib(idime,ipoib) )
        end do

        dista = imbou(iimbo) % fabox(2,idime,iboib) - imbou(iimbo) % fabox(1,idime,iboib)
        imbou(iimbo) % fabox(1,idime,iboib) = imbou(iimbo) % fabox(1,idime,iboib) - abs(dista*toler)
        imbou(iimbo) % fabox(2,idime,iboib) = imbou(iimbo) % fabox(2,idime,iboib) + abs(dista*toler)

     end do

  end do

  do idime = 1,ndime
     imbou(iimbo) % bobox(idime,1) =  1.0e12_rp
     imbou(iimbo) % bobox(idime,2) = -1.0e12_rp
  end do

  do iboib = 1,imbou(iimbo) % nboib
     do idime = 1,ndime
        imbou(iimbo) % bobox(idime,1) = min( imbou(iimbo) % bobox(idime,1) , imbou(iimbo) % fabox(1,idime,iboib) )
        imbou(iimbo) % bobox(idime,2) = max( imbou(iimbo) % bobox(idime,2) , imbou(iimbo) % fabox(2,idime,iboib) )
     end do
  end do

end subroutine ibm_boboib
