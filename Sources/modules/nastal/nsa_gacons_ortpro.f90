subroutine nsa_gacons_ortpro(&
     ielem,pnode,elcod,elunk,elsub,elcon,eldif,elvel,elpre,eltem,elthe,elvis,eldtt)
  !----------------------------------------------------------------------------------
  !
  ! Gather operations 
  !
  !  d phi^a / dt  =  d F_{ia}(phi) / d x_i  =  d F_{ia} / d phi^b   d phi^b / dx_i
  !  
  !  where  a= 1,..,ndime+2  , b= 1,..,ndime+2  ,  i= 1,...,ndime
  !
  !  a,b label equations, i labels space
  !
  !  phi^a is elunk(a, pnode, ncomp_nsa)
  !
  !  d F_{ia} / d phi^b == elcon(a,b,i,pnode)
  !
  !
  !
  !-----------------------------------------------------------------------------------------------------------
  use      def_kintyp, only: ip,rp
  use      def_domain, only: mnode,ndime,lnods,coord
  use      def_master, only: veloc,densi,energ,tempe,press,umome
  use      def_nastal, only: ndofn_nsa,ncomp_nsa,kfl_track_nsa,denss_nsa,eness_nsa,umoss_nsa, &
       dtieq_nsa,cpcoe_nsa,prand_nsa,cvcoe_nsa,rgacv_nsa,kfl_taudi_nsa,ortpr_nsa
  implicit none
  integer(ip)  :: ielem,inode,idime,jdime,kdime,idofn,jdofn,pnode,ipoin
  real(rp)    :: &
       elunk(ndofn_nsa,mnode,ncomp_nsa),elsub(ndofn_nsa,mnode),eltun(ndofn_nsa,mnode), &
       elort(ndofn_nsa,mnode),elcon(ndofn_nsa,ndofn_nsa,ndime,mnode), &
       elcod(ndime,mnode),eldif(ndofn_nsa,ndofn_nsa,ndime,ndime,mnode), &
       eldtt(ndofn_nsa,mnode,2),elvel(ndime,mnode),velno(ndime),elpre(mnode),eltem(mnode), &
       velsq,enepe,visci,dicod, &
       dvite,elvis(mnode),elthe(mnode)

  do inode= 1,pnode
     ipoin= lnods(inode,ielem)
     elunk(ndime+1,inode,1) = densi(ipoin,1)
     elunk(ndime+1,inode,2) = densi(ipoin,2)
     elunk(ndime+1,inode,3) = densi(ipoin,ncomp_nsa)
     elunk(ndime+2,inode,1) = energ(ipoin,1)
     elunk(ndime+2,inode,2) = energ(ipoin,2)
     elunk(ndime+2,inode,3) = energ(ipoin,ncomp_nsa)
     elsub(ndime+1,inode  ) = denss_nsa(ipoin,1)
     elsub(ndime+2,inode  ) = eness_nsa(ipoin,1)
     eltun(ndime+1,inode  ) = elunk(ndime+1,inode,1)
     eltun(ndime+2,inode  ) = elunk(ndime+2,inode,1)
     if (kfl_track_nsa == 1) then
        eltun(ndime+1,inode  ) = elunk(ndime+1,inode,1) + elsub(ndime+1,inode)
        eltun(ndime+2,inode  ) = elunk(ndime+2,inode,1) + elsub(ndime+2,inode)
     end if
     eldtt(ndime+1,inode,1) = dtieq_nsa(2,ipoin,1)
     eldtt(ndime+1,inode,2) = dtieq_nsa(2,ipoin,2)
     eldtt(ndime+2,inode,1) = dtieq_nsa(3,ipoin,1)
     eldtt(ndime+2,inode,2) = dtieq_nsa(3,ipoin,2)

     call nsa_lawvis( -1 , 1 ,elvis(inode),tempe(ipoin,1),dvite)
     elthe(inode) = elvis(inode) * cpcoe_nsa / prand_nsa

     velsq= 0.0_rp
     elpre(inode)            = press(ipoin,1)
     eltem(inode)            = tempe(ipoin,1)
     do idime=1,ndime
        elcod(idime,inode  ) = coord(idime,ipoin  )          
        elunk(idime,inode,1) = umome(idime,ipoin,1)
        elunk(idime,inode,2) = umome(idime,ipoin,2)
        elunk(idime,inode,3) = umome(idime,ipoin,ncomp_nsa)
        elsub(idime,inode  ) = umoss_nsa(idime,ipoin,1)
        eltun(idime,inode  ) = elunk(idime,inode,1)
        if (kfl_track_nsa == 1) then
           eltun(idime,inode  ) = elunk(idime,inode,1) + elsub(idime,inode)
        end if
        elvel(idime,inode  ) = umome(idime,ipoin,1) / densi(ipoin,1)
        velno(idime)         = eltun(idime,inode) / eltun(ndime+1,inode)
        velsq                = velsq + velno(idime)*velno(idime)
        eldtt(idime,inode,1) = dtieq_nsa(1,ipoin,1)
        eldtt(idime,inode,2) = dtieq_nsa(1,ipoin,2)
        do idofn= 1,ndofn_nsa
           do jdofn= 1,ndofn_nsa
              elcon(idofn,jdofn,idime,inode) = 0.0_rp
              do jdime= 1,ndime
                 eldif(idofn,jdofn,idime,jdime,inode) = 0.0_rp
              end do
           end do
        end do
     end do
     enepe = eltun(ndime+2,inode) / eltun(ndime+1,inode) 
     visci = elvis(inode) / eltun(ndime+1,inode)
     dicod = elthe(inode) / cvcoe_nsa / eltun(ndime+1,inode)
     do kdime=1,ndime        
        elcon(ndime+1,kdime  ,kdime,inode)= 1.0_rp
        elcon(kdime  ,ndime+2,kdime,inode)= rgacv_nsa
        elcon(kdime  ,ndime+1,kdime,inode)= rgacv_nsa * 0.5_rp * velsq
        elcon(ndime+2,kdime  ,kdime,inode)= ((1.0_rp + rgacv_nsa) * enepe - rgacv_nsa * 0.5_rp * velsq)
        elcon(ndime+2,ndime+1,kdime,inode)= - velno(kdime) * ((1.0_rp + rgacv_nsa) * enepe - rgacv_nsa * velsq)
        elcon(ndime+2,ndime+2,kdime,inode)= ((1.0_rp + rgacv_nsa) * velno(kdime) )
        eldif(ndime+2,ndime+1,kdime,kdime,inode)= (dicod-visci) * velsq - dicod * enepe
        eldif(ndime+2,ndime+2,kdime,kdime,inode)= dicod
        do idime=1,ndime
           elcon(idime,idime  ,kdime,inode)= elcon(idime,idime  ,kdime,inode) + velno(kdime) 
           elcon(kdime,idime  ,kdime,inode)= elcon(kdime,idime  ,kdime,inode) - rgacv_nsa * velno(idime)
           elcon(idime,kdime  ,kdime,inode)= elcon(idime,kdime  ,kdime,inode) + velno(idime)
           elcon(idime,ndime+1,kdime,inode)= elcon(idime,ndime+1,kdime,inode) - velno(idime) * velno(kdime)
           elcon(ndime+2,idime,kdime,inode)= elcon(ndime+2,idime,kdime,inode) - &
                rgacv_nsa * velno(idime) * velno(kdime)
           eldif(kdime,kdime,idime,idime,inode)= visci
           eldif(kdime,idime,idime,kdime,inode)= eldif(kdime,idime,idime,kdime,inode) + visci
           eldif(kdime,idime,kdime,idime,inode)= eldif(kdime,idime,kdime,idime,inode) - 2.0_rp * visci / 3.0_rp
           eldif(kdime,ndime+1,idime,idime,inode)= - visci * velno(kdime)
           eldif(kdime,ndime+1,idime,kdime,inode)= eldif(kdime,ndime+1,idime,kdime,inode) - visci * velno(idime)
           eldif(kdime,ndime+1,kdime,idime,inode)= eldif(kdime,ndime+1,kdime,idime,inode) &
                + 2.0_rp * visci * velno(idime) / 3.0_rp
           eldif(ndime+2,idime,kdime,kdime,inode)= (visci-dicod) * velno(idime)
           eldif(ndime+2,kdime,kdime,idime,inode)= eldif(ndime+2,kdime,kdime,idime,inode) &
                + visci * velno(idime)
           eldif(ndime+2,kdime,idime,kdime,inode)= eldif(ndime+2,kdime,idime,kdime,inode) - &
                2.0_rp * visci * velno(idime) / 3.0_rp
           eldif(ndime+2,ndime+1,kdime,idime,inode)= eldif(ndime+2,ndime+1,kdime,idime,inode) &
                + 0.5_rp * visci * velno(kdime) * velno(idime)
        end do
     end do     
  end do

end subroutine nsa_gacons_ortpro
