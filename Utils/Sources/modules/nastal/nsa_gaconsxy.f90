subroutine nsa_gaconsxy(&
     ielem,pnode,elcod,elunk,elsub,elcon,eldif,elvel,elmsh,elpre,eltem,&
     elvis,elhcp,elwme,eldtt,elort,elsax)
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
  use      def_master
  use      def_kermod
  use      mod_ker_proper
  use      def_nastal, only: ndofn_nsa,ncomp_nsa,kfl_track_nsa,denss_nsa,eness_nsa,umoss_nsa, &
                             dtieq_nsa,cpcoe_nsa,prand_nsa,kfl_taudi_nsa,ortpr_nsa,mowei_nsa,runiv_nsa,&
                             nevat_nsa,pabdf_nsa,kfl_tisch_nsa

  implicit none
  integer(ip),  intent(in) :: ielem,pnode 
  real(rp),    intent(out) :: elcod(ndime,mnode),elunk(ndofn_nsa,mnode,ncomp_nsa),elsub(ndofn_nsa,mnode), &
                              elort(ndofn_nsa,mnode),elcon(ndofn_nsa,ndofn_nsa,ndime,mnode),              &
                              eldif(ndofn_nsa,ndofn_nsa,ndime,ndime,mnode),eldtt(ndofn_nsa,mnode,2),      &
                              elvel(ndime,mnode),elmsh(ndime,mnode),elpre(mnode),eltem(mnode),elvis(mnode),&
                              elhcp(pnode),elwme(pnode),elsax(nevat_nsa)
  integer(ip)              :: inode,idime,jdime,kdime,idofn,jdofn,ipoin,ievat,itott,dummi,itime_scheme
  real(rp)                 :: dummy(ndime,ndime),auxvi(ndofn_nsa),velno(ndime),elthe(mnode)
  real(rp)                 :: prope_tmp(pnode),velsq,enepe,visci,dicod,dvite,eltun(ndofn_nsa,mnode), &
                              elhcv(pnode),rgacv
 
  ! Initialization of dummy variables for viscosity
  dummy = 0.0_rp
  auxvi = 0.0_rp
  elmsh = 0.0_rp
  !
  ! Properties: viscosity mu, c_p 
  !
  if (kfl_prope /= 0 ) then
     !call ker_proper('VISCO','PNODE',dummi,ielem,prope_tmp,pnode,dummi,dummy,dummy)
     call ker_proper('VISCO','PNODE',dummi,ielem,prope_tmp,pnode)
     elvis(1:pnode) = prope_tmp(1:pnode)
     !call ker_proper('SPHEA','PNODE',dummi,ielem,prope_tmp,pnode,dummi,dummy,dummy)
     call ker_proper('SPHEA','PNODE',dummi,ielem,prope_tmp,pnode)
     elhcp(1:pnode) = prope_tmp(1:pnode)
     elwme  = mowei_nsa
  else
     elhcp  = cpcoe_nsa
     elwme  = mowei_nsa
  endif 

  itime_scheme= ITER_K  ! this is required to avoid forbidden memory acces of the global vectors
  if (kfl_tisch_nsa == 2) itime_scheme= TIME_N_MINUS_1

  do inode= 1,pnode
     !
     ! OJOOOO: SO FAR ONLY BDF 1 AND 2 ARE PROGRAMMED
     !         
     ! There is a "-" because of how the bdf parameters are defined
     !         
     ipoin= lnods(inode,ielem)

!     elunk(ndime+1,inode,ITER_K) =  (- pabdf_nsa(2) * densi(ipoin,ITER_K) - pabdf_nsa(3) * densi(ipoin,itime_scheme))/pabdf_nsa(1) 
!     elunk(ndime+2,inode,ITER_K) =  (- pabdf_nsa(2) * energ(ipoin,ITER_K) - pabdf_nsa(3) * energ(ipoin,itime_scheme))/pabdf_nsa(1)
!     elunk(ndime+1,inode,TIME_N) =  (- pabdf_nsa(2) * densi(ipoin,TIME_N) - pabdf_nsa(3) * densi(ipoin,itime_scheme))/pabdf_nsa(1)
!     elunk(ndime+2,inode,TIME_N) =  (- pabdf_nsa(2) * energ(ipoin,TIME_N) - pabdf_nsa(3) *  energ(ipoin,itime_scheme))/pabdf_nsa(1)
!     elunk(ndime+1,inode,ITER_AUX) =  (- pabdf_nsa(2) * densi(ipoin,ITER_AUX) - pabdf_nsa(3) * densi(ipoin,itime_scheme))/pabdf_nsa(1)
!     elunk(ndime+2,inode,ITER_AUX) =  (- pabdf_nsa(2) * energ(ipoin,ITER_AUX) - pabdf_nsa(3) * energ(ipoin,itime_scheme))/pabdf_nsa(1)
!     do idime=1,ndime
!        elunk(idime,inode,ITER_K)   =  &
!             (- pabdf_nsa(2) * umome(idime,ipoin,ITER_K) - pabdf_nsa(3) * umome(idime,ipoin,itime_scheme))/pabdf_nsa(1)
!        elunk(idime,inode,ITER_AUX) =  &
!             (- pabdf_nsa(2) * umome(idime,ipoin,ITER_AUX) - pabdf_nsa(3) * umome(idime,ipoin,itime_scheme))/pabdf_nsa(1)
!        elunk(idime,inode,TIME_N)   =  &
!             (- pabdf_nsa(2) * umome(idime,ipoin,TIME_N) - pabdf_nsa(3) * umome(idime,ipoin,itime_scheme))/pabdf_nsa(1)
!     end do

     elunk(ndime+1,inode,ITER_K) = densi(ipoin,ITER_K)
     elunk(ndime+2,inode,ITER_K) = energ(ipoin,ITER_K)
     elunk(ndime+1,inode,TIME_N) = densi(ipoin,TIME_N)
     elunk(ndime+2,inode,TIME_N) = energ(ipoin,TIME_N)
     elunk(ndime+1,inode,ITER_AUX) = densi(ipoin,ITER_AUX)
     elunk(ndime+2,inode,ITER_AUX) = energ(ipoin,ITER_AUX)
     do idime=1,ndime
        elunk(idime,inode,ITER_K)   =  &
            umome(idime,ipoin,ITER_K)
        elunk(idime,inode,ITER_AUX) =  &
            umome(idime,ipoin,ITER_AUX)
        elunk(idime,inode,TIME_N)   =  &
            umome(idime,ipoin,TIME_N)
     end do

     elsub(ndime+1,inode  ) = denss_nsa(ipoin,1)
     elsub(ndime+2,inode  ) = eness_nsa(ipoin,1)
     eltun(ndime+1,inode  ) = elunk(ndime+1,inode,1)
     eltun(ndime+2,inode  ) = elunk(ndime+2,inode,1)
     elort(ndime+1,inode  ) = ortpr_nsa(ndime+1,ipoin)
     elort(ndime+2,inode  ) = ortpr_nsa(ndime+2,ipoin)

     do idofn = 1,ndofn_nsa
        ievat = (inode-1) * ndofn_nsa + idofn
        itott = (ipoin-1) * ndofn_nsa + idofn
!!        elsax(ievat) = rhsou_nsa(itott)
     end do

     if (kfl_coupl(ID_NASTAL,ID_CHEMIC) >= 1 ) then       ! Molecular weight of the mixture
       elwme(inode) = wmean(ipoin,1)
     endif

     if (kfl_track_nsa == 1) then
        eltun(ndime+1,inode) = elunk(ndime+1,inode,ITER_K) + elsub(ndime+1,inode)
        eltun(ndime+2,inode) = elunk(ndime+2,inode,ITER_K) + elsub(ndime+2,inode)
     end if
     eldtt(ndime+1,inode,1) = dtieq_nsa(2,ipoin,1)
     eldtt(ndime+1,inode,2) = dtieq_nsa(2,ipoin,2)
     eldtt(ndime+2,inode,1) = dtieq_nsa(3,ipoin,1)
     eldtt(ndime+2,inode,2) = dtieq_nsa(3,ipoin,2)

     !
     ! Properties: viscosity mu, c_p 
     !
     if (kfl_prope /= 0 ) then
        elthe(inode) = elvis(inode) * elhcp(inode) / prand_nsa
     else 
        call nsa_lawvis( -1 , 1 ,elvis(inode),tempe(ipoin,1),dvite) 
        elthe(inode) = elvis(inode) * elhcp(inode) / prand_nsa
     endif

     velsq = 0.0_rp
     elpre(inode)            = press(ipoin,1)
     eltem(inode)            = tempe(ipoin,1)


     do idime= 1,ndime

        elcod(idime,inode  ) = coord(idime,ipoin  )          

        elsub(idime,inode  ) = umoss_nsa(idime,ipoin,1)
        eltun(idime,inode  ) = elunk(idime,inode,ITER_K) 
        elort(idime,inode  ) = ortpr_nsa(idime,ipoin)
        if (kfl_track_nsa == 1) then
           eltun(idime,inode) = elunk(idime,inode,ITER_K) + elsub(idime,inode)
        end if
        elvel(idime,inode  ) = elunk(idime,inode,ITER_K) / elunk(ndime+1,inode,ITER_K)
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

     elhcv(inode) = elhcp(inode) - runiv_nsa / elwme(inode)
     rgacv = runiv_nsa / elwme(inode) / elhcv(inode)

     enepe = eltun(ndime+2,inode) / eltun(ndime+1,inode) 
     visci = elvis(inode) / eltun(ndime+1,inode)
     dicod = elthe(inode) / elhcv(inode) / eltun(ndime+1,inode)
     do kdime=1,ndime        
        elcon(ndime+1,kdime  ,kdime,inode)= 1.0_rp
        elcon(kdime  ,ndime+2,kdime,inode)= rgacv
        elcon(kdime  ,ndime+1,kdime,inode)= rgacv * 0.5_rp * velsq
        elcon(ndime+2,kdime  ,kdime,inode)= ((1.0_rp + rgacv) * enepe - rgacv * 0.5_rp * velsq)
        elcon(ndime+2,ndime+1,kdime,inode)= - velno(kdime) * ((1.0_rp + rgacv) * enepe - rgacv * velsq)
        elcon(ndime+2,ndime+2,kdime,inode)= ((1.0_rp + rgacv) * velno(kdime) )
        eldif(ndime+2,ndime+1,kdime,kdime,inode)= (dicod-visci) * velsq - dicod * enepe
        eldif(ndime+2,ndime+2,kdime,kdime,inode)= dicod
        do idime=1,ndime
           elcon(idime,idime  ,kdime,inode)= elcon(idime,idime  ,kdime,inode) + velno(kdime) 
           elcon(kdime,idime  ,kdime,inode)= elcon(kdime,idime  ,kdime,inode) - rgacv * velno(idime)
           elcon(idime,kdime  ,kdime,inode)= elcon(idime,kdime  ,kdime,inode) + velno(idime)
           elcon(idime,ndime+1,kdime,inode)= elcon(idime,ndime+1,kdime,inode) - velno(idime) * velno(kdime)
           elcon(ndime+2,idime,kdime,inode)= elcon(ndime+2,idime,kdime,inode) - &
                rgacv * velno(idime) * velno(kdime)
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

  !
  ! Mesh velocity
  !     
  if( kfl_coupl(ID_NASTAL,ID_ALEFOR) /= 0 ) then  
     do inode = 1,pnode
        ipoin = lnods(inode,ielem)
        do idime = 1,ndime
           elmsh(idime,inode) = velom(idime,ipoin)
        end do
     end do
  end if



end subroutine nsa_gaconsxy
