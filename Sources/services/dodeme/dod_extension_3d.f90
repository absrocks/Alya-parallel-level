subroutine dod_extension_3d
  use def_kintyp
  use def_parame
  use def_master
  use def_elmtyp
  use def_domain
  use def_dodeme
  use def_kermod
  use mod_memory
  use mod_kdtree
  use mod_dod_extens
  use mod_elmgeo
  use mod_messages, only : livinf
  implicit none
  integer(ip)              :: isubd,jsubd,ielem,jelem,ielem_old,ielem_aux,kelem_aux
  integer(ip)              :: ipoin,ipoin_subd
  integer(ip)              :: inext,jnext,mnext,nnext,npoin_subd
  integer(ip)              :: idime,istat,inode,jnode,ibopo,inodb
  integer(ip)              :: eboun,iboun,jboun,kboun
  integer(ip)              :: mboun,nboun_local,qboun,rboun,sboun
  integer(ip)              :: ii,jj,kk,nedge_ipoin
  integer(ip)              :: jview,iview,kview
  integer(ip)              :: ipoex,jpoex,ipext,npoex,tetex,nutri
  integer(ip)              :: poext,poex1,poex2,poex3,poex4
  integer(ip)              :: poexk,poexq,poexr,poexs
  integer(ip)              :: pobok,poboq,pobos,pobor
  integer(ip)              :: lperm(3,4),lpern(3,4),lpert(3,4)
  integer(ip)              :: cont,flag,tmp,sign
  integer(ip)              :: imatn,ipoiz,kelem,knode,dumm1(2)
  real(rp)                 :: orien(3),norma ,dumm5(3) !,normal(3,3)
  real(rp)                 :: time1,xieta(3)
  real(rp)                 :: q,q1,q2,q_old,kappa,asrad,qaux,aux
  real(rp)                 :: hleng(2),coorcg2(3),hcandi,coorcg(3)
  integer(ip), pointer     :: para(:),list_ordered_boundaries(:),lextp(:,:),lcont(:,:,:)
  integer(ip), pointer     :: pasa(:),indic(:)
  integer(ip), pointer     :: lext2(:) ,lnext(:)
  integer(ip), pointer     :: lnods_aux(:,:),ltype_aux(:),lmate_aux(:),lelez_aux(:),lesub_aux(:)
  integer(ip), allocatable :: lpoex(:,:) !!!Guillaume->si lo pongo pointer, 
!!!no puedo pasarlo como argumento 
!!!a dod_extension_elements_3d
  !!  real(rp),pointer     :: codel(:,:),bocod(:,:)
  real(rp), pointer        :: normx(:),normy(:),normz(:),hcarac(:),coorcg1(:,:)
  real(rp), pointer        :: quali(:,:),lqaux(:),q_max(:),distex(:,:)
  logical(lg)              :: caplim


  call cshder(3_ip) !!!OJOOO-Guillaume
  call livinf(0_ip,'CREATE 3D EXTENSION ELEMENTS',0_ip)

  nelem_dod = 0 
  nullify(para)
  nullify(list_ordered_boundaries)
  nullify(lextp)
  nullify(lcont)
  nullify(pasa)
  nullify(indic)
  nullify(lext2)
  nullify(lnext)
  nullify(lnods_aux)
  nullify(ltype_aux)
  nullify(lmate_aux)
  nullify(lelez_aux)
  nullify(lesub_aux)
  nullify(normx)
  nullify(normy)
  nullify(normz)
  nullify(hcarac)
  nullify(coorcg1)
  nullify(quali)
  nullify(lqaux)
  nullify(q_max)
  nullify(distex)


  !------------------------------------------------------
  !
  ! MNEXT: maximum number of extensions elements
  !
  !------------------------------------------------------

  mnext = 0
  number_fringe_nodes = 0
  do isubd = 1,nsubd
     current_subdomain => subdomain(isubd)
     do ipoin = 1, current_subdomain % npoin
        if( current_subdomain % lsubd_npoin(ipoin) > 0 ) then
           number_fringe_nodes = number_fringe_nodes + 1
           mnext = max(mnext,current_subdomain % extension(ipoin) % number_candidates)
        end if
     end do
  end do

  call memory_alloca(mem_servi(1:2,servi),'LNEXT',    'dod_extension_3d',lnext,mnext)
  !------------------------------------------------------
  !
  ! Allocate memory for fringe nodes type
  !
  !------------------------------------------------------

  call dod_memall(6_ip)

!!$
!!$  !------------------------------------------------------
!!$  !
!!$  ! Para cada punto de contorno:
!!$  ! KFL_INTOP_DOD = 0 in general case
!!$  !               = 1 if an extension node is on the real boundary
!!$  ! +-----oo---+---
!!$  !    i  || j
!!$  !       ||
!!$  !       ++
!!$  !
!!$  !------------------------------------------------------
!!$
  call cputim(time1)
  ielem = 0
  number_fringe_nodes = 0
  orien = 0.0_rp
!!$isubd=2
!!$    current_subdomain => subdomain(isubd)
!!$     npoin_subd        =  current_subdomain % npoin
!!$     mboun = 0
!!$     do ipoin_subd = 1,npoin_subd 
!!$        ipoin    =  current_subdomain % lnper(ipoin_subd)
!!$print*,ipoin
!!$     do ibopo = pbopo_global(ipoin),pbopo_global(ipoin+1)-1
!!$        iboun = lbopo_global(ibopo)
!!$        if(lboch_global(iboun)/=BOHOL)then
!!$           print*,'LOS CONTORNIS', nboun_local,lnodb_global(:,iboun)
!!$        end if
!!$     end do
!!$     end do
  do isubd = 1,nsubd

     current_subdomain => subdomain(isubd)
     npoin_subd        =  current_subdomain % npoin
     mboun = 0
     do ipoin_subd = 1,npoin_subd 
        ipoin    =  current_subdomain % lnper(ipoin_subd)
        mboun     = max(mboun,pbopo_global(ipoin+1)-pbopo_global(ipoin))
     end do
     call memory_alloca(mem_servi(1:2,servi),'DISTEX',   'dod_extension_3d',distex,mboun,mnext)  !!!GUILLAUME:...ESto solo necesito en caso de capa limite con hexass.....como lo hago??
     call memory_alloca(mem_servi(1:2,servi),'HCARAC','dod_extension_3d',hcarac,mboun)
     call memory_alloca(mem_servi(1:2,servi),'COORCG1','dod_extension_3d',coorcg1,3_ip,mboun)
     call memory_alloca(mem_servi(1:2,servi),'LIST_ORDERED_BOUNDARIES','dod_extension_3d',list_ordered_boundaries,mboun)
     call memory_alloca(mem_servi(1:2,servi),'LEXT2','dod_extension_3d',lext2,mboun)
     call memory_alloca(mem_servi(1:2,servi),'PARA','dod_extension_3d',para,mboun)
     call memory_alloca(mem_servi(1:2,servi),'NORMX','dod_extension_3d',normx,mboun)
     call memory_alloca(mem_servi(1:2,servi),'NORMY','dod_extension_3d',normy,mboun)
     call memory_alloca(mem_servi(1:2,servi),'NORMZ','dod_extension_3d',normz,mboun)

     punto_contorno: do ipoin_subd = 1,npoin_subd

        if( current_subdomain % lsubd_npoin(ipoin_subd) > 0 ) then
           caplim = .false.
           number_fringe_nodes =  number_fringe_nodes + 1
           ipoin               =  current_subdomain % lnper(ipoin_subd)
          ! print*,'--------------------ipoin--------------------------------------------------',ipoin,isubd
           jsubd               =  current_subdomain % lsubd_npoin(ipoin_subd)
           neighbor_subdomain  => subdomain(jsubd)
           nnext               =  current_subdomain  % extension(ipoin_subd) % number_candidates 
           do inext=1,nnext
              lnext  (inext)  = neighbor_subdomain % lnper(current_subdomain % extension(ipoin_subd) % candidates(inext))
           end do
           nboun_local         = pbopo_global(ipoin+1)-pbopo_global(ipoin)!!OJOOOO
           nboun_local = 0
           do ibopo = pbopo_global(ipoin),pbopo_global(ipoin+1)-1
              iboun = lbopo_global(ibopo)
              if(lboch_global(iboun)/=BOHOL)then
                 nboun_local = nboun_local + 1  
                 if(ipoin==25748)print*,'LOS CONTORNIS', nboun_local,lnodb_global(:,iboun)
              end if
           end do
           do ibopo = pbopo_global(ipoin),pbopo_global(ipoin+1)-1
              mnodb = max(mnodb,nnode(abs(ltypb_global(lbopo_global(ibopo)))))
           end do
           do jboun = 1,nboun_local
              para(jboun)  = 0_ip
              normx(jboun) = 0.0_rp
              normy(jboun) = 0.0_rp
              normz(jboun) = 0.0_rp
           end do

           !------------------------------------------------------
           ! 
           ! Sorted list of boundaries associated to ipoin 
           !
           !------------------------------------------------------
           nedge_ipoin = size(ledge_npoin_global(ipoin) % l)
if(ipoin==25748)print*,'entro a dod_order_boundaries',nboun_local
           call dod_order_boundaries(&
                1_ip,ipoin,mnodb,nboun_local,npoin_global,nedge_ipoin,ltypb_global,&
                lnodb_global,lboch_global,ledbo_global,ledge_npoin_global(ipoin)%l,&
                list_ordered_boundaries,lperm,iview,0_ip)
if(ipoin==25748)print*,'salgo dod_order_boundaries',nboun_local

           !------------------------------------------------------
           ! 
           ! Selected list of extension's points
           !
           !------------------------------------------------------
           allocate ( lpoex(nnext,nboun_local) , stat = istat )
           allocate ( lextp(nnext,nboun_local) , stat = istat )
           do inext = 1,nnext
              do ii = 1,nboun_local
                 lpoex(inext,ii) = 0
                 lextp(inext,ii) = 0
              end do
           end do
           !
           ! Calculate of the normal to the each boundary
           !
           orien(1:3) = 0.0_rp
           do kboun = 1,nboun_local
              iboun = list_ordered_boundaries(kboun)
              !if(lboch_global(iboun)==BOEXT)then !!!!!OJOOOOOOOOOOOOO-nohacefalta
                 call dod_order_boundaries(&
                      2_ip,ipoin,mnodb,nboun_local,npoin_global,nedge_ipoin,ltypb_global,&
                      lnodb_global,lboch_global,ledbo_global,ledge_npoin_global(ipoin)%l,&
                      list_ordered_boundaries,lperm,iview,kboun)
                 call dod_extension_elements_3d(1_ip,ipoin,lnodb_global(lperm(1,iview),iboun),&
                      lnodb_global(lperm(2,iview),iboun),normx(kboun),         &
                      normy(kboun),normz(kboun),orien,nnext,lpoex(1,kboun),    &
                      0_ip,0_ip,0_ip,0_ip)
                 orien(1)     = orien(1) + normx(kboun)
                 orien(2)     = orien(2) + normy(kboun)
                 orien(3)     = orien(3) + normz(kboun)
                 !if(kboun>1)then
                 !   escalar_prod = normx(kboun)*normx(kboun-1) + normy(kboun)*normy(kboun-1) + normz(kboun)*normz(kboun-1) 
                 !end if

if(ipoin==10023)print*,'NORMALESSSSS',kboun,normz(kboun),normy(kboun),normz(kboun)
              !end if
           end do
           norma = orien(1) * orien(1) + orien(2) * orien(2) + orien(3) * orien(3)
           norma = sqrt(norma)
           if( norma == 0.0_rp ) then
              print*,'norma nulaaaaa',ipoin,orien(1:3),nboun_local
              stop
           else
              do idime = 1,ndime
                 orien(idime) = orien(idime) / norma
              end do
           end if
           !
           ! Eliminate for each boundary the candidates out of the "parabrisas"
           !
           do kboun = 1,nboun_local
              iboun = list_ordered_boundaries(kboun)
              !if(lboch_global(iboun)==BOEXT .or.lboch_global(iboun)==BOEXF )then
                 if(ipoin==7719)print*,'KBOUNNNNNNNNNNNNNN',kboun,lnodb_global(:,iboun),lnodb_global(lperm(1,iview),iboun),lnodb_global(lperm(2,iview),iboun)
                 call dod_order_boundaries(&
                      2_ip,ipoin,mnodb,nboun_local,npoin_global,nedge_ipoin,ltypb_global,&
                      lnodb_global,lboch_global,ledbo_global,ledge_npoin_global(ipoin)%l,&
                      list_ordered_boundaries,lperm,iview,kboun)
                 inode = 0
                 do inext = 1,nnext
                    lpoex(inext,kboun) = lnext(inext)
                 end do
                 call dod_extension_elements_3d(2_ip,ipoin,lnodb_global(lperm(1,iview),iboun),&
                      lnodb_global(lperm(2,iview),iboun),normx(kboun),normy(kboun),&
                      normz(kboun),orien,nnext,lpoex(1,kboun),isubd,jsubd,kboun,nboun_local)
              !end if
           end do
           !
           ! Hago comun denominador de la lista lpoex
           ! lcont(:,1,iboun)= el poexs
           ! lcont(:,2,iboun)= el numero de veces que ha salido
           !
           allocate ( lcont(nnext,2,nboun_local),stat=istat)
           do ii = 1,nnext
              do jj=1,2
                 do kk= 1,nboun_local
                    lcont(ii,jj,kk)=0
                 end do
              end do
           end do
           npoex = 0
           kboun = 0
           flag=0
           do while (kboun<nboun_local)
              kboun = kboun + 1
              do inext =1,nnext
                 cont  = 1
                 lcont(inext,1,kboun) = lpoex(inext,kboun)
                 lcont(inext,2,kboun) = cont
                 jboun = 0
                 do while (jboun<nboun_local )
                    jboun = jboun + 1
                    if(jboun/=kboun) then
                       jnext = 0
                       do while (jnext < nnext)
                          jnext = jnext + 1
                          if(lpoex(inext,kboun) == lpoex(jnext,jboun) .and. lpoex(inext,kboun) /= 0) then
                             cont = cont + 1
                             lcont(inext,2,kboun) = cont
                             lcont(jnext,1,jboun) = lpoex(inext,kboun)
                             if(cont==nboun_local)flag=1
                             jnext = nnext 
                          end if
                       end do
                    end if
                 end do
              end do
              if(flag==1)then
                 do jboun=1,nboun_local
                    do inext =1,nnext 
                       do jnext=1,nnext
                          if( lcont(inext,1,kboun) == lcont(jnext,1,jboun) )then
                             lcont(jnext,2,jboun)=lcont(inext,2,kboun)
                          end if
                       end do
                    end do
                 end do
                 kboun=nboun_local!!!OJO
              end if
           end do
           if(ipoin==7719)print*,'sigooooo-1'
           !
           ! Me creo la lista lextp que recopila los poexs comunes a los ibouns:lextp
           !
           tetex = 0
           do kboun=1,nboun_local
              npoex = 0
              do inext =1,nnext
                 if(lcont(inext,2,kboun)==nboun_local)then
                    npoex = npoex + 1
                    lextp(npoex,kboun) = lcont(inext,1,kboun)
                 end if
              end do
              if(npoex==0) then
                 tetex     = tetex+1
                 ii = 0
                 do while(ii<nboun_local)
                    ii  = ii + 1 
                    inext = 0
                    do while(inext <nnext)
                       inext = inext + 1
                       if(lcont(inext,2,kboun)== nboun_local-ii )then
                          npoex = npoex + 1
                          lextp(npoex,kboun) = lcont(inext,1,kboun)
                       end if
                    end do
                    if(npoex>0) then
                       ii = nboun_local!!!!!!!!!!!!Me voy a por otro contorno
                    end if
                 end do
              end if
           end do
           if(ipoin==7719)then
              print*,'sigooooo-2',tetex
              do kboun=1,nboun_local
                 print*,'kboun----------------------------',lextp(1:10,kboun)
              end do
           end if
           if(tetex>0)then
              do kboun=1,nboun_local
                 para(kboun) = 0
              end do
              do kboun=1,nboun_local
                 if(para(kboun)== 0 )then
                    para(kboun) = kboun
                    ipoex = lextp(1,kboun)
                    do jboun=1,nboun_local
                       if(para(jboun)== 0)then
                          jnext = 1
                          do while(lextp(jnext,jboun)/=0 .and. jnext<nnext)                
                             if(lextp(jnext,jboun)==ipoex) then
                                para(jboun)= kboun 
                             end if
                             jnext = jnext + 1
                          end do
                       end if
                    end do
                 end if
              end do
              allocate ( pasa(nboun_local),stat=istat)
              do jboun=1,nboun_local
                 pasa(jboun)=0
              end do
              tetex = 0
              do jboun = 1,nboun_local
                 if(para(jboun) == jboun)then
                    cont  = 0
                    kboun = 1
                    tetex = tetex + 1
                    do while (kboun<= nboun_local )
                       if(para(kboun) == jboun ) then
                          cont = cont + 1
                          pasa(jboun)= 1
                          para(jboun)= cont
                          if(jboun/=kboun .and. pasa(kboun)==0)para(kboun)= -jboun
                       end if
                       kboun = kboun + 1
                    end do
                 end if
              end do
              deallocate(pasa,stat=istat)

           end if
           if(ipoin==7719)print*,'sigooooo-3'

           ! No hay 1 mismo poex para todos:creo grupos de contornos +o- 
           ! paralelos entre si: primero creo el grupo de contornos que van al mismo poexs
           allocate ( quali(npoex,nboun_local),stat=istat)
           allocate ( q_max(nboun_local),stat=istat)
           allocate ( indic(nboun_local),stat=istat)
           !
           !Detecto si hay capa limite en los hexas
           !
           hcandi =10000000.0_rp
           do kboun=1,nboun_local
              eboun = list_ordered_boundaries(kboun)
              if( ltypb_global(eboun) == 12)then
                 hcarac(kboun) = 0.0_rp
                 hleng(1) = 0.0_rp
                 hleng(2) = 0.0_rp
                 if(ipoin==448180)print*,'2',lnodb_global(1,eboun),lnodb_global(2,eboun),lnodb_global(3,eboun),lnodb_global(4,eboun)
                 do idime=1,3
                    hleng(1) = hleng(1)+ ((coord(idime,lnodb_global(2,eboun))-coord(idime,lnodb_global(1,eboun)))*&
                         (coord(idime,lnodb_global(2,eboun))-coord(idime,lnodb_global(1,eboun))))
                    hleng(2) = hleng(2)+((coord(idime,lnodb_global(4,eboun))-coord(idime,lnodb_global(1,eboun)))*&
                         (coord(idime,lnodb_global(4,eboun))-coord(idime,lnodb_global(1,eboun))))
                 end do

                 if(hleng(1)==0.0_rp .or.hleng(2)==0.0_rp)then
                    print*,'contorno con nodo 4 y 1 iguales o 2 y 1 iguales!'
                    stop
                 end if
                 hleng(1) = sqrt(hleng(1))
                 hleng(2) = sqrt(hleng(2))
                 if(hleng(1)<hleng(2))then
                    aux = hleng(1) 
                    hleng(1) = hleng(2)
                    hleng(2) = aux
                 end if
                 do idime=1,3
                    coorcg1(idime,kboun) = 0.0_rp
                    coorcg2(idime)  = 0.0_rp
                 end do
                 do inodb=1,4
                    do idime=1,ndime
                       coorcg1(idime,kboun) = coorcg1(idime,kboun) + coord(idime,lnodb_global(inodb,eboun))
                    end do
                 end do
                 !if(ipoin==448180)print*,'4',hleng(1),hleng(2),coorcg1(:,kboun)
                 do idime=1,ndime
                    coorcg1(idime,kboun) = coorcg1(idime,kboun)/4.0_rp
                 end do
                 if(hleng(1)/hleng(2) > 0.9_rp*sqrt(hleng(1)*hleng(2)) )  caplim = .true.
                 jelem   = lboel_global(eboun)
                 if(jelem<1_ip)then
                    print*,'no encuentro elemento para comprobacion de capalim'
                    stop
                 end if
                 do inode=1,8_ip
                    cont = 0_ip
                    do jnode=1,4_ip
                       if(lnodb_global(jnode,eboun) /= lnods(inode,jelem)) cont = cont + 1
                    end do
                    if(cont==4)then
                       do idime=1,ndime
                          coorcg2(idime)=coorcg2(idime)+coord(idime,lnods(inode,jelem))
                          !if(ipoin==448180)print*,'lnods del jelem',idime,lnods(inode,jelem)
                       end do
                    end if
                 end do
                 !if(ipoin==448180)print*,coorcg2
                 do idime=1,ndime
                    coorcg2(idime) = coorcg2(idime)/4.0_rp
                    hcarac(kboun) = hcarac(kboun) + ( coorcg2(idime)-coorcg1(idime,kboun) )*( coorcg2(idime)-coorcg1(idime,kboun) )
                 end do
                 hcarac(kboun) = hcarac(kboun)**0.5_rp
                 !if(ipoin==448180)print*,hcarac(kboun)
                 do inext=1,npoex
                    distex(kboun,inext) = 0.0_rp
                 end do
              end if
              hcandi = min(hcarac(kboun),hcandi)
              if(ipoin==448180)print*,hcandi
           end do
           !---------------------------------------------------------
           !
           ! Construyo los tetras a base del contorno+poex
           !
           !---------------------------------------------------------
           !
           ! Un solo poexs: ningun otro elem de ext asociado
           !
           if(tetex == 0) then
              ielem_aux = 0
              nutri = 0
              do kboun=1,nboun_local
                 if(ipoin==10023 .or. ipoin==25748)print*,'tetex-0',kboun
                 eboun = list_ordered_boundaries(kboun)
                 call dod_order_boundaries(&
                      2_ip,ipoin,mnodb,nboun_local,npoin_global,nedge_ipoin,ltypb_global,&
                      lnodb_global,lboch_global,ledbo_global,ledge_npoin_global(ipoin)%l,&
                      list_ordered_boundaries,lperm,iview,kboun)

                 !if(ipoin==10023 .or. ipoin==25748)print*,'tetex-0',lnodb_global(lperm(2,iview),eboun),lnodb_global(lperm(1,iview),eboun)
                 q = 1000
                 do inext = 1,npoex
                    ipext = lextp(inext,kboun)
                    if(ipext/=0) then
                       call dod_extens_qual3d(2_ip,ipoin,lnodb_global(lperm(1,iview),eboun),&  
                            lnodb_global(lperm(2,iview),eboun),ipext,kappa,asrad,q,sign,dumm5)
                       quali(inext,kboun)=kappa
                 if(ipoin==7719)print*,'lleno-quali',kboun,ipext,quali(inext,kboun),sign
                       if( ltypb_global(eboun) == 12)then
                          if (caplim)then !me extiendo desde un elemento alargado

                             do idime=1,ndime
                                distex(kboun,inext) = distex(kboun,inext) + (coorcg1(idime,kboun) - coord(idime,ipext))* (coorcg1(idime,kboun) - coord(idime,ipext))
                             end do
                             distex(kboun,inext)=distex(kboun,inext)**0.5_rp
                             distex(kboun,inext)=distex(kboun,inext)-hcandi

                             call dod_extens_qual3d(2_ip,ipoin,lnodb_global(lperm(2,iview),eboun),&  
                                  lnodb_global(lperm(3,iview),eboun),ipext,kappa,asrad,q,sign,dumm5)
                             quali(inext,kboun) = quali(inext,kboun)+kappa
                             quali(inext,kboun) = quali(inext,kboun) * 0.5_rp

                          else
                             call dod_extens_qual3d(2_ip,ipoin,lnodb_global(lperm(2,iview),eboun),&  
                                  lnodb_global(lperm(3,iview),eboun),ipext,kappa,asrad,q,sign,dumm5)
                             quali(inext,kboun) = quali(inext,kboun)+kappa
                             quali(inext,kboun) = quali(inext,kboun) * 0.5_rp
                          end if
                       end if
                    else
                       quali(inext,kboun)=100000.0_rp
                    end if
                 end do
              end do
              !
              !choose the poexts in function on optimal qualitly
              !   QUIERO ORDENARLO DE MENOR A MAYORRRRRR

              if(ipoin==448180)print*,'b',caplim
              if(caplim)then
!!$REPASAR ORDENARLOOOOOO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 do kboun=1,nboun_local
                    do inext=1,npoex-1
                       if(lextp(inext,kboun)/=0)then
                          do jnext=inext+1,npoex
                             if(distex(kboun,inext)>distex(kboun,jnext) )then
                                tmp = lextp(inext,kboun)
                                lextp(inext,kboun)= lextp(jnext,kboun)
                                lextp(jnext,kboun)= tmp
                                qaux = distex(kboun,inext)
                                distex(kboun,inext) = distex(kboun,jnext) 
                                distex(kboun,jnext) = qaux
                             end if
                          end do
                       end if
                    end do
                 end do
              else
                 do kboun=1,nboun_local
                    do inext=1,npoex-1
                       if(lextp(inext,kboun)/=0)then
                          do jnext=inext+1,npoex
                             if(quali(inext,kboun)>quali(jnext,kboun) )then
                                tmp = lextp(inext,kboun)
                                lextp(inext,kboun)= lextp(jnext,kboun)
                                lextp(jnext,kboun)= tmp
                                qaux = quali(inext,kboun)
                                quali(inext,kboun)= quali(jnext,kboun)
                                quali(jnext,kboun)= qaux
                             end if
                          end do
                       end if
                    end do
                 end do
              end if
              do kboun=1,nboun_local
                 poext=lextp(1,kboun)   !cogiendo la menor de cada kboun....
                 !poext=lextp(npoex,kboun)   !cogiendo la MALLOR OJOOOO de cada kboun....
                 indic(kboun) = 1
                 do jboun=1,nboun_local
                    do inext=1,npoex
                       if(lextp(inext,jboun)==poext)then
                          if(inext>indic(kboun))indic(kboun) = inext !miro cual me da 
                          !if(inext<indic(kboun))indic(kboun) = inext !miro cual me da 
                       end if                                        !indic mayor:i.e.
                    end do                                           !la peor calidad
                 end do
              end do
              aux=nboun_local+1
              do kboun=1,nboun_local
                 if(indic(kboun)<aux)then !!ESTO ES LO QUE TIENE QUE SER!!
                    !if(indic(kboun)>aux)then!!!!!!!!!!!OJOOOOOOOOOOOOOOO PARA PILLAR MAL CALIDAD
                    aux = indic(kboun)
                    poext = lextp(1,kboun)
                 end if
              end do

              !!OJOOO..GUILLAUME......
              !if( ltypb_global(eboun) == 12)then
              !   kelem = nboun_local
              !   knode = 5
              !end if
              !if( ltypb_global(eboun) == 10)then
              !   kelem   =  nboun_local
              !   knode   =  4
              !end if
              !lo aloco todo a knode=5 por si hay piramides...
              kelem = nboun_local
              knode = 5
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LNODS','dod_extension',lpext(number_fringe_nodes) % lnods,knode,kelem,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LTYPE','dod_extension',lpext(number_fringe_nodes) % ltype,kelem      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LMATE','dod_extension',lpext(number_fringe_nodes) % lmate,kelem      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LELEZ','dod_extension',lpext(number_fringe_nodes) % lelez,kelem      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LESUB','dod_extension',lpext(number_fringe_nodes) % lesub,kelem      ,'DO_NOT_INITIALIZE')
              do kboun =1,nboun_local
                 eboun = list_ordered_boundaries(kboun)
                 call dod_order_boundaries(&
                      2_ip,ipoin,mnodb,nboun_local,npoin_global,nedge_ipoin,ltypb_global,&
                      lnodb_global,lboch_global,ledbo_global,ledge_npoin_global(ipoin)%l,&
                      list_ordered_boundaries,lperm,iview,kboun)
                 if( ltypb_global(eboun) == 12)then
                    imatn          = lmatn_dod(poext)
                    ipoiz          = lpoiz_dod(poext)
                    lpext(number_fringe_nodes) % ltype(kboun) = PYR05
                    lpext(number_fringe_nodes) % lmate(kboun) = imatn
                    lpext(number_fringe_nodes) % lelez(kboun) = ipoiz
                    lpext(number_fringe_nodes) % lesub(kboun) = isubd
                    poex1 = lnodb_global(lperm(1,iview),eboun)
                    poex2 = lnodb_global(lperm(2,iview),eboun)
                    poex3 = lnodb_global(lperm(3,iview),eboun)
                    poex4 = poext
                    ! call dod_qauxp(ipoin,poex1,poex2,poex3,poext)
                    !call dod_extens_qual3d(3_ip,ipoin,poex1,poex2,poex3,kappa,asrad,q,sign,dumm5)
                    lpext(number_fringe_nodes) % lnods(1,kboun) = ipoin
                    lpext(number_fringe_nodes) % lnods(2,kboun) = poex1
                    lpext(number_fringe_nodes) % lnods(3,kboun) = poex2                 
                    lpext(number_fringe_nodes) % lnods(4,kboun) = poex3                 
                    lpext(number_fringe_nodes) % lnods(5,kboun) = poext 
                    nelem_dod = nelem_dod + 1 

                 end if
                 if( ltypb_global(eboun) == 10)then
                    imatn   =  lmatn_dod(poext)
                    ipoiz   =  lpoiz_dod(poext)
                    !!poex1   =  lnodb_global(lperm(2,iview),eboun)
                    !!poex2   =  lnodb_global(lperm(1,iview),eboun)
                    poex1   =  lnodb_global(lperm(1,iview),eboun)   !OJOOOOOOOOOOOOOO (en estenosis me doy cuenta...)
                    poex2   =  lnodb_global(lperm(2,iview),eboun)   !OJOOOOOOOOOOOOOO
                    poex3   =  poext
                    lpext(number_fringe_nodes) % ltype(kboun) = TET04
                    lpext(number_fringe_nodes) % lmate(kboun) = imatn
                    lpext(number_fringe_nodes) % lelez(kboun) = ipoiz
                    lpext(number_fringe_nodes) % lesub(kboun) = isubd
                    !call dod_quaux(ipoin,poex1,poex2,poex3)
                    if(ipoin==7719 )print*,'entrooooooo',ipoin,poex1,poex2,poex3,q
                    call dod_extens_qual3d(3_ip,ipoin,poex1,poex2,poex3,kappa,asrad,q,sign,dumm5)
                    if(ipoin==10023 .or. ipoin==25748 )print*,'salggoooooo',ipoin,poex1,poex2,poex3,q,sign
                    if(sign ==1_ip )print*,'CAMBIOOOOOOOOOO DE SIGNOOOOOOOOOOO',ipoin,poex1,poex2,poex3,eboun
!.and. lboch_global(eboun)/=BOEXF )then
                      ! print*,'CAMBIOOOOOOOOOO DE SIGNOOOOOOOOOOO',ipoin,poex1,poex2,poex3,eboun
                      ! stop
                    !end if
                    lpext(number_fringe_nodes) % lnods(1,kboun) = ipoin
                    lpext(number_fringe_nodes) % lnods(2,kboun) = poex1
                    lpext(number_fringe_nodes) % lnods(3,kboun) = poex2                 
                    lpext(number_fringe_nodes) % lnods(4,kboun) = poex3                 
                    nelem_dod = nelem_dod + 1 
                    if(ipoin==7719 )print*,lpext(number_fringe_nodes) % lnods(1:4,kboun)
                 end if
!!$                 do inode = 1,nnode(lpext(number_fringe_nodes) % ltype(kboun))
!!$                    ipoex = lpext(number_fringe_nodes) % lnods(inode,kboun)
!!$                    do idime=1,ndime
!!$                       coorcg(idime) = coorcg(idime) + coord(idime,ipoex)
!!$                    end do
!!$                 end do
!!$                 call elsest(&
!!$                      2_ip,jsubd,ielse,mnode,ndime,neighbor_subdomain % npoin,neighbor_subdomain % nelem,&
!!$                      nnode(1:),neighbor_subdomain % lnods,neighbor_subdomain % ltype,&
!!$                      ltopo,neighbor_subdomain % coord,coorcg,relse,jelem,&
!!$                      shape_dod,deriv_dod,xieta,dumm1)
!!$                 lpext(number_fringe_nodes) % lmate(kboun) = lmate(jelem)

                 if(lpext(number_fringe_nodes) % lnods(1,kboun)==0 .or. lpext(number_fringe_nodes) % lnods(2,kboun) ==0.or. &
                      lpext(number_fringe_nodes) % lnods(3,kboun)==0 .or. lpext(number_fringe_nodes) % lnods(4,kboun) ==0)then
                    print*,ipoin,'cerrrrrrrrrrrrrrrrrrro',lpext(number_fringe_nodes) % lnods(1:4,kboun)
                    stop
                 end if

              end do
              !if(ipoin==7719  )stop

           else
              if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'tetex no 0'

              !
              ! Mas de un poex's: necesita de mas elems de exts para cerrar:
              !            -tetras con ipoin-contornos-poexs :ielem_aux
              !            -tetras con ipoin-poexs: 'catalan': nutri

              do kboun=1,nboun_local
                 if(para(kboun)>0)then
                    poext=lextp(1,kboun)
                    lext2(kboun) = poext
                 else if(para(kboun)< 0) then
                    poext = lextp(1,abs(para(kboun)))
                    lext2(kboun)= poext
                 end if
              end do
              ielem_aux = 0
              do kboun=1,nboun_local
                 if(kboun>1) then 
                    if(lext2(kboun)/=lext2(kboun-1))then
                       ielem_aux = ielem_aux + 1
                    end if
                 end if
                 if(kboun==nboun_local)then
                    if(lext2(kboun)/=lext2(1))then
                       ielem_aux = ielem_aux + 1
                    end if
                 end if
              end do
              if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux',ielem_aux

              if(ielem_aux==2)then
                 ipoex = lext2(1)
                 kboun=1
                 do while (lext2(kboun)==ipoex .and. kboun<nboun_local)
                    kboun = kboun +  1
                    jpoex = lext2(kboun)
                 end do
                 call dod_order_boundaries(&
                      2_ip,ipoin,mnodb,nboun_local,npoin_global,nedge_ipoin,ltypb_global,&
                      lnodb_global,lboch_global,ledbo_global,ledge_npoin_global(ipoin)%l,&
                      list_ordered_boundaries,lpern,kview,kboun)
                 if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux2-a'
                 ii = kboun
                 do while (lext2(ii)==jpoex .and. ii<nboun_local)
                    ii = ii + 1
                 end do
                 if(abs(kboun - ii) == 1) then
                    qboun = kboun
                 else
                    qboun = ii
                 end if
                 call dod_order_boundaries(&
                      2_ip,ipoin,mnodb,nboun_local,npoin_global,nedge_ipoin,ltypb_global,&
                      lnodb_global,lboch_global,ledbo_global,ledge_npoin_global(ipoin)%l,&
                      list_ordered_boundaries,lpert,jview,qboun)
                 poex1=lext2(1)
                 poex2=lext2(kboun)
                 poex3=lnodb_global(lperm(2,kview),list_ordered_boundaries(kboun))
                 if( ltypb_global(list_ordered_boundaries(kboun)) == 12)&
                      poex3 =lnodb_global(lperm(3,kview),list_ordered_boundaries(kboun))
                 if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux2-b'

                 if (qboun==kboun .or. lext2(1)/=lext2(nboun_local)) then
                    poex4 = lnodb_global(lpert(1,jview),list_ordered_boundaries(qboun))
                    if( ltypb_global(list_ordered_boundaries(qboun)) == 12)&
                         poex4 =lnodb_global(lperm(2,jview),list_ordered_boundaries(qboun))
                    if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux2-c'
                 else 
                    poex4 = lnodb_global(lpert(2,jview),list_ordered_boundaries(qboun))
                    if( ltypb_global(list_ordered_boundaries(qboun)) == 12)&
                         poex4 =lnodb_global(lperm(3,jview),list_ordered_boundaries(qboun))
                    if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux2-d'
                 end if
                 !Quizas con el grafo de lados pueda localizar los puntos de contorno con los que debo formar tetras facilmente...
                 !porque si mezclo tris y quads no se bien que nodo de contorno es el compartido
!!$                 if (nboun_local==3)poex4 = lnodb(lpert(1,jview),list_ordered_boundaries(ipoin,qboun))
                 ! call dod_qualit(ipoin,poex1,poex3,poex2,q1,sign)!!!!comprobar
                 call dod_extens_qual3d(2_ip,ipoin,poex1,poex3,poex2,kappa,asrad,q,sign,orien)
                 if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux2-e'
                 ! call dod_qualit(ipoin,poex1,poex2,poex4,q2,sign)!!!!...
                 call dod_extens_qual3d(2_ip,ipoin,poex1,poex2,poex4,kappa,asrad,q,sign,orien)
!!$                 !call dod_condex(ipoin,poex1,poex3,poex2,q1,asrad)!!!!comprobar
!!$                 !call dod_condex(ipoin,poex1,poex2,poex4,q2,asrad)!!!!...
!!$                 call dod_orient(ipoin,poex3,poex1,poex2,orien,corre)
                 if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux2-f'
                 call dod_extens_qual3d(4_ip,ipoin,poex3,poex1,poex2,kappa,asrad,q,sign,orien)
!!$                 call dod_orient(ipoin,poex4,poex1,poex2,orien,corre)

                 call dod_extens_qual3d(4_ip,ipoin,poex4,poex1,poex2,kappa,asrad,q,sign,orien)
                 if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux2-g'
                 if(q1>50.0_rp .or. q2>50.0_rp .or. sign==0_ip)then
                    q=max(q1,q2)
                    q_old= q
                    ii = 0
                    do inext = 1,nnext
                       if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux2-h-inext',inext
                       ipoex = lpoex(inext,1)
                       do jnext =1,nnext
                          if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux2-h-jnext',jnext
                          jpoex = lpoex(jnext,kboun)
                          if(ipoex/=0 .and. jpoex/=0 .and. ipoex/=jpoex)then
                             if(nboun_local==3)then  !!!!OJOOOOO
                                if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux2-i-nboun_local3',jnext
!!$                                call dod_qualit(ipoin,ipoex,&
!!$                                     lnodb_global(lpern(2,kview),list_ordered_boundaries(kboun)),&
!!$                                     jpoex,q1,sign)
                                call dod_extens_qual3d(2_ip,ipoin,ipoex,&
                                     lnodb_global(lpern(2,kview),list_ordered_boundaries(kboun)),&
                                     jpoex,kappa,asrad,q1,sign,orien)
!!$                                call dod_qualit(ipoin,ipoex,jpoex,&
!!$                                     lnodb_global(lpert(1,kview),list_ordered_boundaries(kboun)),q2,sign)
                                call dod_extens_qual3d(2_ip,ipoin,ipoex,jpoex,&
                                     lnodb_global(lpert(1,kview),list_ordered_boundaries(kboun)), &
                                     kappa,asrad,q2,sign,orien)
                             else
                                if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux2-i-nboun_localNO3-a',jnext,ipoin,ipoex,lnodb(lpern(2,kview),list_ordered_boundaries(kboun)),jpoex
                                call dod_extens_qual3d(2_ip,ipoin,ipoex,&
                                     lnodb_global(lpern(2,kview),list_ordered_boundaries(kboun)),&
                                     jpoex,kappa,asrad,q1,sign,orien)
                                if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux2-i-nboun_localNO3-b',jnext,ipoin,ipoex,jpoex,lnodb(lpert(2,jview),list_ordered_boundaries(qboun))
                                call dod_extens_qual3d(2_ip,ipoin,ipoex,jpoex,& 
                                     lnodb_global(lpert(2,jview),list_ordered_boundaries(qboun)),kappa,asrad,q2,sign,orien)
                             end if
                             q = max(q1,q2)
                             if(q<q_old) then
                                lext2(1)     = ipoex !aqui modifico lext2 para conseguir la mejor calidad
                                lext2(kboun) = jpoex 
                                q_old = q
                                ii = ii + 1
                             end if
                          end if
                       end do
                    end do
                    if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux2-j'

                    ielem_aux = 0
                    do kboun=1,nboun_local
                       if(kboun>1) then 
                          if(lext2(kboun)/=lext2(kboun-1))then
                             ielem_aux= ielem_aux + 1
                          end if
                       end if
                       if(kboun==nboun_local)then
                          if(lext2(kboun)/=lext2(1))then
                             ielem_aux= ielem_aux + 1
                          end if
                       end if
                    end do
                    tetex = 1
                    ipoex=lext2(1)
                    ii = 0
                    do while(ii<=nboun_local-1)
                       ii = ii + 1
                       flag = 0
                       do jj=1,ii-1
                          if(lext2(ii)==lext2(jj))flag=1
                       end do
                       if(lext2(ii)/=ipoex .and. flag==0)then
                          tetex = tetex + 1
                          ipoex = lext2(ii)
                       end if
                    end do
                    flag = 0
                    do jj =1,nboun_local-1
                       if(lext2(nboun_local)==lext2(jj))flag=1  
                    end do
                    if(flag==0)then
                       tetex = tetex + 1
                    end if
                 end if
              end if
              if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'ielem_aux2-k'
              if(ielem_aux>2)then
                 allocate ( lqaux(nboun_local),stat=istat)
                 do kboun=1,nboun_local
                    lqaux(kboun) = 0
                 end do
                 flag = 0
                 do while(flag==0)
                    ielem_old = ielem_aux
                    ielem_aux = 0
                    do kboun=1,nboun_local
                       eboun = list_ordered_boundaries(kboun)
                       call dod_order_boundaries(&
                            2_ip,ipoin,mnodb,nboun_local,npoin_global,nedge_ipoin,ltypb_global,&
                            lnodb_global,lboch_global,ledbo_global,ledge_npoin_global(ipoin)%l,&
                            list_ordered_boundaries,lperm,iview,kboun)
                       if(kboun>1 ) then 
                          if(lext2(kboun)/=lext2(kboun-1))then
                             ielem_aux= ielem_aux + 1
                             poex1 = lext2(kboun)
                             poex2 = lext2(kboun-1)
                             poex3 = lnodb_global(lperm(2,iview),eboun)
                             if( ltypb_global(eboun) == 12)poex3 =lnodb_global(lperm(3,iview),eboun)
                             ! call dod_qualit(ipoin,poex1, poex2 ,poex3,lqaux(ielem_aux),sign)
                             call dod_extens_qual3d(2_ip,ipoin,poex1,poex2,poex3,&
                                  kappa,asrad,lqaux(ielem_aux),sign,orien)
                          end if
                       end if
                       if(kboun==nboun_local)then
                          if(lext2(kboun)/=lext2(1))then
                             ielem_aux              = ielem_aux + 1
                             poex1 = lnodb_global(lperm(1,iview),eboun)
                             if( ltypb_global(eboun) == 12)poex1 =lnodb_global(lperm(2,iview),eboun)
                             poex2 = lext2(1)    
                             poex3 = lext2(kboun)  
                             call dod_extens_qual3d(2_ip,ipoin,poex1,poex2,poex3,&
                                  kappa,asrad,lqaux(ielem_aux),sign,orien)
                             ! call dod_qualit(ipoin,poex1, poex2 ,poex3,lqaux(ielem_aux),sign)
                          end if
                       end if
                    end do
                    jj = 0
                    do while(jj<ielem_aux)
                       jj  = jj  + 1
                       if(lqaux(jj)>50)then
                          q  =  lqaux(jj)
                          jj = ielem_aux
                          cont=0
                          poex1=lext2(1)
                          do ii=1,nboun_local
                             eboun = list_ordered_boundaries(ii)
                             call dod_order_boundaries(&
                                  2_ip,ipoin,mnodb,nboun_local,npoin_global,nedge_ipoin,ltypb_global,&
                                  lnodb_global,lboch_global,ledbo_global,ledge_npoin_global(ipoin)%l,&
                                  list_ordered_boundaries,lperm,iview,ii)
                             if(ii>1 ) then 
                                if(lext2(ii)/=lext2(ii-1))then
                                   cont = cont + 1
                                   if(cont ==1)then
                                      kboun =  ii
                                      pobok =  lnodb_global(lperm(2,iview),list_ordered_boundaries(ii))
                                      if( ltypb_global(list_ordered_boundaries(ii)) == 12)pobok =  lnodb_global(lperm(3,iview),list_ordered_boundaries(ii))
                                      poexk =  lext2(ii)
                                   else if(cont==2)then
                                      qboun =  ii
                                      poboq =  lnodb_global(lperm(2,iview),list_ordered_boundaries(ii))
                                      if( ltypb_global(list_ordered_boundaries(ii)) == 12)poboq =  lnodb_global(lperm(3,iview),list_ordered_boundaries(ii))
                                      poexq =  lext2(ii)
                                   else if(cont==3 .and.ielem_old==3)then
                                      rboun =  ii
                                      pobor =  lnodb_global(lperm(1,iview),list_ordered_boundaries(ii))
                                      if( ltypb_global(list_ordered_boundaries(ii)) == 12)pobor =  lnodb_global(lperm(2,iview),list_ordered_boundaries(ii))
                                      poexr =  lext2(ii)
                                   else if(cont==3 .and.ielem_old==4)then
                                      rboun =  ii
                                      pobor =  lnodb_global(lperm(2,iview),list_ordered_boundaries(ii))
                                      if( ltypb_global(list_ordered_boundaries(ii)) == 12)pobor =  lnodb_global(lperm(3,iview),list_ordered_boundaries(ii))
                                      poexr =  lext2(ii)
                                   else if(cont==4 .and.ielem_old==4)then
                                      sboun =  ii
                                      pobos =  lnodb_global(lperm(1,iview),list_ordered_boundaries(ii))
                                      if( ltypb_global(list_ordered_boundaries(ii)) == 12)pobos =  lnodb_global(lperm(2,iview),list_ordered_boundaries(ii))
                                      poexs =  lext2(ii)
                                   end if
                                end if
                             end if
                             if(ii==nboun_local)then
                                if(lext2(kboun)/=lext2(1))then
                                   if(ielem_old== 3)then
                                      rboun =  ii
                                      pobor =  lnodb_global(lperm(1,iview),list_ordered_boundaries(ii))
                                      if( ltypb_global(list_ordered_boundaries(ii)) == 12)pobor =  lnodb_global(lperm(2,iview),list_ordered_boundaries(ii))
                                      poexr =  lext2(ii)
                                      sboun =  0
                                      pobos =  0
                                   else if(ielem_old==4)then
                                      sboun =  ii
                                      pobos =  lnodb_global(lpern(1,iview),list_ordered_boundaries(ii))
                                      if( ltypb_global(list_ordered_boundaries(ii)) == 12)pobos =  lnodb_global(lperm(2,iview),list_ordered_boundaries(ii))
                                      poexs =  lext2(ii)
                                   end if
                                end if
                             end if
                          end do
                          call dod_extension_correc(q,tetex,ielem_aux,nboun_local,nnext,ipoin,&
                               pobok,poboq,pobor,pobos,poex1,poexk,&
                               poexq,poexr,qboun,kboun,rboun,sboun,lpoex,lext2)
                          ielem_aux = 0
                          do kboun=1,nboun_local
                             if(kboun>1) then 
                                if(lext2(kboun)/=lext2(kboun-1))then
                                   ielem_aux  = ielem_aux + 1
                                end if
                             end if
                             if(kboun==nboun_local)then
                                if(lext2(kboun)/=lext2(1))then
                                   ielem_aux  = ielem_aux + 1
                                end if
                             end if
                          end do
                          if(ielem_aux==ielem_old)then
                             flag=1
                          else
                             ielem_old=ielem_aux
                          end if

                       end if
                    end do
                    flag = 1
                 end do
                 tetex = 1
                 ipoex=lext2(1)
                 ii = 0
                 do while(ii<=nboun_local-1)
                    ii = ii + 1
                    flag = 0
                    do jj=1,ii-1
                       if(lext2(ii)==lext2(jj))flag=1
                    end do
                    if(lext2(ii)/=ipoex .and. flag==0)then
                       tetex = tetex + 1
                       ipoex = lext2(ii)
                    end if
                 end do
                 flag = 0
                 do ii =1,nboun_local-1
                    if(lext2(nboun_local)==lext2(jj))flag=1  
                 end do
                 if(flag==0)then
                    tetex = tetex + 1
                 end if
                 deallocate(lqaux,stat=istat)
              end if
              kelem   =  nboun_local
              knode   =  4
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LNODS','dod_extension',lpext(number_fringe_nodes) % lnods,knode,kelem,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LTYPE','dod_extension',lpext(number_fringe_nodes) % ltype,kelem      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LMATE','dod_extension',lpext(number_fringe_nodes) % lmate,kelem      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LELEZ','dod_extension',lpext(number_fringe_nodes) % lelez,kelem      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LESUB','dod_extension',lpext(number_fringe_nodes) % lesub,kelem      ,'DO_NOT_INITIALIZE')
              do kboun =1,nboun_local
                 eboun = list_ordered_boundaries(kboun)
                 call dod_order_boundaries(&
                      2_ip,ipoin,mnodb,nboun_local,npoin_global,nedge_ipoin,ltypb_global,&
                      lnodb_global,lboch_global,ledbo_global,ledge_npoin_global(ipoin)%l,&
                      list_ordered_boundaries,lperm,iview,kboun)
                 if( ltypb_global(eboun) == 12)then
                    poex1 = lnodb_global(lperm(1,iview),eboun)
                    poex2 = lnodb_global(lperm(2,iview),eboun)
                    poex3 = lnodb_global(lperm(3,iview),eboun)
                    poex4 = lext2(kboun)
                    imatn          = lmatn_dod(poex4)
                    ipoiz          = lpoiz_dod(poex4)
                    lpext(number_fringe_nodes) % ltype(kboun) = PYR05
                    lpext(number_fringe_nodes) % lmate(kboun) = imatn
                    lpext(number_fringe_nodes) % lelez(kboun) = ipoiz
                    lpext(number_fringe_nodes) % lesub(kboun) = isubd
                    ! call dod_qauxp(ipoin,poex1,poex2,poex3,poext)
                    !call dod_extens_qual3d(3_ip,ipoin,poex1,poex2,poex3,kappa,asrad,q,sign,dumm5)
                    lpext(number_fringe_nodes) % lnods(1,kboun) = ipoin
                    lpext(number_fringe_nodes) % lnods(2,kboun) = poex1
                    lpext(number_fringe_nodes) % lnods(3,kboun) = poex2                 
                    lpext(number_fringe_nodes) % lnods(4,kboun) = poex3                 
                    lpext(number_fringe_nodes) % lnods(5,kboun) = poex4                 
                    nelem_dod = nelem_dod + 1 
                 end if
                 if( ltypb_global(eboun) == 10)then
                    !poex1 = lnodb_global (lperm(2,iview),eboun)
                    !poex2 = lnodb_global (lperm(1,iview),eboun)
                    poex1 = lnodb_global (lperm(1,iview),eboun) !!!!!!!OJOOOOOOO-estenosis
                    poex2 = lnodb_global (lperm(2,iview),eboun) !!!!!!!OJOOOOOOO
                    poex3 = lext2(kboun)
                    imatn          = lmatn_dod(poex3)
                    ipoiz          = lpoiz_dod(poex3)
                    lpext(number_fringe_nodes) % ltype(kboun) = TET04
                    lpext(number_fringe_nodes) % lmate(kboun) = imatn
                    lpext(number_fringe_nodes) % lelez(kboun) = ipoiz
                    lpext(number_fringe_nodes) % lesub(kboun) = isubd
                    call dod_extens_qual3d(3_ip,ipoin,poex1,poex2,poex3,kappa,asrad,q,sign,dumm5)
                    if(sign ==1_ip) print*,'CAMBIOOOOO SIGNO-1',ipoin,poex1,poex2,poex3,eboun
!.and. lboch_global(eboun)/=BOEXF )then
!                       print*,'CAMBIOOOOO SIGNO-1',ipoin,poex1,poex2,poex3,eboun
!                       stop
!                    end if
                    lpext(number_fringe_nodes) % lnods(1,kboun) = ipoin
                    lpext(number_fringe_nodes) % lnods(2,kboun) = poex1
                    lpext(number_fringe_nodes) % lnods(3,kboun) = poex2                 
                    lpext(number_fringe_nodes) % lnods(4,kboun) = poex3                 
                    nelem_dod = nelem_dod + 1 
                    if(lpext(number_fringe_nodes) % lnods(1,kboun)==0 .or. lpext(number_fringe_nodes) % lnods(2,kboun)==0 .or. &
                         lpext(number_fringe_nodes) % lnods(3,kboun)==0 .or. lpext(number_fringe_nodes) % lnods(4,kboun)==0 )then
                       print*,ipoin,'cerrrrrrrrrrrrrrrrrrro-b'
                       stop
                    end if
                 end if
!!$                 do inode = 1,nnode(lpext(number_fringe_nodes) % ltype(kboun))
!!$                    ipoex = lpext(number_fringe_nodes) % lnods(inode,kboun)
!!$                    do idime=1,ndime
!!$                       coorcg(idime) = coorcg(idime) + coord(idime,ipoex)
!!$                    end do
!!$                 end do
!!$                 call elsest(&
!!$                      2_ip,jsubd,ielse,mnode,ndime,neighbor_subdomain % npoin,neighbor_subdomain % nelem,&
!!$                      nnode(1:),neighbor_subdomain % lnods,neighbor_subdomain % ltype,&
!!$                      ltopo,neighbor_subdomain % coord,coorcg,relse,jelem,&
!!$                      shape_dod,deriv_dod,xieta,dumm1)
!!$                 lpext(number_fringe_nodes) % lmate(kboun) = lmate(jelem)
              end do

           end if
           if (tetex<2)then!!!este primer if lo puedo sacar....
              ielem_aux = 0
              nutri = 0
           else if(tetex>=2)then
              !
              !numero de tetras entre ipoin-2exts-bopo 
              !
              ielem_aux = 0
              do kboun=1,nboun_local
                 if(kboun>1) then 
                    if(lext2(kboun)/=lext2(kboun-1))then
                       ielem_aux= ielem_aux + 1
                    end if
                 end if
                 if(kboun==nboun_local)then
                    if(lext2(kboun)/=lext2(1))then
                       ielem_aux              = ielem_aux + 1
                    end if
                 end if
              end do
              !
              !numero de tetras para cerrar la tapa con catalan
              !
              nutri = 0
              tetex = 1
              ipoex=lext2(1)
              ii = 0
              do while(ii<=nboun_local-1)
                 ii = ii + 1
                 flag = 0
                 do jj=1,ii-1
                    if(lext2(ii)==lext2(jj))flag=1
                 end do
                 if(lext2(ii)/=ipoex .and. flag==0)then
                    tetex = tetex + 1
                    ipoex = lext2(ii)
                 end if
              end do
              flag = 0
              do ii =1,nboun_local-1
                 if(lext2(nboun_local)==lext2(jj))flag=1  
              end do
              if(flag==0)then
                 tetex = tetex + 1
              end if

              if(tetex>2) nutri = tetex - 2
              knode = 5
              kelem_aux = ielem_aux+nutri
              call memory_alloca(mem_servi(1:2,servi),'LNODS_AUX','dod_extension_3d',lnods_aux,knode,nboun_local,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LTYPE_AUX','dod_extension_3d',ltype_aux,nboun_local      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LMATE_AUX','dod_extension_3d',lmate_aux,nboun_local      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LELEZ_AUX','dod_extension_3d',lelez_aux,nboun_local      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LESUB_AUX','dod_extension_3d',lesub_aux,nboun_local      ,'DO_NOT_INITIALIZE')
              do ielem = 1, nboun_local
                 knode = nnode(lpext(number_fringe_nodes)%ltype(ielem))
                 do inode = 1,knode 
                    lnods_aux(inode,ielem) = lpext(number_fringe_nodes) % lnods(inode,ielem)
                 end do
                 ltype_aux(ielem)=lpext(number_fringe_nodes)%ltype(ielem)
                 lmate_aux(ielem)=lpext(number_fringe_nodes)%lmate(ielem)
                 lelez_aux(ielem)=lpext(number_fringe_nodes)%lelez(ielem)
                 lesub_aux(ielem)=lpext(number_fringe_nodes)%lesub(ielem)
              end do
              call memory_deallo(mem_servi(1:2,servi),'LPEXT(NUMBER_FRINGE_NODES) % LNODS','dod_extension_3d',lpext(number_fringe_nodes) % lnods)
              call memory_deallo(mem_servi(1:2,servi),'LPEXT(NUMBER_FRINGE_NODES) % LTYPE','dod_extension_3d',lpext(number_fringe_nodes) % ltype)
              call memory_deallo(mem_servi(1:2,servi),'LPEXT(NUMBER_FRINGE_NODES) % LMATE','dod_extension_3d',lpext(number_fringe_nodes) % lmate)
              call memory_deallo(mem_servi(1:2,servi),'LPEXT(NUMBER_FRINGE_NODES) % LELEZ','dod_extension_3d',lpext(number_fringe_nodes) % lelez)
              call memory_deallo(mem_servi(1:2,servi),'LPEXT(NUMBER_FRINGE_NODES) % LESUB','dod_extension_3d',lpext(number_fringe_nodes) % lesub)

              kelem = nboun_local + kelem_aux
              knode = 5

              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LNODS','dod_extension',lpext(number_fringe_nodes) % lnods,knode,kelem,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LTYPE','dod_extension',lpext(number_fringe_nodes) % ltype,kelem      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LMATE','dod_extension',lpext(number_fringe_nodes) % lmate,kelem      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LELEZ','dod_extension',lpext(number_fringe_nodes) % lelez,kelem      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LESUB','dod_extension',lpext(number_fringe_nodes) % lesub,kelem      ,'DO_NOT_INITIALIZE')
              if(nutri==0)print*,'--------------------ALOKAAAAAAAAAADA-4------',ipoin,nboun_local,ielem_aux,kelem_aux,tetex,nutri
              do ielem =1,nboun_local
                 knode = nnode(ltype_aux(ielem))
                 do inode=1,knode
                    lpext(number_fringe_nodes)%lnods(inode,ielem)=lnods_aux(inode,ielem)
                 end do
                 lpext(number_fringe_nodes)%ltype(ielem)=ltype_aux(ielem)
                 lpext(number_fringe_nodes)%lmate(ielem)=lmate_aux(ielem)
                 lpext(number_fringe_nodes)%lelez(ielem)=lelez_aux(ielem)
                 lpext(number_fringe_nodes)%lesub(ielem)=lesub_aux(ielem)
                 !nelem_dod = nelem_dod + 1 
              end do
              call memory_deallo(mem_servi(1:2,servi),'LNODS_AUX','dod_extension_3d',lnods_aux)
              call memory_deallo(mem_servi(1:2,servi),'LTYPE_AUX','dod_extension_3d',ltype_aux)
              call memory_deallo(mem_servi(1:2,servi),'LMATE_AUX','dod_extension_3d',lmate_aux)
              call memory_deallo(mem_servi(1:2,servi),'LELEZ_AUX','dod_extension_3d',lelez_aux)
              call memory_deallo(mem_servi(1:2,servi),'LESUB_AUX','dod_extension_3d',lesub_aux)
           end if  !!!OJOOOO
           if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'TETEXXXX-1',tetex,ielem_aux

           if (tetex ==2 .and. ielem_aux<3) then
              ielem_aux = 0
              ipoex = lext2(1)
              kboun=1
              do while (lext2(kboun)==ipoex .and. kboun<nboun_local)
                 kboun = kboun +  1
                 jpoex = lext2(kboun)
              end do
              call dod_order_boundaries(&
                   2_ip,ipoin,mnodb,nboun_local,npoin_global,nedge_ipoin,ltypb_global,&
                   lnodb_global,lboch_global,ledbo_global,ledge_npoin_global(ipoin)%l,&
                   list_ordered_boundaries,lpern,kview,kboun)
              ii = kboun
              do while (lext2(ii)==jpoex .and. ii<nboun_local)
                 ii = ii + 1
              end do
              if(abs(kboun - ii) == 1 .and. ii/=nboun_local) then
                 qboun = kboun
              else
                 qboun = ii
              end if
              call dod_order_boundaries(&
                   2_ip,ipoin,mnodb,nboun_local,npoin_global,nedge_ipoin,ltypb_global,&
                   lnodb_global,lboch_global,ledbo_global,ledge_npoin_global(ipoin)%l,&
                   list_ordered_boundaries,lpert,jview,qboun)
              ielem_aux = ielem_aux + 1
              poex1 = lext2(1)  
              poex2 = lnodb_global(lpern(2,kview),list_ordered_boundaries(kboun))
              if( ltypb_global(list_ordered_boundaries(kboun)) == 12)poex2 = lnodb_global(lperm(3,kview),list_ordered_boundaries(kboun))
              poex3 = lext2(kboun)
              imatn          = lmatn_dod(poex1)
              ipoiz          = lpoiz_dod(poex1)
              lpext(number_fringe_nodes) % ltype(nboun_local + ielem_aux) = TET04
              lpext(number_fringe_nodes) % lmate(nboun_local + ielem_aux) = imatn
              lpext(number_fringe_nodes) % lelez(nboun_local + ielem_aux) = ipoiz
              lpext(number_fringe_nodes) % lesub(nboun_local + ielem_aux) = isubd
              call dod_extens_qual3d(3_ip,ipoin,poex1,poex2,poex3,kappa,asrad,q,sign,dumm5)
              lpext(number_fringe_nodes) % lnods(1,nboun_local + ielem_aux) = ipoin
              lpext(number_fringe_nodes) % lnods(2,nboun_local + ielem_aux) = poex1 
              lpext(number_fringe_nodes) % lnods(3,nboun_local + ielem_aux) = poex2
              lpext(number_fringe_nodes) % lnods(4,nboun_local + ielem_aux) = poex3
              nelem_dod = nelem_dod + 1 

              if(lpext(number_fringe_nodes) % lnods(1,nboun_local + ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(2,nboun_local + ielem_aux)==0 .or. &
                   lpext(number_fringe_nodes) % lnods(3,nboun_local + ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(4,nboun_local + ielem_aux)==0  )print*,'BOUCEAUUUUUUUXXXXXXXXXXX',ielem_aux,ipoin,tetex
!!$
!!$              do inode = 1,nnode(lpext(number_fringe_nodes) % ltype(nboun_local + ielem_aux))
!!$                 ipoex = lpext(number_fringe_nodes) % lnods(inode,nboun_local + ielem_aux)
!!$                 do idime=1,ndime
!!$                    coorcg(idime) = coorcg(idime) + coord(idime,ipoex)
!!$                 end do
!!$              end do
!!$              call elsest(&
!!$                   2_ip,jsubd,ielse,mnode,ndime,neighbor_subdomain % npoin,neighbor_subdomain % nelem,&
!!$                   nnode(1:),neighbor_subdomain % lnods,neighbor_subdomain % ltype,&
!!$                   ltopo,neighbor_subdomain % coord,coorcg,relse,jelem,&
!!$                   shape_dod,deriv_dod,xieta,dumm1)
!!$              lpext(number_fringe_nodes) % lmate(nboun_local + ielem_aux) = lmate(jelem)

              poex1 = lext2(1)  
              poex2 = lext2(kboun)
              if (qboun==kboun  .or. lext2(1)/=lext2(nboun_local)) then
                 poex3 = lnodb_global(lpert(1,jview),list_ordered_boundaries(qboun))
                 if( ltypb_global(list_ordered_boundaries(qboun)) == 12)poex3 = lnodb_global(lperm(2,jview),list_ordered_boundaries(qboun))
              else 
                 poex3 = lnodb_global(lpert(2,jview),list_ordered_boundaries(qboun))
                 if( ltypb_global(list_ordered_boundaries(qboun)) == 12)poex3 = lnodb_global(lperm(3,jview),list_ordered_boundaries(qboun))
              end if
              ielem_aux = ielem_aux + 1

              imatn          = lmatn_dod(poex1)
              ipoiz          = lpoiz_dod(poex1)
              lpext(number_fringe_nodes) % ltype(nboun_local + ielem_aux) = TET04
              lpext(number_fringe_nodes) % lmate(nboun_local + ielem_aux) = imatn
              lpext(number_fringe_nodes) % lelez(nboun_local + ielem_aux) = ipoiz
              lpext(number_fringe_nodes) % lesub(nboun_local + ielem_aux) = isubd
              call dod_extens_qual3d(3_ip,ipoin,poex1,poex2,poex3,kappa,asrad,q,sign,dumm5)
              lpext(number_fringe_nodes) % lnods(1,nboun_local + ielem_aux) = ipoin
              lpext(number_fringe_nodes) % lnods(2,nboun_local + ielem_aux) = poex1 
              lpext(number_fringe_nodes) % lnods(3,nboun_local + ielem_aux) = poex2
              lpext(number_fringe_nodes) % lnods(4,nboun_local + ielem_aux) = poex3
              nelem_dod = nelem_dod + 1 

              if(lpext(number_fringe_nodes) % lnods(1,nboun_local + ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(2,nboun_local + ielem_aux)==0 .or. &
                   lpext(number_fringe_nodes) % lnods(3,nboun_local + ielem_aux)==0 .or.lpext(number_fringe_nodes) % lnods(4,nboun_local + ielem_aux)==0  )print*,'BOUCEAUUUUUUUXXXXXXXXXXX',ielem_aux,ipoin,tetex
!!$              do inode = 1,nnode(lpext(number_fringe_nodes) % ltype(nboun_local + ielem_aux))
!!$                 ipoex = lpext(number_fringe_nodes) % lnods(inode,nboun_local + ielem_aux)
!!$                 do idime=1,ndime
!!$                    coorcg(idime) = coorcg(idime) + coord(idime,ipoex)
!!$                 end do
!!$              end do
!!$              call elsest(&
!!$                   2_ip,jsubd,ielse,mnode,ndime,neighbor_subdomain % npoin,neighbor_subdomain % nelem,&
!!$                   nnode(1:),neighbor_subdomain % lnods,neighbor_subdomain % ltype,&
!!$                   ltopo,neighbor_subdomain % coord,coorcg,relse,jelem,&
!!$                   shape_dod,deriv_dod,xieta,dumm1)
!!$              lpext(number_fringe_nodes) % lmate(nboun_local + ielem_aux) = lmate(jelem)

              if(ipoin==10023 .or. ipoin==2468  .or. ipoin==25748)print*,'TETEXXXX',tetex,ielem_aux
           else if(tetex>2 .or. ielem_aux>2 )then
              ielem_aux = 0

              do kboun=1,nboun_local
                 eboun = list_ordered_boundaries(kboun)
                 call dod_order_boundaries(&
                      2_ip,ipoin,mnodb,nboun_local,npoin_global,nedge_ipoin,ltypb_global,&
                      lnodb_global,lboch_global,ledbo_global,ledge_npoin_global(ipoin)%l,&
                      list_ordered_boundaries,lperm,iview,kboun)

                 if(kboun>1 ) then 
                    if(lext2(kboun)/=lext2(kboun-1))then
                       ielem_aux= ielem_aux + 1

                       poex1 = lext2(kboun)
                       poex2 = lext2(kboun-1)
                       poex3 = lnodb_global(lperm(2,iview),eboun)
                       if( ltypb_global(eboun) == 12)poex3 = lnodb_global(lperm(3,iview),eboun)
                       imatn          = lmatn_dod(poex1)
                       ipoiz          = lpoiz_dod(poex1)
                       lpext(number_fringe_nodes) % ltype(nboun_local + ielem_aux) = TET04
                       lpext(number_fringe_nodes) % lmate(nboun_local + ielem_aux) = imatn
                       lpext(number_fringe_nodes) % lelez(nboun_local + ielem_aux) = ipoiz
                       lpext(number_fringe_nodes) % lesub(nboun_local + ielem_aux) = isubd
                       call dod_extens_qual3d(3_ip,ipoin,poex1,poex2,poex3,kappa,asrad,q,sign,dumm5)
                       lpext(number_fringe_nodes) % lnods(1,nboun_local + ielem_aux) = ipoin
                       lpext(number_fringe_nodes) % lnods(2,nboun_local + ielem_aux) = poex1 
                       lpext(number_fringe_nodes) % lnods(3,nboun_local + ielem_aux) = poex2
                       lpext(number_fringe_nodes) % lnods(4,nboun_local + ielem_aux) = poex3
                       nelem_dod = nelem_dod + 1 

                       if(  lpext(number_fringe_nodes) % lnods(1,nboun_local + ielem_aux)==0 .or. &
                            lpext(number_fringe_nodes) % lnods(2,nboun_local + ielem_aux)==0 .or. &
                            lpext(number_fringe_nodes) % lnods(3,nboun_local + ielem_aux)==0 .or. &
                            lpext(number_fringe_nodes) % lnods(4,nboun_local + ielem_aux)==0 )then
                          print*,ipoin,'cerrrrrrrrrrrrrrrrrrro-b'
                          stop
                       end if
!!$                       do inode = 1,nnode(lpext(number_fringe_nodes) % ltype(nboun_local + ielem_aux))
!!$                          ipoex = lpext(number_fringe_nodes) % lnods(inode,nboun_local + ielem_aux)
!!$                          do idime=1,ndime
!!$                             coorcg(idime) = coorcg(idime) + coord(idime,ipoex)
!!$                          end do
!!$                       end do
!!$                       call elsest(&
!!$                            2_ip,jsubd,ielse,mnode,ndime,neighbor_subdomain % npoin,neighbor_subdomain % nelem,&
!!$                            nnode(1:),neighbor_subdomain % lnods,neighbor_subdomain % ltype,&
!!$                            ltopo,neighbor_subdomain % coord,coorcg,relse,jelem,&
!!$                            shape_dod,deriv_dod,xieta,dumm1)
!!$                       lpext(number_fringe_nodes) % lmate(nboun_local + ielem_aux) = lmate(jelem)

                    end if
                 end if
                 if(kboun==nboun_local)then
                    if(lext2(kboun)/=lext2(1))then
                       ielem_aux              = ielem_aux + 1
                       poex1 = lnodb_global(lperm(1,iview),eboun)
                       if( ltypb_global(eboun) == 12)poex1 = lnodb_global(lperm(2,iview),eboun)
                       poex2 = lext2(1)    
                       poex3 = lext2(kboun)  
                       imatn          = lmatn_dod(poex2)
                       ipoiz          = lpoiz_dod(poex2)
                       lpext(number_fringe_nodes) % ltype(nboun_local + ielem_aux) = TET04
                       lpext(number_fringe_nodes) % lmate(nboun_local + ielem_aux) = imatn
                       lpext(number_fringe_nodes) % lelez(nboun_local + ielem_aux) = ipoiz
                       lpext(number_fringe_nodes) % lesub(nboun_local + ielem_aux) = isubd
                       call dod_extens_qual3d(3_ip,ipoin,poex1,poex2,poex3,kappa,asrad,q,sign,dumm5)
                       lpext(number_fringe_nodes) % lnods(1,nboun_local + ielem_aux) = ipoin
                       lpext(number_fringe_nodes) % lnods(2,nboun_local + ielem_aux) = poex1 
                       lpext(number_fringe_nodes) % lnods(3,nboun_local + ielem_aux) = poex2
                       lpext(number_fringe_nodes) % lnods(4,nboun_local + ielem_aux) = poex3
                       nelem_dod = nelem_dod + 1 

                       if(  lpext(number_fringe_nodes) % lnods(1,nboun_local + ielem_aux)==0 .or. &
                            lpext(number_fringe_nodes) % lnods(2,nboun_local + ielem_aux)==0 .or. &
                            lpext(number_fringe_nodes) % lnods(3,nboun_local + ielem_aux)==0 .or. &
                            lpext(number_fringe_nodes) % lnods(4,nboun_local + ielem_aux)==0 )then
                          print*,ipoin,'cerrrrrrrrrrrrrrrrrrro-b'
                          stop
                       end if
!!$                       do inode = 1,nnode(lpext(number_fringe_nodes) % ltype(nboun_local + ielem_aux))
!!$                          ipoex = lpext(number_fringe_nodes) % lnods(inode,nboun_local + ielem_aux)
!!$                          do idime=1,ndime
!!$                             coorcg(idime) = coorcg(idime) + coord(idime,ipoex)
!!$                          end do
!!$                       end do
!!$                       call elsest(&
!!$                            2_ip,jsubd,ielse,mnode,ndime,neighbor_subdomain % npoin,neighbor_subdomain % nelem,&
!!$                            nnode(1:),neighbor_subdomain % lnods,neighbor_subdomain % ltype,&
!!$                            ltopo,neighbor_subdomain % coord,coorcg,relse,jelem,&
!!$                            shape_dod,deriv_dod,xieta,dumm1)
!!$                       lpext(number_fringe_nodes) % lmate(nboun_local + ielem_aux) = lmate(jelem)
!!$
                    end if
                 end if
              end do
              if(nutri>0)then
                 print*,'TENGO QUE AÑADIR EL DOD_CATALA...LA VIDA ES ASI',ipoin,nutri,ielem_aux
                 !stop
                 call dod_concat(ipoin,isubd,nboun_local,lext2,nboun_local,ielem_aux)
              end if
!!!a los de concat no les he modificado lo del lmate!!!!!!!
              if(  lpext(number_fringe_nodes) % lnods(1,nboun_local+ielem_aux)==0 .or. &
                   lpext(number_fringe_nodes) % lnods(2,nboun_local+ielem_aux)==0 .or. &
                   lpext(number_fringe_nodes) % lnods(3,nboun_local+ielem_aux)==0 .or. &
                   lpext(number_fringe_nodes) % lnods(4,nboun_local+ielem_aux)==0 )then
                 print*,ipoin,'cerrrrrrrrrrrrrrrrrrro-b'
                 stop
              end if

           end if
           deallocate(indic,stat=istat)
           deallocate(q_max,stat=istat)
           deallocate(quali,stat=istat)
           deallocate(lpoex,stat=istat)
           deallocate(lextp,stat=istat)
           deallocate(lcont,stat=istat)
        end if
!!!OJOOO end if
!if(ipoin==25748)stop
     end do punto_contorno
     call memory_deallo(mem_servi(1:2,servi),'NORMZ','dod_extension_3d',normz)
     call memory_deallo(mem_servi(1:2,servi),'NORMY','dod_extension_3d',normy)
     call memory_deallo(mem_servi(1:2,servi),'NORMX','dod_extension_3d',normx)
     call memory_deallo(mem_servi(1:2,servi),'PARA','dod_extension_3d',para)
     call memory_deallo(mem_servi(1:2,servi),'LEXT2' ,'dod_extension_3d',lext2)           
     call memory_deallo(mem_servi(1:2,servi),'LIST_ORDERED_BOUNDARIES','dod_extension_3d',list_ordered_boundaries)
     call memory_deallo(mem_servi(1:2,servi),'DISTEX',   'dod_extension_3d',distex)  
     call memory_deallo(mem_servi(1:2,servi),'HCARAC','dod_extension_3d',hcarac)
     call memory_deallo(mem_servi(1:2,servi),'COORCG1','dod_extension_3d',coorcg1)

  end do
  !stop
  call memory_deallo(mem_servi(1:2,servi),'LNEXT' ,'dod_extension_3d',lnext)

  call memory_deallo(mem_servi(1:2,servi),'LINVP_NBOUN','dod_memall',linvp_nboun)
end subroutine dod_extension_3d

