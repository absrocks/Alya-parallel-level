subroutine ibm_updtss()
  !-----------------------------------------------------------------------
  !****f* Nastin/ibm_updtss
  ! NAME 
  !    ibm_updtss
  ! DESCRIPTION
  !    This routine computes the time step size for the incompressible NS
  !    equation. 
  ! USED BY
  !    ibm_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_immbou
  implicit none 
  integer(ip) :: iimbo,iboun,pblty,pnodb,pgaub,inodb,idime,igaub
  integer(ip) :: kgaus,inode,ipoin,pnode,jelem,pelty
  real(rp)    :: xfact,eucta,baloc(9),bocod(ndime,mnoib)
  real(rp)    :: gbpos(3),x(3),v(3),h,xjacm(9),vnorm,gpdet
  real(rp)    :: elcod(ndime,mnode)

  !
  ! OJO: Quitar esto
  !
  return

  if( kfl_timei_ibm /= 0 .and. kfl_coibm /= 0 ) then

     xfact     = 1.0_rp/real(ndime)
     x(3)      = 0.0_rp
     v(3)      = 0.0_rp
     dtcri_ibm = 1.0e9_rp
     kgaus     = 0

     if( INOTMASTER ) then
        do iimbo = 1,nimbo

           do iboun = 1,imbou(iimbo)%nboib
              !
              ! Boundary properties and dimensions
              !
              pblty    = imbou(iimbo)%ltyib(iboun)
              pnodb    = nnode(pblty)
              pgaub    = ngaib(pblty)
              !
              ! BOCOD: Gather operations
              !
              do inodb = 1,pnodb
                 ipoin = imbou(iimbo)%lnoib(inodb,iboun)
                 do idime = 1,ndime
                    bocod(idime,inodb) = imbou(iimbo)%cooib(idime,ipoin)
                 end do
              end do

              do igaub = 1,pgaub

                 call bouder(&
                      pnodb,ndime,ndimb,elmar(pblty)%derib(1,1,igaub),& 
                      bocod,baloc,eucta) 
                 kgaus = kgaus + 1
                 !
                 ! Gauss point coordinates
                 !
                 do idime = 1,ndime
                    gbpos(idime) = 0.0_rp
                    do inodb = 1,pnodb
                       gbpos(idime) = gbpos(idime) &
                            + elmar(pblty)%shaib(inodb,igaub) * bocod(idime,inodb)
                    end do
                 end do

                 x(1)     = gbpos(1) - imbou(iimbo) % posil(    1,1)
                 x(2)     = gbpos(2) - imbou(iimbo) % posil(    2,1)
                 x(ndime) = gbpos(ndime) - imbou(iimbo) % posil(ndime,1)
                 v(1)     = 0.0_rp
                 v(2)     = 0.0_rp
                 call vecpro(imbou(iimbo)%veloa,x,v,3_ip)            
                 v(1)     = v(1)     + imbou(iimbo) % velol(    1,1)
                 v(2)     = v(2)     + imbou(iimbo) % velol(    2,1)
                 v(ndime) = v(ndime) + imbou(iimbo) % velol(ndime,1)           
                 vnorm    = sqrt( v(1)*v(1) + v(2)*v(2) + v(3)*v(3) ) + epsilon(1.0_rp)
                 !
                 ! Element properties
                 !             
                 jelem    = lgaib(kgaus)
                 if( jelem /= 0 ) then
                    pelty    = ltype(jelem)
                    if( pelty > 0 ) then
                       pnode = nnode(pelty)
                       do inode = 1,pnode
                          ipoin = lnods(inode,jelem)
                          do idime = 1,ndime
                             elcod(idime,inode) = coord(idime,ipoin)
                          end do
                       end do
                       call jacdet(&
                            ndime,pnode,elcod,elmar(pelty)%dercg,&
                            xjacm,gpdet)
                       h = ( elmar(pelty) % weicg * gpdet ) ** xfact
                       dtcri_ibm = min(dtcri_ibm,h/vnorm)
                    end if
                 end if
              end do
           end do
        end do
     end if
     !
     ! Parall: Look for minimum over all subdomains (dtmin)
     !
     call pararr('MIN',0_ip,1_ip,dtcri_ibm)
     !
     ! Assign 1/dt
     !
     if( dtcri_ibm /= 0.0_rp ) dtinv_ibm = 1.0_rp / ( dtcri_ibm * safet_ibm )
     if( kfl_timco     == 1 )  dtinv     = max(dtinv,dtinv_ibm)

  end if

end subroutine ibm_updtss
