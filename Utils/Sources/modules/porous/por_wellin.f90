!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_wellin.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Obtains the Well Index
!> @details Obtains the Well transmisibilities(Toi) and then the Well Index
!> @} 
!------------------------------------------------------------------------

subroutine por_wellin()
  use def_kintyp, only :  ip,rp
  use def_parame, only :  pi
  use def_domain, only :  npoin,ltype,lnnod,ngaus,llapl
  use def_domain, only :  lorde,ltopo,lnods,elmar,ndime,vmass,mnode
  use def_domain, only :  ntens,mgaus,coord,r_dom,c_dom,nelem
  use def_porous, only :  perme_por,winde_por,iwell_por,tywel_por
  use def_porous, only :  nwell_por,wmass_por
  use def_master, only :  ID_POROUS,mem_modul,modul,INOTMASTER
  use def_master, only :  gisca,pard1,pard2,pard3,NPOIN_TYPE,lninv_loc
  use def_master, only :  current_zone
  use mod_memory, only :  memory_alloca,memory_deallo
  use mod_matrix, only :  matrix_assrhs
  implicit none
  integer(ip)          :: pelty,pnode,kelem,inode,idime,jnode,jpoin
  integer(ip)          :: pgaus,plapl,porde,ptopo,izdom,kpoin
  integer(ip)          :: ielem,igaus,ipoin,iwell,jwell
  real(rp)             :: permx,smope,dista,distz,auxii,re,h3
  real(rp)             :: radiu
  real(rp)             :: elrhs(3_ip,mnode)
  real(rp)             :: elcod(ndime,mnode)
  real(rp)             :: gphes(ntens,mnode,mgaus)   
  real(rp)             :: gpvol(mgaus)               
  real(rp)             :: gpcar(ndime,mnode,mgaus)  
  real(rp),   pointer  :: rhsid(:,:)
  real(rp),   pointer  :: tmass(:)

  if( nwell_por > 0 ) then

     allocate( tmass(nwell_por) )

     if( INOTMASTER ) then

        nullify(rhsid)
        call memory_alloca(mem_modul(1:2,modul),'RHSID','por_wellin',rhsid,3_ip,npoin)

        elements: do ielem = 1,nelem
           !
           ! Element dimensions
           !
           pelty = ltype(ielem)
           if( pelty > 0 ) then
              pnode = lnnod(ielem)
              pgaus = ngaus(pelty)
              plapl = llapl(pelty)
              porde = lorde(pelty)
              ptopo = ltopo(pelty)
              !
              ! Gather operations
              !
              elcod(1:ndime,1:pnode) = coord(1:ndime,lnods(1:pnode,ielem))
              !
              ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
              !
              gphes = 0.0_rp
              call elmcar(&
                   pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
                   elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                   gphes,ielem)

              permx = perme_por(1,ielem)
              elrhs = 0.0_rp

              do jnode = 1,pnode
                 do igaus = 1,pgaus
                    elrhs(3,jnode) = elrhs(3,jnode) + elmar(pelty) % shape(jnode,igaus) * permx * gpvol(igaus)   ! Smooth perm_x
                 end do
                 jpoin = lnods(jnode,ielem)
                 if( iwell_por(jpoin) /= 0 ) then
                    !
                    ! 1) sum (Toi), 2) sum (Toi * ln(ri)) , 3) Smooth perm_x
                    !              
                    do inode = 1,pnode
                       if( inode /= jnode ) then
                          dista = sqrt( ( (elcod(1,inode)-elcod(1,jnode))**2_ip ) &
                               &      + ( (elcod(2,inode)-elcod(2,jnode))**2_ip ) )
                          distz = sqrt( (  elcod(3,inode)-elcod(3,jnode))**2_ip )
                          if ( dista > distz ) then                                                              ! Eliminate nodes directly on top and below 
                             do igaus = 1,pgaus
                                do idime = 1,2_ip                                                                ! only x and y
                                   auxii = gpcar(idime,jnode,igaus)*gpcar(idime,inode,igaus)*gpvol(igaus)*permx
                                   elrhs(1,jnode) = elrhs(1,jnode) + auxii                                       ! I directly sum from nodes conected to jnode
                                   elrhs(2,jnode) = elrhs(2,jnode) + ( auxii * log (dista) )   
                                end do
                             end do
                          end if
                       end if
                    end do
                 end if
              end do
           end if

           call matrix_assrhs(3_ip,1_ip,pnode,npoin,lnods(1:pnode,ielem),elrhs,rhsid)

        end do elements
        !
        ! Periodicity and Parall service
        !
        call rhsmod(3_ip,rhsid)
        !
        ! Loop over edges IPOIN-JPOIN
        ! Check which other subdomains owns it. The one with lowest rank keeps it
        !
        call memgen(1_ip,npoin,0_ip)
        do ipoin = 1,npoin
           iwell = iwell_por(ipoin)
           if( iwell /= 0 ) then
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 jwell = iwell_por(jpoin)         
                 if( iwell == jwell .and. ipoin < jpoin ) then
                    !
                    ! Check if I should keep this edge
                    !                         
                    pard1 = ipoin
                    pard2 = jpoin
                    pard3 = 1
                    call Parall(431_ip)
                    if( pard3 == 1 ) then
                       !
                       ! I own this edge
                       !
                       gisca(ipoin) = iwell
                       gisca(jpoin) = jwell
                    end if
                 end if
              end do
           end if
        end do
        !
        ! WINDE_POR: Compute edge lengths
        !
        do ipoin = 1,npoin
           winde_por(ipoin) = 0.0_rp        
        end do
        do ipoin = 1,npoin
           iwell = gisca(ipoin)
           if( iwell /= 0 ) then
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 jwell = gisca(jpoin)              
                 if( iwell == jwell .and. lninv_loc(jpoin) < lninv_loc(ipoin) ) then
                    dista = 0.0_rp
                    do idime = ndime,ndime
                       dista = dista + ( coord(idime,ipoin) - coord(idime,jpoin) ) ** 2
                    end do
                    dista            = sqrt(dista)
                    winde_por(ipoin) = winde_por(ipoin) + dista
                    winde_por(jpoin) = winde_por(jpoin) + dista
                 end if
              end do
           end if
        end do
        !
        ! Compute total mass for each well
        !
        do iwell = 1,nwell_por
           tmass(iwell) = 0.0_rp
        end do
        do ipoin = 1,npoin
           iwell = gisca(ipoin)
           if( iwell /= 0 ) &
                tmass(iwell) = tmass(iwell) + winde_por(ipoin)
        end do

     end if
     !
     ! Total mass
     !
     call pararr('SUM',0_ip,nwell_por,tmass)     
     
     if( INOTMASTER ) then
        !
        ! Mass matrix
        !
        do ipoin = 1,npoin
           iwell = gisca(ipoin)
           if( iwell > 0 ) &
                wmass_por(ipoin) = winde_por(ipoin) / tmass(iwell) 
        end do         
        !
        ! Parallel exchange
        !
        call pararr('SLX',NPOIN_TYPE,npoin,winde_por)
        !
        ! Compute well index
        !
        do ipoin = 1,npoin
           iwell = iwell_por(ipoin)
           if( iwell > 0 ) then
              h3               = 0.5_rp * winde_por(ipoin)
              smope            = rhsid(3,ipoin) / vmass(ipoin)
              re               = exp( (rhsid(2,ipoin) - (2.0_rp*pi*smope*h3) ) / rhsid(1,ipoin) )
              radiu            = tywel_por(iwell) % radiu
              winde_por(ipoin) = 2.0_rp * pi * smope * h3 / log ( re / radiu )
           end if
        end do

        call memgen(3_ip,npoin,0_ip)
        call memory_deallo(mem_modul(1:2,modul),'RHSID','por_wellin',rhsid )

     end if

     deallocate( tmass )

  end if

end subroutine por_wellin
