subroutine chm_bouset(ibsec,ibset)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_bouset
  ! NAME 
  !    chm_bouset
  ! DESCRIPTION
  !    This routine computes variables on a boundary set W.
  !    The variable are:
  !    1. setfl: set flux of class i  =  int_S k grad(Ci).n ds
  ! USES
  !    bouder
  !    chenor
  ! USED BY
  !    chm_outset
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_ker_proper

  implicit none
  integer(ip), intent(in)  :: ibsec,ibset
  real(rp),    pointer     :: setsu(:)
  real(rp),    pointer     :: set_conce(:)
  real(rp),    pointer     :: set_mass_flux(:)
  real(rp),    pointer     :: set_avgco(:)
  real(rp)                 :: baloc(ndime,ndime)
  real(rp)                 :: bocod(ndime,mnodb)
  real(rp)                 :: bocon(nspec_chm, mnodb)
  real(rp)                 :: elcod(ndime,mnode)
  real(rp)                 :: gbcon(nspec_chm, mnodb)
  real(rp)                 :: bovel(ndime,mnodb)
  real(rp)                 :: gbvel(ndime,mgaus)
  integer(ip)              :: ielem,inode,ipoin,nn
  integer(ip)              :: igaus,idime,igaub,iboun,inodb,pblty
  integer(ip)              :: kboun,ispec,dummi
  integer(ip)              :: pnodb,pmate,pnode,pelty,pgaus,pgaub
  real(rp)                 :: eucta,dsurf,xfact
  real(rp)                 :: gbden(mgaus),mcon
  real(rp)                 :: dummr(ndime,mnode)

  if( INOTMASTER ) then

     !----------------------------------------------------------------------
     !
     ! Initialization
     !
     !----------------------------------------------------------------------

     nn            =  postp(1) % nvabs + 1
     setsu         => vbset( nn:nn , ibset )             ! Surface           
     set_mass_flux => vbset(  1: 1 , ibset )             ! Mass flux 
     set_conce     => vbset(  2: 9 , ibset )             ! Mean species mass fraction     
!     set_avgco     => vbset( 10:18 , ibset )             ! Mean mass fraction

     setsu         = 0.0_rp
     set_conce     = 0.0_rp
     set_mass_flux = 0.0_rp
!     set_avgco     = 0.0_rp
     !
     ! Loop over elements
     !
     boundaries: do iboun = 1,nboun

        if( lbset(iboun) == ibsec ) then

           !----------------------------------------------------------------
           !
           ! Element properties and dimensions and gather
           !
           !----------------------------------------------------------------

           pblty = ltypb(iboun)
           pnodb = nnode(pblty)
           pgaub = ngaus(pblty)
           pmate = 1
           
           do inodb = 1,pnodb
              ipoin = lnodb(inodb,iboun)
              do ispec = 1,nspec_chm
                 bocon(ispec,inodb) = conce(ipoin,ispec,1)
              enddo
              do idime = 1,ndime
                 bocod(idime,inodb) = coord(idime,ipoin)
                 bovel(idime,inodb) = veloc(idime,ipoin,1) !!!DEFINE
              end do
           end do
          
           ielem = lelbo(iboun)
           pelty = ltype(ielem)
           if( nmate > 1 ) pmate = lmate(ielem)
           
           if (pelty > 0) then
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 do idime = 1,ndime
                    elcod(idime,inode) = coord(idime,ipoin)
                 end do
              end do
              
              gbcon = 0.0_rp
              gbvel = 0.0_rp
              do igaub = 1,pgaub
                 do inodb = 1,pnodb
                    do ispec = 1,nspec_chm
                       gbcon(ispec, igaub) = gbcon(ispec,igaub) + elmar(pblty)%shape(inodb,igaub) * bocon(ispec, inodb)
                    enddo
                    do idime = 1,ndime
                       gbvel(idime,igaub) = gbvel(idime,igaub) + elmar(pblty)%shape(inodb,igaub) * bovel(idime,inodb)
                    end do
                 end do
              end do
              
              if( kfl_prope /= 0 ) then
                 call ker_proper('DENSI','PGAUB',dummi,iboun,gbden)
              else
                 call runend('CHEMIC GOT THIS FAR WITHOUT KERMOD')
              endif

              !----------------------------------------------------------------
              !
              ! Loop over Gauss points
              !
              !----------------------------------------------------------------           

              gauss_points: do igaub = 1,pgaub        

                 call bouder(&
                      pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&
                      bocod,baloc,eucta)
                 call chenor(pnode,baloc,bocod,elcod)
                 dsurf = elmar(pblty)%weigp(igaub)*eucta 
                 setsu = setsu + dsurf
                 
                 !-------------------------------------------------------------
                 !
                 ! Mass flux  
                 !
                 !-------------------------------------------------------------

                 if( postp(1) % npp_setsb(1) /= 0 ) then
                   do idime=1,ndime
                      set_mass_flux = set_mass_flux + gbden(igaub) * gbvel(idime,igaub) * baloc(idime,ndime) * dsurf
                   end do
                 end if

                 !-------------------------------------------------------------
                 !
                 ! Species mass fraction 
                 !
                 !-------------------------------------------------------------
 
                 if( postp(1) % npp_setsb(2) /= 0 ) then
                    do ispec = 1, nspec_chm
                       do idime=1,ndime
                          mcon = gbden(igaub)*gbcon(ispec,igaub)
                          set_conce(ispec) = set_conce(ispec) + mcon * gbvel(idime,igaub)*baloc(idime,ndime) * dsurf
                       end do
                    enddo
                 end if
                                  
              end do gauss_points
           end if
        end if

     end do boundaries
 
  end if 

end subroutine chm_bouset
