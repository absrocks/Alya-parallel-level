subroutine chm_elmusr()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_elmusr
  ! NAME 
  !    chm_elmusr
  ! DESCRIPTION
  !    This routine 
  ! USES
  !    chm_elmgat
  !    elmder
  ! USED BY
  !    chm_outset
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_chemic
  implicit none
  integer(ip), save :: ipass=0
  integer(ip)       :: pnode,pgaus,pelty,ispec,ispev,nvalu,iodei,iodes
  integer(ip)       :: ielem,igaus,idime,inode,ipoin,ivalu,iodev
  real(rp)          :: xjacm(ndime,ndime) 
  real(rp)          :: elcod(ndime,mnode)
  real(rp)          :: elcon(mnode,nspec_chm)
  real(rp)          :: gpcon(mgaus,nspec_chm)
  real(rp)          :: gpvol(mgaus),gpdet,T,rweig

  if( kfl_model_chm == 1 ) then
     !
     ! Initialization
     !
     if (wprob_chm == 'COIN1' ) then
        nvalu = 3
     else if (wprob_chm == 'INVN1' ) then
        nvalu = 6
     else if (wprob_chm == 'INTE1' ) then
        nvalu = 4
        iodei = nodes_chm / 2
        iodev = iodei + 1
     else 
        return
     end if

     do ivalu = 1,nvalu
        xvael_chm(ivalu) = 0.0_rp
     end do

     if( INOTMASTER ) then
        !
        ! Loop over elements
        !
        elements: do ielem = 1,nelem
           ! 
           ! Element properties and dimensions
           !
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           !
           ! Gather operations
           !
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              do ispec = 1,nspec_chm
                 elcon(inode,ispec) = conce(ipoin,ispec,1)
              end do
              do idime = 1,ndime
                 elcod(idime,inode) = coord(idime,ipoin)
              end do
           end do
           !
           ! Gauss point values
           !
           do igaus = 1,pgaus
              do ispec = 1,nspec_chm           
                 gpcon(igaus,ispec) = 0.0_rp
              end do
           end do
           do ispec = 1,nspec_chm           
              do igaus = 1,pgaus
                 do inode = 1,pnode
                    gpcon(igaus,ispec) = gpcon(igaus,ispec) &
                         +elmar(pelty)%shape(inode,igaus)   &
                         *elcon(inode,ispec)
                 end do
              end do
           end do
           !
           ! dV:=GPVOL=|J|*wg
           !
           do igaus = 1,pgaus    
              call jacdet(&
                   ndime,pnode,elcod,elmar(pelty)%deriv(1,1,igaus),&
                   xjacm,gpdet)
              gpvol(igaus) = elmar(pelty)%weigp(igaus)*gpdet 
           end do
           !
           ! Output
           !
           if( wprob_chm == 'COIN1' ) then
              do igaus = 1,pgaus    
                 xvael_chm(1) = xvael_chm(1) + gpvol(igaus) * gpcon(igaus,1)                         ! I
                 xvael_chm(2) = xvael_chm(2) + gpvol(igaus) * gpcon(igaus,2)                         ! I2
                 do ispec = 2,nspec_chm
                    xvael_chm(3) = xvael_chm(3) + gpvol(igaus) * real(ispec) * gpcon(igaus,ispec)    ! Sum^n_i=2 n*In 
                 end do
              end do
           else if( wprob_chm == 'INVN1' ) then
              do igaus = 1,pgaus    
                 xvael_chm(1) = xvael_chm(1) + gpvol(igaus) * gpcon(igaus,1)                         ! I
                 xvael_chm(2) = xvael_chm(2) + gpvol(igaus) * gpcon(igaus,2)                         ! V
                 xvael_chm(3) = xvael_chm(3) + gpvol(igaus) * gpcon(igaus,nclas_chm+1)               ! I2
                 xvael_chm(4) = xvael_chm(4) + gpvol(igaus) * gpcon(igaus,nclas_chm+nodes_chm/2+1)   ! V2
                 ispev = 2+nodes_chm/2
                 do ispec = 3,2+nodes_chm/2
                    ispev    = ispev + 1
                    xvael_chm(5) = xvael_chm(5) + gpvol(igaus) * real(ispec-1) * gpcon(igaus,ispec)  ! Sum^n_i=2 n*In 
                    xvael_chm(6) = xvael_chm(6) + gpvol(igaus) * real(ispec-1) * gpcon(igaus,ispev)  ! Sum^n_i=2 n*Vn 
                 end do
              end do
           else if( wprob_chm == 'INTE1' ) then
              do igaus = 1,pgaus    
                 xvael_chm(1) = xvael_chm(1) + gpvol(igaus) * gpcon(igaus,1)                         ! I
                 xvael_chm(2) = xvael_chm(2) + gpvol(igaus) * gpcon(igaus,2)                         ! V
                 xvael_chm(3) = xvael_chm(3) + gpvol(igaus) * gpcon(igaus,3)                         ! I2
                 xvael_chm(4) = xvael_chm(4) + gpvol(igaus) * &
                      &                (gpcon(igaus,1)+gpcon(igaus,2)+2.0_rp*gpcon(igaus,3)) ! I + V + I2
                 rweig = 2.0_rp
                 ispec = nclas_chm
                 do iodes = 1,iodei
                    ispec = ispec + 1
                    rweig = rweig + 1.0_rp
                    xvael_chm(4) = xvael_chm(4) + gpvol(igaus) * rweig * gpcon(igaus,ispec)          ! + Sum^n_i=3 n*In 
                 end do

                 rweig = 1.0_rp
                 ispec = nclas_chm + iodei
                 do iodes = 1,iodev
                    ispec = ispec + 1
                    rweig = rweig + 1.0_rp
                    xvael_chm(4) = xvael_chm(4) + gpvol(igaus) * rweig * gpcon(igaus,ispec)          ! + Sum^n_i=3 n*In 
                 end do

              end do
           end if

        end do elements

     end if
     !
     ! Sum over subdomains
     !
     call pararr('SUM',0_ip,nvalu,xvael_chm)
     !
     ! Output
     !
     if( INOTSLAVE ) then
        if( ipass == 0 ) then
           if( wprob_chm == 'COIN1' ) then
              write(lun_resu1_chm,4) 
           else if( wprob_chm == 'INVN1' ) then
              write(lun_resu1_chm,1) 
           else if( wprob_chm == 'INTE1' ) then
              write(lun_resu1_chm,3) 
           end if
           ipass = 1
        end if
        if( wprob_chm == 'INTE1' ) then
           call chm_usrtem(T)
           write(lun_resu1_chm,2) cutim,T,(xvael_chm(ivalu),ivalu=1,nvalu)
        else
           write(lun_resu1_chm,2) cutim,(xvael_chm(ivalu),ivalu=1,nvalu)
        end if
     end if

  end if

4 format('         Time','            I','           I2','          nIn')
1 format('         Time','            I','            V','           I2','           V2','          nIn','          nVn')
3 format('         Time','  Temperature','            I','            V','           I2','      nIn+nVn')
2 format(101(1x,e12.6))

end subroutine chm_elmusr
