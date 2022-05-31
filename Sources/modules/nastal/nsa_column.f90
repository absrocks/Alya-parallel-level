!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_column.f90
!> @author  Simone Marras
!> @date    16/11/1966
!> @brief   
!> @details 
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_column()
  use def_kintyp
  use def_master
  use def_domain
  implicit none
  type type_column
     integer(ip), pointer :: lnode(:)
     integer(ip), pointer :: lposi(:)
     integer(ip)          :: nnode
     integer(ip)          :: iglob
     integer(ip)          :: ilast
     real(rp)             :: s
  end type type_column
  integer(ip)                :: icolo_glo,icolo_loc,ipoin,ncolo_loc,inode,pnode
  integer(ip)                :: icolo_max,comcont,jnext,ntota_glo,itota_glo
  integer(ip)                :: inext,itota,ntota,knext,lnext,kposi,iposi_max
  integer(ip)                :: itota_old
  integer(ip), pointer       :: itest(:)   
  integer(ip), pointer       :: my_paris(:)   
  real(rp),    pointer       :: my_parrs(:)
  integer(ip), pointer       :: my_parin(:)   
  real(rp),    pointer       :: my_parre(:)
  type(type_column), pointer :: colum(:)
  !
  ! Fill in type
  !

  if (nzone > 1) call runend("NSA_INITIAL_CONDITIONS: THIS SUB IS NOT PREPARED TO RUN WITH ZONES.")

  if( INOTMASTER ) then
     call memgen(1_ip,npoin,0_ip)
     iposi_max = 0
     icolo_max = 0
     do ipoin = 1,npoin
        icolo_glo = int(xfiel(1) % a(1,ipoin,1),ip)
        icolo_max = max(icolo_max,icolo_glo)
        iposi_max = max(iposi_max,int(xfiel(1) % a(2,ipoin,1),ip))
        gisca(icolo_glo) = gisca(icolo_glo) + 1
     end do
  end if
  call parari('MAX',0_ip,1_ip,iposi_max)
  call parari('MAX',0_ip,1_ip,icolo_max)

  if( INOTMASTER ) then
     ncolo_loc = 0
     do icolo_glo = 1,npoin
        if( gisca(icolo_glo) > 0 ) ncolo_loc = ncolo_loc + 1 
     end do
     allocate( colum(ncolo_loc) )
     do icolo_loc = 1,ncolo_loc
        colum(icolo_loc) % nnode = 0
     end do

     icolo_loc = 0
     do icolo_glo = 1,npoin
        if( gisca(icolo_glo) > 0 ) then
           icolo_loc        = icolo_loc + 1
           pnode            = gisca(icolo_glo)
           gisca(icolo_glo) = icolo_loc
           allocate( colum(icolo_loc) % lnode(max(1_ip,pnode)) )
           allocate( colum(icolo_loc) % lposi(max(1_ip,pnode)) )
        end if
     end do

     do ipoin = 1,npoin
        if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then
           icolo_glo                       = int(xfiel(1) % a(1,ipoin,1),ip)
           icolo_loc                       = gisca(icolo_glo)
           colum(icolo_loc) % nnode        = colum(icolo_loc) % nnode + 1
           colum(icolo_loc) % iglob        = icolo_glo
           inode                           = colum(icolo_loc) % nnode
           colum(icolo_loc) % lnode(inode) = ipoin
           colum(icolo_loc) % lposi(inode) = int(xfiel(1) % a(2,ipoin,1),ip)
        end if
     end do

    do icolo_loc = 1,ncolo_loc
        pnode = colum(icolo_loc) % nnode
        colum(icolo_loc) % s = 0.0_rp
        if( pnode > 0 ) call hsort1(2_ip,pnode,colum(icolo_loc) % lposi,colum(icolo_loc) % lnode)
     end do

     allocate( my_paris(npoin) )
     allocate( my_parrs(npoin) )
     allocate( my_parin(npoin) )
     allocate( my_parre(npoin) )
     paris => my_paris   ! Snd
     parrs => my_parrs   ! Snd
     parin => my_parin   ! Rcv
     parre => my_parre   ! Rcv
     !
     ! NTOTA: my number of contributions
     !
     ntota = 0
     do icolo_loc = 1,ncolo_loc
        pnode = colum(icolo_loc) % nnode
        ntota = ntota + pnode
        colum(icolo_loc) % ilast = 1
     end do
     ntota_glo = ntota
     allocate( itest(npoin) )
     do ipoin = 1,npoin
        itest(ipoin) = 0
     end do
  end if
  
  call parari('SUM',0_ip,1_ip,ntota_glo)

  comcont   =  0
  itota     =  0
  inext     =  0
  itota_old = -1

  do while( comcont /= -1 )

     comcont = comcont + 1

     if( INOTMASTER ) then

        inext = 0

        do icolo_loc = 1,ncolo_loc

           icolo_glo = colum(icolo_loc) % iglob
           pnode     = colum(icolo_loc) % nnode
           inode     = 1
           igene     = 0

           do while( inode <= pnode )
              if( colum(icolo_loc) % ilast == colum(icolo_loc) % lposi(inode) ) then
                 ipoin = colum(icolo_loc) % lnode(inode)
                 if( itest(ipoin) == 0 ) then
                    itota                    = itota + 1
                    igene                    = 1
                    colum(icolo_loc) % ilast = colum(icolo_loc) % lposi(inode) + 1
                    itest(ipoin)             = 1
                    call cacapipi(ipoin,colum(icolo_loc) % s)
                 end if
              end if
              inode = inode + 1
           end do
           if( igene == 1 ) then
              inext          = inext + 1
              gisca(inext)   = icolo_glo
              inext          = inext + 1
              gisca(inext)   = colum(icolo_loc) % ilast 
              parrs(inext/2) = colum(icolo_loc) % s
           end if

        end do
        itota_glo = itota
     end if

     call parari('SUM',0_ip,1_ip,itota_glo)

     if( INOTMASTER .and. itota_glo == itota_old ) then
        print*,'je passe'
        do icolo_loc = 1,ncolo_loc
           icolo_glo      = colum(icolo_loc) % iglob
           inext          = inext + 1
           gisca(inext)   = icolo_glo
           inext          = inext + 1
           gisca(inext)   = colum(icolo_loc) % ilast 
           parrs(inext/2) = colum(icolo_loc) % s
        end do
     end if
     itota_old = itota_glo

if(imaster) print*,'a=',itota_glo , ntota_glo
if(comcont==20 ) call runend('SWSW')

     if( itota_glo == ntota_glo ) then
        !
        ! All columns are finished
        !
        comcont = -1     

     else if( ISLAVE ) then
        
        paris => my_paris          ! Snd int
        parrs => my_parrs          ! Snd real
        parin => my_parin          ! Rcv int
        parre => my_parre          ! Rcv real


!!! OJO CON ESTA COSA QUE ANTES ESTABA PERO NECESITA A DEF PARAL... NO SE COMO HAY QUE HACERLO

!!        do kfl_desti_par = 1,nneig
!           npari    = 0
!           nparr    = 0
!           npasi    = 1
!           npari    = 1
!           paris(1) = inext
!           call par_slaves(4_ip)
!           do inode = 1,inext
!              paris(inode) = gisca(inode)
!           end do
!           jnext    = parin(1)
!           npasi    = inext
!           npasr    = inext/2
!           npari    = jnext
!           nparr    = jnext/2
!           call par_slaves(4_ip)
!           knext = 0
!           do lnext = 1,jnext/2
!              knext     = knext + 1
!              icolo_glo = parin(knext)
!              knext     = knext + 1
!              kposi     = parin(knext)    
!          
!              icolo_loc = 1
!              loop_icolo: do while( colum(icolo_loc) % iglob /= icolo_glo )
!                 if( icolo_loc >= ncolo_loc ) exit loop_icolo
!                 icolo_loc = icolo_loc + 1
!              end do loop_icolo
!              if( colum(icolo_loc) % iglob == icolo_glo ) then
!                 colum(icolo_loc) % ilast = kposi  
!                 colum(icolo_loc) % s     = parre(lnext)
!              end if
!           end do
!        end do
     end if

  end do

  call memgen(0_ip,icolo_max,0_ip)
  if( INOTMASTER ) then
     do icolo_loc = 1,ncolo_loc
        icolo_glo = colum(icolo_loc) % iglob
        pnode     = colum(icolo_loc) % nnode
        if( pnode > 0 ) then
           if( colum(icolo_loc) % lposi(pnode) == iposi_max ) gesca(icolo_glo) = colum(icolo_loc) % s
        end if
     end do
  end if

  call pararr('SUM',0_ip,icolo_max,gesca)
  if( INOTSLAVE ) print*,'FINAL=',gesca
  call memgen(2_ip,icolo_max,0_ip)
  if( INOTMASTER ) then
     call memgen(3_ip,npoin,0_ip)
     allocate( my_paris(npoin) )
     allocate( my_parrs(npoin) )
     allocate( my_parin(npoin) )
     allocate( my_parre(npoin) )
  end if
  !if( INOTMASTER ) print*,kfl_paral,colum(:) % s
  call runend('SWSW')

end subroutine nsa_column

subroutine cacapipi(ipoin,s)
  use def_kintyp
  use def_domain
  implicit none
  integer(ip) :: ipoin
  real(rp)    :: s

  s=s+coord(2,ipoin)

end subroutine cacapipi
