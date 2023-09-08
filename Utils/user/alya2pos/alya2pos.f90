!-----------------------------------------------------------------------
!> @addtogroup Alya2pos
!> @{
!> @file    alya2pos.f90
!> @author  Mariano Vazquez
!> @brief   Postprocess tool
!> @date    16/11/1966
!> @details Postprocess tool
!> @}
!-----------------------------------------------------------------------
program alya2pos
  use def_kintyp
  use def_elmtyp
  use def_inpout
  use mod_output, only : output_gid_mesh
  use mod_maths
  implicit none
  real(rp),    pointer     :: gesca(:),gesc2(:)
  real(rp),    pointer     :: gevec(:,:),geve2(:,:),geve4(:,:,:)
  integer(ip), pointer     :: gisca(:),gisc2(:)
  integer(ip), pointer     :: givec(:,:)
  integer(ip)              :: ielem,ipoin,idime,inode,mnode,kpoin,jpoin,p1,p2,p3
  integer(ip)              :: nelem,npoin,ndime,ii,ieles,ipois,npart,nboun
  integer(ip)              :: mnodb,inodb,ipart,dummi,pdime,iboun,ibous
  integer(ip)              :: kboun,jboun,iblty,ielty,mpoin,iesta,iesto,igaus
  integer(ip)              :: npoin_total,nelem_total,nboun_total,idim1,idim2
  integer(ip)              :: npoin_2,nelem_2,nboun_2,mgaus,idim3
  integer(4)               :: ihead,ioerr
  integer(ip)              :: mjpoi
  real(rp),    pointer     :: coord_loc(:,:)
  real(rp)                 :: dummr,vec(3,3),ni,nj,nk,nn
  character(8)             :: chtim,varia
  character(80)            :: wstlb
  integer(ip)              :: ittim,timst,flag,mpoin_2
  real(rp)                 :: rttim
  logical(lg)              :: lopen,lgexi

  integer(ip)              :: kfl_markm,tag1,tag2
  integer(ip)              :: kfl_elimi
  integer(ip)              :: kfl_order
  integer(ip)              :: kfl_inter
  integer(ip)              :: kfl_bound
  integer(ip)              :: kfl_multi
  character(150)           :: forma
  logical(lg)              :: flag_coh

  integer(ip)              :: kfl_conve
  integer(ip)              :: npart_par
  real(rp)                 :: cutim
  character(5)             :: wwwww(10)
  character(8)             :: wwww8(10)
  integer(4)               :: iiiii(10)
  real(8)                  :: rrrrr(10)
  character(150)           :: namda,filna
  integer(ip), pointer     :: npoin_par(:)
  integer(ip), pointer     :: nelem_par(:)
  integer(ip), pointer     :: nboun_par(:)
  integer(ip), pointer     :: lsubd(:)

  integer(ip), pointer     :: markm(:)
  integer(ip), pointer     :: lnods(:,:)
  integer(ip), pointer     :: ltype(:)
  integer(ip), pointer     :: ltype_vtk(:)
  integer(ip), pointer     :: lelch(:)
  integer(ip), pointer     :: lesub(:)
  real(rp),    pointer     :: coord(:,:)
  integer(ip), pointer     :: lninv(:)
  integer(ip), pointer     :: leinv(:)
  integer(ip), pointer     :: leinv_cpy(:)
  integer(ip), pointer     :: lbinv(:)
  integer(ip), pointer     :: lbinv_cpy(:)
  integer(ip), pointer     :: ltypb(:)
  integer(ip), pointer     :: lnodb(:,:)
  integer(ip), pointer     :: lboch(:)
  real(rp),    pointer     :: gesc3(:),geve3(:,:)
  integer(ip), pointer     :: gisc3(:)
  integer(ip), pointer     :: lmate(:)
  character(13)            :: elemt

  integer(ip)              :: kfl_onlys,kauxi,jstep,kfl_stlbo
  integer(ip)              :: max_onlys, kfl_infil, kfl_field
  integer(ip), pointer     :: kfl_ppste(:)

  real(rp)                 :: small(3)
  real(rp)                 :: large(3)
  integer(ip)              :: resol(3),i

  ! restart
  character(5)               :: wopos(10)
  character(150)             :: wonam(10)
  integer(ip)                :: iwopo,kwopo
  integer(ip)                :: kfl_resta, kfofi,ndofn,idofn,kfl_ensbi
  integer(ip)                :: filt_msh_ensi,iavtk,lnods_vtk(90)

  !----------------------------------------------------------------------
  !
  ! Read argument + default options
  !
  !----------------------------------------------------------------------

  ! Nullify all pointers
  nullify(gesca)
  nullify(gesc2)
  nullify(gevec)
  nullify(geve2)
  nullify(geve4)
  nullify(gisca)
  nullify(gisc2)
  nullify(givec)
  nullify(coord_loc)
  nullify(npoin_par)
  nullify(nelem_par)
  nullify(nboun_par)
  nullify(lsubd)
  nullify(markm)
  nullify(lnods)
  nullify(ltype)
  nullify(ltype_vtk)
  nullify(lelch)
  nullify(lesub)
  nullify(coord)
  nullify(lninv)
  nullify(leinv)
  nullify(lbinv)
  nullify(leinv_cpy)
  nullify(lbinv_cpy)
  nullify(ltypb)
  nullify(lboch)
  nullify(lnodb)
  nullify(lmate)
  nullify(gesc3)
  nullify(geve3)
  nullify(gisc3)
  nullify(kfl_ppste)

  kfl_field= 0
  kfl_infil= 1
  call GETARG(1_4,namda)
  if( trim(namda) == '' ) then
     print*,'Wrong problem name'
  else if ( trim(namda) == '-p' ) then
     call GETARG(2_4,namda)
     forma= trim(namda)
     call GETARG(3_4,namda)
     if (trim(forma) /= 'vtk') then
        if (trim(forma) /= 'gid') then
           if (trim(forma) /= 'vu') then
              if (trim(forma) /= 'txt') then
                 if (trim(forma) /= 'ensight') then
                    if (trim(forma) /= 'zfem') then
                       if (trim(forma) == 'visit') then
                          forma = 'ensight'
                       else
                          print*,'Wrong postprocess program.'
                          print*,'   Must be gid , vu , ensight, visit, txt o zfem.'
                          stop
                       end if
                    end if
                 end if
              end if
           end if
        end if
     end if

     kfl_infil = 0

  else if ( trim(namda) == '-h' ) then
     print*,'--| '
     print*,'--| ALYA  alya2pos POSTPROCESS FORMATS CONVERTER'
     print*,'--| ALYA     Usage:'
     print*,'--| ALYA     Create a file yourproblem.post.alyadat '
     print*,'--| ALYA     as follows (for instance):'
     print*,'$-------------------------------------------------------------------'
     print*,'DATA'
     print*,'$  FORMAT:                   gid'
     print*,'  FORMAT:                   ensight'
     print*,'  MARK_ELEMENTS:            SUBDOMAIN'
     print*,'  ELIMINATE_BOUNDARY_NODES: Off'
     print*,'  MULTIPLE_FILE:            OFF'
     print*,'  BOUNDARY:                 ON'
     print*,'  SUBDOMAINS, ALL'
     print*,'  END_SUBDOMAINS'
     print*,'END_DATA'
     print*,'$-------------------------------------------------------------------'
     print*,'--| '
     print*,'--| ALYA     If the file is not created, postprocess files will be in GiD format'
     print*,'--| '
     print*,'--| ALYA     Finally, execute it as follows:   alya2pos.x yourproblem'
     print*,'--| ALYA   Done.'
     print*,'--| '

     stop

  end if

  kfl_markm = 3     ! 3=subdomain 0=element type, 4=element chacateristic, 6=dodeme
  kfl_elimi = 0     ! =1 eliminate boundary nodes
  kfl_bound = 1     ! =1 to postprocess boundary
  kfl_inter = 1     ! =1 to postprocess interior mesh
  kfl_multi = 0     ! =1 for multifile
  if (kfl_infil == 1) forma     = 'gid' ! Output format: gid, vu, filtre_vu, only when an alyadat file is to be read
  kfl_onlys = 0     ! /=0 ,only postprocess some steps independent of what alyafil says - avoid touching alyafil
  max_onlys = 3000000
  kfl_stlbo = 0
  kfl_resta = 0     ! restart
  kfl_order = 0     ! Original numbering for elements
  
  !----------------------------------------------------------------------
  !
  ! Header
  !
  !----------------------------------------------------------------------

  call livinf(0_ip,' ',0_ip)
  call livinf(1_ip,'CONVERT PROBLEM '//trim(namda),0_ip)
  call livinf(0_ip,' ',0_ip)

  !----------------------------------------------------------------------
  !
  ! Check ENDIAN
  !
  !----------------------------------------------------------------------

  kfl_conve = 0
  ihead = 1234
  open(unit=15,file=trim(namda)//'-COORD.post.alyabin',form='unformatted')
  read(15,err=1111) ihead
  if( ihead /= 1234 ) then
     call livinf(1_ip,'Converting ENDIAN',0_ip)
     goto 1111
  end if
  goto 2222
1111 kfl_conve = 1
2222 continue
  close(unit=15)
  if( ihead /= 1234 ) kfl_conve = 1

  !----------------------------------------------------------------------
  !
  ! Run data: subdomain dimensions
  !
  !----------------------------------------------------------------------

  open(unit=8,file=trim(namda)//'.post.alyapar',form='formatted')
  read(8,*) npart_par
  allocate( nelem_par(npart_par) )
  allocate( npoin_par(npart_par) )
  allocate( nboun_par(npart_par) )
  allocate( lsubd(npart_par) )
  do ipart = 1,npart_par
     read(8,*) dummi,nelem_par(ipart),npoin_par(ipart),nboun_par(ipart)
     lsubd(ipart) = 1
  end do
  close(8)

  !----------------------------------------------------------------------
  !
  ! Read data file (only if kfl_infil = 1)
  !
  !----------------------------------------------------------------------

  open(unit=10,file=trim(namda)//'.post.alyadat',status='old',form='formatted',iostat=ioerr)

  if( ioerr == 0 ) then

     lispa =  0      ! 0 passes through ecoute
     lisda = 10      ! Temporary data file
     lisre =  6      ! Results file
     endst =  1      ! Stop Alya of end of file found

     if (kfl_infil == 1) then

        call ecoute('ALYA2P')
        if( words(1) /= 'DATA ' ) call runend('ALYA2POS: WRONG DATA CARD')
        do while( words(1) /= 'ENDDA' )
           if( words(1) == 'FORMA' ) then
              if( words(2) == 'GID  ' ) forma = 'gid'
              if( words(2) == 'VU   ' ) forma = 'vu'
              if( words(2) == 'ENSIG' ) forma = 'ensight'
              if( words(2) == 'TXT  ' ) forma = 'txt'
              if( words(2) == 'VISIT' ) forma = 'ensight'
              if( words(2) == 'PARAV' ) forma = 'ensight'
              if( words(2) == 'FILTR' ) forma = 'filtre_vu'
              if( words(2) == 'ALYA ' ) forma = 'alya'
              if( words(2) == 'VTK  ' ) forma = 'vtk'
              if( words(2) == 'ZFEM ' ) forma = 'zfem'
              if( words(2) == 'FIELD' ) then
                 forma = 'ensight'
                 kfl_field=1
              end if
           else if( words(1) == 'ORDER' ) then
              if( words(2) == 'YES  ' ) kfl_order = 1
              if( words(2) == 'ON   ' ) kfl_order = 1
           else if( words(1) == 'INTER' ) then
              if( words(2) == 'OFF  ' ) kfl_inter = 0
              if( words(2) == 'NO   ' ) kfl_inter = 0
           else if( words(1) == 'MARKE' ) then
              if( words(2) == 'TYPE ' ) kfl_markm = 0
              if( words(2) == 'SUBDO' ) kfl_markm = 3
              if( words(2) == 'PARTI' ) kfl_markm = 3
              if( words(2) == 'CHARA' ) kfl_markm = 4
              if( words(2) == 'MATER' ) kfl_markm = 5
              if( words(2) == 'DODEM' ) kfl_markm = 6
              if( words(2) == 'BEATR' ) kfl_markm = 7
              if( words(2) == 'LBOCH' ) kfl_markm = 8
           else if( words(1) == 'STL  ' ) then
              if( words(2) == 'ON   ' ) kfl_stlbo = 1
              if( words(2) == 'YES  ' ) kfl_stlbo = 1
              if( exists('BINAR') )     kfl_stlbo = 2
           else if( words(1) == 'ELIMI' ) then
              if( words(2) == 'ON   ' ) kfl_elimi = 1
              if( words(2) == 'YES  ' ) kfl_elimi = 1
              if( words(2) == 'NO   ' ) kfl_elimi = 0
              if( words(2) == 'OFF  ' ) kfl_elimi = 0
           else if( words(1) == 'BOUND' ) then
              if( words(2) == 'ON   ' ) kfl_bound = 1
              if( words(2) == 'YES  ' ) kfl_bound = 1
              if( words(2) == 'NO   ' ) kfl_bound = 0
              if( words(2) == 'OFF  ' ) kfl_bound = 0
           else if( words(1) == 'MULTI' ) then
              if( words(2) == 'ON   ' ) kfl_multi = 1
              if( words(2) == 'YES  ' ) kfl_multi = 1
              if( words(2) == 'NO   ' ) kfl_multi = 0
              if( words(2) == 'OFF  ' ) kfl_multi = 0
           else if( words(1) == 'SUBDO' ) then
              if( words(2) == 'ALL  ') then
                 do while( words(1) /= 'ENDSU' )
                    call ecoute('ALYA2P')
                 end do
              else
                 do ipart = 1,npart_par
                    lsubd(ipart) = 0
                 end do
                 call ecoute('ALYA2P')
                 do while( words(1) /= 'ENDSU' )
                    ipart = int(param(1))
                    if( ipart < 1 .or. ipart > npart_par ) call runend('ALYA2POS: WRONG SUBDOMAIN NUMBER')
                    lsubd(ipart) = 1
                    call ecoute('ALYA2P')
                 end do
              end if
           else if( words(1) == 'ONLYS' ) then
              if( words(2) == 'ALL  ') then
                 kfl_onlys = 0
                 do while( words(1) /= 'ENDSU' )
                    call ecoute('ALYA2P')
                 end do
              else
                 kfl_onlys = 1
                 allocate (kfl_ppste(0:max_onlys))
                 do kauxi = 0,max_onlys
                    kfl_ppste(kauxi) = 0
                 end do
                 call ecoute('ALYA2P')
                 do while( words(1) /= 'ENDON' )
                    kauxi= int(param(1))
                    if( kauxi < 1 .or. kauxi > max_onlys ) call runend('ALYA2POS: WRONG STEP IN ONLY_STEPS')
                    kfl_ppste(kauxi) = 1
                    call ecoute('ALYA2P')
                 end do
              end if

!!!!!!! CREATE RESTART FILES FROM GIVEN FIELDS

           else if( words(1) == 'RESTA' ) then
              kfl_resta= 1
              kfofi= 0                        ! binary format
              if(exists('ASCII')) kfofi= 1    ! ascii format
              iwopo= 0
              call ecoute('ALYA2P')
              do while( words(1) /= 'ENDRE' )
                 iwopo= iwopo+1
                 wopos(iwopo) = words(1)
                 wonam(iwopo) = adjustl(trim(words(2)))
                 call ecoute('ALYA2P')
              end do
              kwopo = iwopo
           end if
           call ecoute('ALYA2P')
        end do

     end if
     close(unit=10)

  else

     call livinf(1_ip,'DATA FILE IS MISSING: '//trim(namda)//'.post.alyadat. CHOOSE DEFAULT OPTIONS',0_ip)

  end if

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------
  call elmtyp()
  do ii = 1,nelty
     lexis(ii) = 0
     lbxis(ii) = 0
  end do
  jttim = -1
  if( npart_par == 1 ) kfl_elimi = 0
  if( npart_par == 1 ) lsubd = 1
  nboun_total = 0

  !----------------------------------------------------------------------
  !
  ! Geometry
  !
  !----------------------------------------------------------------------

  open(unit=49,file=trim(namda)//'.post.alyafil',form=  'formatted')

  do ii = 15,26

     if( ii == 15 ) then
        !
        ! LTYPE
        !
        if( kfl_conve == 0 ) then
           open(unit=15,file=trim(namda)//'-LTYPE.post.alyabin',form='unformatted')
        else
           open(unit=15,file=trim(namda)//'-LTYPE.post.alyabin',form='unformatted',convert='BIG_ENDIAN')
        end if
        call reahed(ii,npart_par,wwww8,iiiii,rrrrr)

        ieles = 0
        nelem_total = iiiii(2)
        call livinf(2_ip,'LTYPE',0_ip)
        allocate(  ltype(nelem_total) )

        do ipart = 1,npart_par
           if( lsubd(ipart) == 1 ) then
              read(ii) nelem
              read(ii) ( ltype(ielem),ielem=ieles+1,ieles+nelem)
              do ielem = ieles+1,ieles+nelem
                 ielty = abs(ltype(ielem))
                 lexis(ielty) = lexis(ielty) + 1
              end do
              ieles = ieles + nelem_par(ipart)
           else
              read(ii) nelem
              read(ii) ( dummi,ielem=1,nelem)
           end if
        end do

        close(unit=15)
        call livinf(3_ip,'LTYPE',0_ip)

     else if( ii == 16 ) then
        !
        ! LNODS
        !
        if( kfl_conve == 0 ) then
           open(unit=16,file=trim(namda)//'-LNODS.post.alyabin',form='unformatted')
        else
           open(unit=16,file=trim(namda)//'-LNODS.post.alyabin',form='unformatted',convert='BIG_ENDIAN')
        end if
        call reahed(ii,npart_par,wwww8,iiiii,rrrrr)

        call livinf(2_ip,'LNODS',0_ip)
        ieles = 0
        ipois = 0
        nelem_total = iiiii(2)
        mnode = iiiii(1)
        allocate(  lnods(mnode,nelem_total) )

        do ipart = 1,npart_par
           if( lsubd(ipart) == 1 ) then
              read(ii) nelem
              read(ii) ( (lnods(inode,ielem),inode=1,mnode),ielem=ieles+1,ieles+nelem)
              do ielem = ieles+1,ieles+nelem
                 do inode = 1,nnode(ltype(ielem))
                    lnods(inode,ielem) = lnods(inode,ielem) + ipois
                 end do
              end do
              ieles = ieles + nelem_par(ipart)
              ipois = ipois + npoin_par(ipart)
           else
              read(ii) nelem
              read(ii) ( (dummi,inode=1,mnode),ielem=1,nelem)
           end if

        end do

        close(unit=16)
        call livinf(3_ip,'LNODS',0_ip)

     else if( ii == 17 ) then
        !
        ! COORD
        !
        if( kfl_conve == 0 ) then
           open(unit=17,file=trim(namda)//'-COORD.post.alyabin',form='unformatted')
        else
           open(unit=17,file=trim(namda)//'-COORD.post.alyabin',form='unformatted',convert='BIG_ENDIAN')
        end if
        call reahed(ii,npart_par,wwww8,iiiii,rrrrr)

        ipois = 0
        npoin_total = iiiii(2)
        ndime = iiiii(1)
        call livinf(2_ip,'COORD',0_ip)
        allocate(  coord_loc(ndime,npoin_total) )

        do ipart = 1,npart_par
           if( lsubd(ipart) == 1 ) then
              read(ii) npoin
              read(ii) ( (coord_loc(idime,ipoin),idime=1,ndime),ipoin=ipois+1,ipois+npoin)
              ipois = ipois + npoin_par(ipart)
           else
              read(ii) npoin
              read(ii) ( (dummi,idime=1,ndime),ipoin=1,npoin)
           end if
        end do

        close(unit=17)
        call livinf(3_ip,'COORD',0_ip)

     else if( ii == 18 .and. npart_par > 1 ) then
        !
        ! LNINV
        !
        if( kfl_conve == 0 ) then
           open(unit=18,file=trim(namda)//'-LNINV.post.alyabin',form='unformatted')
        else
           open(unit=18,file=trim(namda)//'-LNINV.post.alyabin',form='unformatted',convert='BIG_ENDIAN')
        end if
        call reahed(ii,npart_par,wwww8,iiiii,rrrrr)

        ipois = 0
        npoin_total = iiiii(2)
        call livinf(2_ip,'LNINV',0_ip)
        allocate(  lninv(npoin_total) )

        do ipoin = 1,npoin_total
           lninv(ipoin) = 0
        end do
        do ipart = 1,npart_par
           if( lsubd(ipart) == 1 ) then
              read(ii) npoin
              read(ii) ( lninv(ipoin),ipoin=ipois+1,ipois+npoin)
              ipois = ipois + npoin_par(ipart)
           else
              read(ii) npoin
              read(ii) ( dummi,ipoin=1,npoin)
           end if
        end do

        close(unit=18)
        call livinf(3_ip,'LNINV',0_ip)

     else if( ii == 19 ) then
        !
        ! LELCH
        !
        if( kfl_conve == 0 ) then
           open(unit=19,file=trim(namda)//'-LELCH.post.alyabin',form='unformatted')
        else
           open(unit=19,file=trim(namda)//'-LELCH.post.alyabin',form='unformatted',convert='BIG_ENDIAN')
        end if
        call reahed(ii,npart_par,wwww8,iiiii,rrrrr)

        ieles = 0
        nelem_total = iiiii(2)
        call livinf(2_ip,'LELCH',0_ip)
        allocate(  lelch(nelem_total) )

        do ipart = 1,npart_par
           read(ii) nelem
           read(ii) ( lelch(ielem),ielem=ieles+1,ieles+nelem)
           do ielem = ieles+1,ieles+nelem
              lelch(ielem) = lelch(ielem)+10
           end do
           ieles = ieles + nelem_par(ipart)
        end do

        close(unit=19)
        call livinf(3_ip,'LELCH',0_ip)

     else if( ii == 20 .and. kfl_bound == 1 ) then
        !
        ! LTYPB
        !
        if( kfl_conve == 0 ) then
           open(unit=20,file=trim(namda)//'-LTYPB.post.alyabin',form='unformatted',status='old',iostat=ioerr)
        else
           open(unit=20,file=trim(namda)//'-LTYPB.post.alyabin',form='unformatted',status='old',iostat=ioerr,convert='BIG_ENDIAN')
        end if
        if( ioerr /= 0 ) then
           kfl_bound = 0
           goto 998
        end if
        call reahed(ii,npart_par,wwww8,iiiii,rrrrr)

        ibous = 0
        nboun_total = iiiii(2)
        call livinf(2_ip,'LTYPB',0_ip)
        allocate(  ltypb(nboun_total) )

        do ipart = 1,npart_par
           if( lsubd(ipart) == 1 ) then
              read(ii) nboun
              read(ii) ( ltypb(iboun),iboun=ibous+1,ibous+nboun)
              do iboun = ibous+1,ibous+nboun
                 iblty = abs(ltypb(iboun))
                 lbxis(iblty) = lbxis(iblty) + 1
              end do
              ibous = ibous + nboun_par(ipart)
           else
              read(ii) nboun
              read(ii) ( dummi,iboun=1,nboun)
           end if
        end do

        close(unit=20)
        call livinf(3_ip,'LTYPB',0_ip)

     else if( ii == 21 .and. kfl_bound == 1 ) then
        !
        ! LNODB
        !
        if( kfl_conve == 0 ) then
           open(unit=21,file=trim(namda)//'-LNODB.post.alyabin',form='unformatted')
        else
           open(unit=21,file=trim(namda)//'-LNODB.post.alyabin',form='unformatted',convert='BIG_ENDIAN')
        end if
        call reahed(ii,npart_par,wwww8,iiiii,rrrrr)

        call livinf(2_ip,'LNODB',0_ip)
        ibous = 0
        ipois = 0
        nboun_total = iiiii(2)
        mnodb = iiiii(1)
        allocate(  lnodb(mnodb,nboun_total) )

        do ipart = 1,npart_par
           if( lsubd(ipart) == 1 ) then
              read(ii) nboun
              read(ii) ( (lnodb(inodb,iboun),inodb=1,mnodb),iboun=ibous+1,ibous+nboun)
              do iboun = ibous+1,ibous+nboun
                 do inodb = 1,nnode(ltypb(iboun))
                    lnodb(inodb,iboun) = lnodb(inodb,iboun) + ipois
                 end do
              end do
              ibous = ibous + nboun_par(ipart)
              ipois = ipois + npoin_par(ipart)
           else
              read(ii) nboun
              read(ii) ( (dummi,inodb=1,mnodb),iboun=1,nboun)
           end if

        end do

        !write(*,*)'nboun_total=',nboun_total
        !write(*,*)'mnodb=',mnodb
        !write(*,*)''
        !do iboun=1,nboun_total
        !   write(*,*)'lnodb(inodb,iboun)',(lnodb(inodb,iboun),inodb=1,mnodb)
        !end do

        close(unit=21)
        call livinf(3_ip,'LNODB',0_ip)

     else if( ii == 22 ) then
        !
        ! LESUB
        !
        lgexi = .false.
        inquire(file=trim(namda)//'-LESUB.post.alyabin',exist=lgexi)
        if( lgexi ) then
           if( kfl_conve == 0 ) then
              open(unit=22,file=trim(namda)//'-LESUB.post.alyabin',form='unformatted')
           else
              open(unit=22,file=trim(namda)//'-LESUB.post.alyabin',form='unformatted',convert='BIG_ENDIAN')
           end if

           call reahed(ii,npart_par,wwww8,iiiii,rrrrr)

           ieles = 0
           nelem_total = iiiii(2)
           call livinf(2_ip,'LESUB',0_ip) 
           allocate(  lesub(nelem_total) )

           do ipart = 1,npart_par
              if( lsubd(ipart) == 1 ) then
                 read(ii) nelem
                 read(ii) ( lesub(ielem),ielem=ieles+1,ieles+nelem)
                 ieles = ieles + nelem_par(ipart)
              else
                 read(ii) nelem
                 read(ii) ( dummi,ielem=1,nelem)
              end if
           end do
           
           !do ipart = 1,npart_par
           !   read(ii) nelem
           !   read(ii) ( lesub(ielem),ielem=ieles+1,ieles+nelem)
           !   !print*,lesub(ieles+1:ieles+nelem)
           !   ieles = ieles + nelem_par(ipart)
           !end do
           !do ipart = 1,npart_par
           !   read(ii) nelem
           !   allocate(gesca(nelem))
           !   read(ii) ( gesca(ielem),ielem=1,nelem)
           !   do ielem = 1,nelem
           !      lesub(ieles+ielem) = int(gesca(ielem),rp)
           !   end do
           !   deallocate(gesca)
           !   ieles = ieles + nelem_par(ipart)
           !end do

           close(unit=22)
           call livinf(3_ip,'LESUB',0_ip)

        end if

     else if( ii == 23 ) then
        !
        ! LEINV
        !
        if( kfl_conve == 0 ) then
           open(unit=23,file=trim(namda)//'-LEINV.post.alyabin',form='unformatted',status='old',iostat=ioerr)
        else
           open(unit=23,file=trim(namda)//'-LEINV.post.alyabin',form='unformatted',status='old',convert='BIG_ENDIAN',iostat=ioerr)
        end if
        if( forma /= 'gid' ) call livinf(1_ip,'WARNING! ELEMENT NUMBERING IS NOT THE ORIGINAL',0_ip)
        if( ioerr /= 0 ) then
           goto 997
        end if

        call reahed(ii,npart_par,wwww8,iiiii,rrrrr)

        ieles = 0
        nelem_total = iiiii(2)
        call livinf(2_ip,'LEINV',0_ip)
        allocate(  leinv(nelem_total) )

        do ipart = 1,npart_par
           read(ii) nelem
           read(ii) ( leinv(ielem),ielem=ieles+1,ieles+nelem)
           ieles = ieles + nelem_par(ipart)
        end do

        close(unit=23)
        call livinf(3_ip,'LEINV',0_ip)

997     continue

     else if( ii == 24 ) then
        !
        ! LMATE
        !
        if( kfl_conve == 0 ) then
           open(unit=24,file=trim(namda)//'-LMATE.post.alyabin',form='unformatted',status='old',iostat=ioerr)
        else
           open(unit=24,file=trim(namda)//'-LMATE.post.alyabin',form='unformatted',status='old',convert='BIG_ENDIAN',iostat=ioerr)
        end if

        !call livinf(1_ip,'WARNING! ELEMENT NUMBERING IS NOT THE ORIGINAL',0_ip)
        if( ioerr /= 0 ) then
           call livinf(1_ip,'WARNING! MATERIAL FILE NOT AVAILABLE',0_ip)
           goto 1999
        end if

        call reahed(ii,npart_par,wwww8,iiiii,rrrrr)

        ieles = 0
        nelem_total = iiiii(2)
        call livinf(2_ip,'LMATE',0_ip)
        allocate(  lmate(nelem_total) )

        do ipart = 1,npart_par
           read(ii) nelem
           read(ii) ( lmate(ielem),ielem=ieles+1,ieles+nelem)
           ieles = ieles + nelem_par(ipart)
        end do

        close(unit=24)
        call livinf(3_ip,'LMATE',0_ip)

1999    continue

     else if( ii == 25 .and. kfl_bound == 1 ) then
        !
        ! LBOCH
        !
        if( kfl_conve == 0 ) then
           open(unit=25,file=trim(namda)//'-LBOCH.post.alyabin',form='unformatted',status='old',iostat=ioerr)
        else
           open(unit=25,file=trim(namda)//'-LBOCH.post.alyabin',form='unformatted',status='old',iostat=ioerr,convert='BIG_ENDIAN')
        end if
        if( ioerr /= 0 ) then
           goto 998
        end if
        call reahed(ii,npart_par,wwww8,iiiii,rrrrr)

        ibous = 0
        nboun_total = iiiii(2)
        call livinf(2_ip,'LBOCH',0_ip)
        allocate(  lboch(nboun_total) )

        do ipart = 1,npart_par
           if( lsubd(ipart) == 1 ) then
              read(ii) nboun
              read(ii) ( lboch(iboun),iboun=ibous+1,ibous+nboun)
              ibous = ibous + nboun_par(ipart)
           else
              read(ii) nboun
              read(ii) ( dummi,iboun=1,nboun)
           end if
        end do

        close(unit=25,iostat=ioerr)
        if( ioerr /= 0 ) call runend('COULD NOT CLOSE FILE')
        call livinf(3_ip,'LBOCH',0_ip)

     else if( ii == 26 .and. kfl_bound == 1 ) then
        !
        ! LBINV
        !
        if( kfl_conve == 0 ) then
           open(unit=26,file=trim(namda)//'-LBINV.post.alyabin',form='unformatted',status='old',iostat=ioerr)
        else
           open(unit=26,file=trim(namda)//'-LBINV.post.alyabin',form='unformatted',status='old',iostat=ioerr,convert='BIG_ENDIAN')
        end if
        if( ioerr /= 0 ) then
           goto 998
        end if
        call reahed(ii,npart_par,wwww8,iiiii,rrrrr)

        ibous = 0
        nboun_total = iiiii(2)
        call livinf(2_ip,'LBINV',0_ip)
        allocate(  lbinv(nboun_total) )
        
        do ipart = 1,npart_par
           if( lsubd(ipart) == 1 ) then
              read(ii) nboun
              read(ii) ( lbinv(iboun),iboun=ibous+1,ibous+nboun)
              ibous = ibous + nboun_par(ipart)
           else
              read(ii) nboun
              read(ii) ( dummi,iboun=1,nboun)
           end if
        end do

        close(unit=26,iostat=ioerr)
        if( ioerr /= 0 ) call runend('COULD NOT CLOSE FILE')
        call livinf(3_ip,'LBINV',0_ip)

     end if

998  continue

  end do

  !----------------------------------------------------------------------
  !
  ! Eliminate repeated nodes
  !
  !----------------------------------------------------------------------

  if( kfl_elimi == 1 .and. npart_par > 0 ) then

     allocate( gisca(npoin_total) )
     do ipoin = 1,npoin_total
        gisca(ipoin) = 0
     end do
     kpoin = 0
     mjpoi = 0
     do ipoin = 1,npoin_total
        jpoin = lninv(ipoin)
        if(jpoin>mjpoi) mjpoi=jpoin  !debugging herbert
        if( jpoin < 0 ) then
           write(*,*) 'Wrong mesh'
           stop
        else if( gisca(jpoin) == 0 ) then
           kpoin = kpoin + 1
           gisca(jpoin) = lninv(ipoin)
        end if
     end do
     ! debugging herbert
     !write(*,*)'mjpoi',mjpoi
     do jpoin = 1,mjpoi
        if(gisca(jpoin)==0) write(*,*)'gisca(jpoin)==0 at point',&
             jpoin,'it does not belong to any element: problem with the mesh!!'
     end do

     ! end debugging herbert
     deallocate( gisca )

     npoin_total = mjpoi

     allocate(  coord(ndime,npoin_total) )

     kpoin = 0
     do ipart = 1,npart_par
        do ipoin = 1,npoin_par(ipart)
           kpoin = kpoin + 1
           jpoin = lninv(kpoin)
           do idime = 1,ndime
              coord(idime,jpoin) = coord_loc(idime,kpoin)
           end do
        end do
     end do

     deallocate( coord_loc )

     do ielem = 1,nelem_total
        do inode = 1,nnode(abs(ltype(ielem)))
           ipoin = lnods(inode,ielem)
           lnods(inode,ielem) = lninv(ipoin)
        end do
     end do

     do iboun = 1,nboun_total
        do inodb = 1,nnode(abs(ltypb(iboun)))
           ipoin = lnodb(inodb,iboun)
           lnodb(inodb,iboun) = lninv(ipoin)
        end do
     end do

  else

     coord => coord_loc

  end if

  !----------------------------------------------------------------------
  !
  ! Boundary filter
  !
  !----------------------------------------------------------------------

  goto 88
  jboun = 0
  iboun = 0
  do ipart = 1,npart_par
     kboun = 0
     do ibous = 1,nboun_par(ipart)
        dummi = 0
        iboun = iboun + 1
        do inodb = 1,nnode(abs(ltypb(iboun)))
           ipoin = lnodb(inodb,iboun)
           if(        coord(1,ipoin) > -0.2_rp .and. coord(1,ipoin) < 3.1_rp &
                .and. coord(2,ipoin) > -0.3_rp .and. coord(2,ipoin) < 0.3_rp &
                .and. coord(3,ipoin) > -0.3_rp .and. coord(3,ipoin) < 0.3_rp ) then
              dummi = 1
           end if
        end do
        if( dummi == 1 ) then
           kboun = kboun + 1
           jboun = jboun + 1
           do inodb = 1,nnode(abs(ltypb(iboun)))
              lnodb(inodb,jboun) = lnodb(inodb,iboun)
           end do
        end if
     end do
     nboun_par(ipart) = kboun
  end do
88 continue

  !----------------------------------------------------------------------
  !
  ! Output geometry
  !
  !----------------------------------------------------------------------

  npoin_2 = 0
  nelem_2 = 0
  nboun_2 = 0
  do ipart = 1,npart_par
     if( lsubd(ipart) == 1 ) then
        npoin_2 = npoin_2 + npoin_par(ipart)
        nelem_2 = nelem_2 + nelem_par(ipart)
        nboun_2 = nboun_2 + nboun_par(ipart)
     end if
  end do

  !
  ! Initialization flag for cohesive element detection
  !
  flag_coh = .False.
  do ielem = 1,nelem_total
     if( lelch(ielem) == 17_ip) then
        flag_coh = .True.
        exit
     end if
  end do

  if( kfl_elimi == 1 ) npoin_2 = npoin_total

  if( forma == 'vtk' ) then
     open(unit=100,file=trim(namda)//'.msh.vtu',form='formatted')
     write (6,*) '     WARNING!! VTK FORMAT ONLY EXPORTS THE MESH'
     write (100,*) &
          '    <VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
     write (100,*) &
          '      <UnstructuredGrid>0'
     write (100,*) &
          '        <Piece NumberOfPoints="',&
          trim(intost(npoin)),&
          '" NumberOfCells="',&
          trim(intost(nelem)),&
          '">'
     write (100,*) &
          '          <Points>'
     write (100,*) &
          '            <DataArray type="Float32" NumberOfComponents="',&
          trim(intost(ndime)),&
          '" Format="ascii">'
     do ipoin=1,npoin
        write(100,  1020) coord(ipoin,1:ndime)
     end do
     write (100,*) &
          '            </DataArray>'
     write (100,*) &
          '          </Points>'

     allocate(ltype_vtk(nelem))

     iavtk= 0
     do ielem = 1,nelem
        if (ltype(ielem) == 30) then
           ltype_vtk(ielem)=10
        else if (ltype(ielem) == -34) then
           ltype_vtk(ielem)=12
        else
           write(6,*) 'VTK TYPE NOT YET MATCHED: ',ltype(ielem)
           stop
        end if
        iavtk= iavtk + nnode(ltype(ielem))
     end do
     write (100,*) &
          '          <Cells>'
     ! lnods
     write (100,*) &
          '            <DataArray type="Int32" Name="connectivity" Format="ascii">'
     do ielem = 1,nelem
        ! recall that vtk starts numbering by 0
        lnods_vtk=lnods(1:nnode(ltype(ielem)),ielem)
        lnods_vtk=lnods_vtk - 1
        write(100,1030) lnods_vtk(1:nnode(ltype(ielem)))
     end do
     write (100,*) &
          '            </DataArray>'
     ! nnode
     write (100,*) &
          '            <DataArray type="Int32" Name="offsets" Format="ascii">'
     iavtk=0
     do ielem = 1,nelem
        iavtk= iavtk + nnode(ltype(ielem))
        write(100,1030) iavtk
     end do
     write (100,*) &
          '            </DataArray>'
     ! vtk types
     write (100,*) &
          '            <DataArray type="Int32" Name="types" Format="ascii">'
     do ielem = 1,nelem
        write(100,*) ltype_vtk(ielem)
     end do
     write (100,*) &
          '            </DataArray>'
     write (100,*) &
          '          </Cells>'
     write (100,*) &
          '     </Piece>'
     write (100,*) &
          '  </UnstructuredGrid>'
     write (100,*) &
          '</VTKFile>'

     write (6,*) '     MESH EXPORTED!'
     stop

1020 format(10(1x,e16.8E3))
1030 format(20(1x,i8))

  else if( forma == 'gid' ) then
     !
     ! Gid format
     !
     !
     ! Create mesh file (.post.msh)
     !
     open(unit=100,file=trim(namda)//'.post.msh',form='formatted')
     call livinf(1_ip,'CONVERT GEOMETRY',0_ip)

     if( kfl_markm == 0 ) then
        markm => ltype
     else if( kfl_markm == 1 .or. kfl_markm == 3 ) then
        allocate( markm(nelem_2) )
        ieles = 0
        do ipart = 1,npart_par
           if( lsubd(ipart) == 1 ) then
              do ielem = 1,nelem_par(ipart)
                 ieles = ieles + 1
                 markm(ieles) = ipart
              end do
           end if
        end do
     else if( kfl_markm == 4 ) then
        markm => lelch
     else if( kfl_markm == 5 ) then
        markm => lmate
     else if( kfl_markm == 6 .or. kfl_markm == 7 ) then
        markm => lesub
        if( .not. associated(lesub) ) call runend('CANNOT POSTPROCESS USING DODEME')
     else if ( kfl_markm == 8 ) then
        markm => lboch
     end if

     ipois = 0
     ieles = 0
     if( kfl_inter == 1 ) then
        !
        ! Reroder elements (mainly for debugging purpose)
        !
        if( kfl_order == 1 ) then
           allocate(leinv_cpy(nelem_total))
           leinv_cpy(1:nelem_total) = leinv(1:nelem_total) 
           call maths_heapsort_int(2_ip,mnode,nelem_total,leinv,lnods)
        end if        
        call output_gid_mesh(mnode,npoin_2,nelem_2,100_ip,lexis,ltype,lnods, &
             coord,markm,lelch,leinv,ieles,ipois,ndime,ndime,kfl_markm,   &
             flag_coh,namda)
        !if( kfl_order == 1 ) then
        !   leinv(1:nelem_total) = leinv_cpy(1:nelem_total)
        !   deallocate(leinv_cpy)
        !end if
        
     end if
     if( kfl_bound == 1 ) then
        !if( associated(lboch) ) then
        !   markm => lboch
        !else
        markm => ltypb
        !end if
        ipois =  0
        if( kfl_order == 1 ) then
           allocate(lbinv_cpy(nboun_total))
           lbinv_cpy(1:nboun_total) = lbinv(1:nboun_total) 
           call maths_heapsort_int(2_ip,mnodb,nboun_total,lbinv,lnodb)
        end if        
        call output_gid_mesh(mnodb,npoin,nboun_2,100_ip,lbxis,ltypb,lnodb,  &
             coord,markm,lboch,lbinv,ieles,ipois,ndime,ndime-1,kfl_markm, &
             flag_coh,trim(namda)//'_BOU')
     end if

     if( kfl_markm > 0 ) deallocate( markm )

     close(100)

  else if ( forma == 'alya' ) then

     goto 99999
     open(unit=100,file=trim(namda)//'.geo.dat',form='formatted')
     write(100,'(a)') 'TYPES, ALL=TET04'
     write(100,'(a)') 'END_TYPES'
     write(100,'(a)') 'COORDINATES'
     do ipoin = 1,npoin_2
        write(100,'(i11,3(1x,e13.6))') ipoin,coord(1,ipoin),coord(2,ipoin),coord(3,ipoin)
     end do
     write(100,'(a)') 'END_COORDINATES'
     write(100,'(a)') 'ELEMENTS'
     do ielem = 1,nelem_2
        write(100,'(i11,20(1x,i11))') ielem,(lnods(inode,ielem),inode=1,nnode(ltype(ielem)))
     end do
     write(100,'(a)') 'END_ELEMENTS'
     write(100,'(a)') 'BOUNDARIES'
     do iboun = 1,nboun_2
        write(100,'(i11,20(1x,i11))') iboun,(lnodb(inodb,iboun),inodb=1,nnode(ltypb(iboun)))
     end do
     write(100,'(a)') 'END_BOUNDARIES'
99999 continue

  else if ( forma == 'vu' ) then
     !
     ! VU format
     !
     open(unit=100,file=trim(namda)//'.msh.vu',form='formatted')
     open(unit=101,file=trim(namda)//'.msh.vubin',form='unformatted')

     if( kfl_bound == 1 ) then
        call vu_msh(&
             kfl_bound,mnode,mnodb,npoin_2,nelem_2,nboun_2,100_ip,&
             101_ip,lexis,ltype,lnods,lbxis,ltypb,lnodb,lelch,coord,&
             ndime,namda,kfl_markm,npart_par,lsubd,nelem_par)
     else
        call vu_msh(&
             kfl_bound,mnode,mnodb,npoin_2,nelem_2,nboun_2,100_ip,&
             101_ip,lexis,ltype,lnods,dummi,dummi,dummi,lelch,coord,&
             ndime,namda,kfl_markm,npart_par,lsubd,nelem_par)
     end if

     close(100)
     close(101)

  else if ( forma == 'ensight' ) then
     !
     ! Ensight format
     !
     kfl_ensbi = 0_ip
     !    kfl_ensbi = 1_ip
     open(unit=100,file=trim(namda)//'.ensi.case',form='formatted')

     if ( kfl_ensbi == 0 ) then
        open(unit=101,file=trim(namda)//'.ensi.geo', form='formatted')
        if (kfl_markm == 4) open(unit=202,file=trim(namda)//'.ensi.LELCH', form='formatted')
        if( kfl_bound == 1 ) then
           call ensmsh(&
                kfl_bound,mnode,mnodb,npoin_2,nelem_2,nboun_2,100_ip,&
                101_ip,lexis,ltype,leinv,lnods,lbxis,ltypb,lnodb,lelch,coord,&
                ndime,namda,kfl_markm,npart_par,lsubd,nelem_par)
        else
           call ensmsh(&
                kfl_bound,mnode,mnodb,npoin_2,nelem_2,nboun_2,100_ip,&
                101_ip,lexis,ltype,leinv,lnods,dummi,dummi,dummi,lelch,coord,&
                ndime,namda,kfl_markm,npart_par,lsubd,nelem_par)
        end if
     else

!!! ONLY FOR INTEL (recordtype)

        !        open(unit=101,file=trim(namda)//'.ensi.geo',&
        !             &     form='unformatted'                   ,&
        !             &     recordtype='stream'                  ,&
        !             &     action='write'                       ,&
        !             &     convert='LITTLE_ENDIAN'                 ,&
        !             &     access='sequential'                   )

        open(unit=101,file=trim(namda)//'.ensi.geo',&
             &     form='unformatted'                   ,&
             &     action='write'                       ,&
             &     convert='LITTLE_ENDIAN'                 ,&
             &     access='sequential'                   )

        if( kfl_bound == 1 ) then
           call ensmsh_bin(&
                kfl_bound,mnode,mnodb,npoin_2,nelem_2,nboun_2,100_ip,&
                101_ip,lexis,ltype,lnods,lbxis,ltypb,lnodb,lelch,coord,&
                ndime,namda,kfl_markm,npart_par,lsubd,nelem_par)
        else
           call ensmsh_bin(&
                kfl_bound,mnode,mnodb,npoin_2,nelem_2,nboun_2,100_ip,&
                101_ip,lexis,ltype,lnods,dummi,dummi,dummi,lelch,coord,&
                ndime,namda,kfl_markm,npart_par,lsubd,nelem_par)
        end if
     end if

     !close(100)
     close(101)

     !------------------------------
     !flag for writing the mesh of the filter with ensignt
     filt_msh_ensi=0
     !------------------------------

  else if ( forma == 'zfem' ) then
     !
     ! ZFEM de coco (ONLY FOR WEDGE ELEMENTS)
     !
     tipoe_ens= -1
     open(unit=100,file=trim(namda)//'.points.zfem',form='formatted')
     write(100,'(a)') 'CONTINUOUS'
     write(100,'(a)') 'NODS'
     write(100,'(a)') 'POINTS'
     write(100,    *) npoin_2
     do ipoin = 1,npoin_2
        write(100,'(3(1x,e13.6))') coord(1,ipoin),coord(2,ipoin),coord(3,ipoin)
     end do
     close(100)
     open(unit=100,file=trim(namda)//'.elems.zfem',form='formatted')
     write(100,'(a)') 'ELEMENTS'
     write(100,'(a)') 'PRISM'
     write(100,'(a)') 'ELEMENTS'
     write(100,    *) nelem_2
     do ielem = 1,nelem_2
        write(100,'(20(1x,i7))') (lnods(inode,ielem),inode=1,nnode(ltype(ielem)))
     end do
     write(100,'(a)') 'END'
     close(100)

  end if

  !----------------------------------------------------------------------
  !
  ! STL
  !
  !----------------------------------------------------------------------

  if( kfl_stlbo >= 1 .and. ndime == 3 .and. kfl_bound==1) then
     if (.not. associated(gevec)) allocate(gevec(3,1)) ! Temporary because wristl needs something
     call wristl(kfl_stlbo,namda,-1_ip,200,coord,0,gevec,npoin_total,lnodb,mnodb,ltypb,nboun_total)
     deallocate(gevec)
  end if

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  if( ndime == 2 ) then
     iesta = 10
     iesto = 29
  else if( ndime == 3 ) then
     iesta = 30
     iesto = 50
  end if
  !----------------------------------------------------------------------
  !
  ! Results
  !
  !----------------------------------------------------------------------

111 continue

  if (kfl_resta == 1) then
!!!!!  RESTART:

     ! read-write files
     allocate(gesca(npoin_total))
     allocate(gevec(ndofn,npoin_total))
     do iwopo= 1,kwopo
        if (kfofi == 1) then  ! ascii

           open(unit=15,file=trim(wonam(iwopo)),form='unformatted')
           read(15,*) ndofn
           if (ndofn == 1) then
              do ipoin= 1,npoin_total
                 read(15,*) gesca(ipoin)
              end do
           else
              do idofn= 1,ndofn
                 do ipoin= 1,npoin_total
                    read(15,*) gevec(idofn,ipoin)
                 end do
              end do
           end if

           close(15)

        else     ! binary

        end if

     end do
     deallocate(gesca)
     deallocate(gevec)

  end if

  if( forma == 'gid' ) then
     open(unit=11,file=trim(namda)//'.post.res',form='formatted')
     write(11,1) 'GiD Post Results File 1.0'
     write(11,1) ' '

  end if

  do

     read(49,'(a)',end=999) filna
     kauxi = len(trim(namda))+1
     if(filna(kauxi:kauxi)=='/') filna = filna(kauxi+1:len(filna))
     if( kfl_onlys /= 0 ) then
        kauxi = index(filna,'-',BACK = .TRUE.)
        read(filna(kauxi+1:kauxi+8),*) jstep
        if ( (kfl_ppste(jstep)/=0) .or. (jstep>max_onlys)) then
           kauxi = 1
        else
           kauxi = 0
        end if
     else
        kauxi = 1
     end if
     if( ( trim(filna) /= '' ) .and. ( kauxi == 1) )then
        if( kfl_conve == 0 ) then
           open(unit=10,file=trim(filna),form='unformatted')
        else
           open(unit=10,file=trim(filna),form='unformatted',convert='BIG_ENDIAN')
        end if
        ii = 10
        call livinf(1_ip,'CONVERT RESULT FILE: '//trim(filna),0_ip)
        call reahed(ii,npart_par,wwww8,iiiii,rrrrr)
        wwwww(1) = wwww8(1)(1:5)
        wwwww(2) = wwww8(2)(1:5)
        wwwww(3) = wwww8(3)(1:5)
        wwwww(4) = wwww8(4)(1:5)
        wwwww(5) = wwww8(5)(1:5)
        wwwww(6) = wwww8(6)(1:5)
        wwwww(7) = wwww8(7)(1:5)
        wwwww(8) = wwww8(8)(1:5)
        wwwww(9) = wwww8(9)(1:5)

        pdime       = int (iiiii(1),ip)
        npoin_total = int (iiiii(2),ip)
        ittim       = int (iiiii(4),ip)
        rttim       = real(rrrrr(1),rp)

        if( wwwww(4) == 'SCALA' ) then

           !-------------------------------------------------------------
           !
           ! Scalar
           !
           !-------------------------------------------------------------

           if( forma == 'gid' ) then
              if( wwwww(5) == 'NELEM' .or. wwwww(5) == 'NBOUN' ) then
                 do ielty=iesta,iesto
                    if(lexis(ielty)/=0) then
                       if(ndime==2) then
                          if(     nnode(ielty)==2) then
                             elemt='Linear'
                          elseif(nnode(ielty)==3.or.nnode(ielty)==6.or.nnode(ielty)==7) then
                             elemt='Triangle'
                          else
                             elemt='Quadrilateral'
                          end if
                       else
                          if(nnode(ielty)==4.or.nnode(ielty)==10) then
                             elemt='Tetrahedra'
                          else if(nnode(ielty)==8.or.nnode(ielty)==20.or.nnode(ielty)==27) then
                             elemt='Hexahedra'
                          else if(nnode(ielty)==6.or.nnode(ielty)==15) then
                             elemt='Prism'
                          end if 
                       end if

                       write(11,81) 'GaussPoints '//'GP'//' Elemtype '//trim(elemt)
                       write(11,85)  'Number of Gauss Points: ',1
                       write(11,81) 'Natural Coordinates: Internal'
                       write(11,81) 'End GaussPoints'
                    end if
                 end do
                 write(11,82) wwwww(3),'ALYA',rrrrr(1),'Scalar','GP'
              else
                 write(11,2)  wwwww(3),'ALYA',rrrrr(1),'Scalar'
              end if
              write(11,3) wwwww(3)
              write(11,1) 'Values'
           else if( forma == 'vu' .or. forma == 'filtre_vu'.or. forma == 'ensight' .or. forma=='zfem') then
              mpoin=0
              allocate(gisc3(npoin_total))
              allocate(gesc3(npoin_total))
              do ipoin = 1,npoin_total
                 gisc3(ipoin) = 0_ip
                 gesc3(ipoin) = 0_rp
              end do
           end if

           if( wwwww(6) == 'INTEG' ) then

              !if( npoin /= 0 ) then

              if( wwwww(9) == 'FILTE' ) then
                 !
                 ! Int(:) with filter
                 !
                 ipois = 0
                 if( kfl_elimi == 0 ) then
                    do ipart = 1,npart_par
                       if( lsubd(ipart) == 1 ) then
                          read(ii) npoin
                          if( npoin > 0 ) then
                             allocate( gisca(npoin) )
                             allocate( gesca(npoin) )
                             read(ii) ( gisca(ipoin), ipoin=1,npoin )
                             read(ii) ( gesca(ipoin), ipoin=1,npoin )
                             if( forma == 'gid' ) then
                                do ipoin = 1,npoin
                                   write(11,4) ipois+gisca(ipoin),gesca(ipoin)
                                end do
                             else if( forma == 'vu' .or. forma == 'filtre_vu'.or. forma == 'ensight' .or. forma == 'zfem')then
                                do ipoin = 1,npoin
                                   gisc3(ipois+gisca(ipoin))=ipois+gisca(ipoin)
                                   gesc3(ipois+gisca(ipoin))=gesca(ipoin)
                                   !write(*,*)'gisc3,gesc3',gisc3,gesc3
                                end do
                             end if
                             deallocate( gisca )
                             deallocate( gesca )
                          end if
                          ipois = ipois + npoin_par(ipart)
                       else
                          read(ii) npoin
                          if( npoin > 0 ) then
                             read(ii) ( dummi,ipoin=1,npoin)
                             read(ii) ( dummr,ipoin=1,pdime*npoin)
                          end if
                       end if
                    end do

                    if ( forma == 'vu' ) then
                       call vu_res(&
                            ittim,npoin_total,nelem_total,gesc3,dummi,102_ip,103_ip,lexis,&
                            lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                            wwwww(7),rttim,kfl_markm,kfl_multi,kfl_bound)

                    else if (forma == 'filtre_vu' ) then
                       call vu_filter(&
                            kfl_bound,mnode,mnodb,npoin_total,nelem_total,nboun,104_ip,105_ip,&
                            lexis,ltype,lnods,dummi,dummi,dummi,dummi,coord,&
                            ndime,namda,kfl_markm,dummi,dummi,dummi,&
                            ittim,wwwww(3),wwwww(4),wwwww(6),wwwww(7),&
                            rttim,kfl_multi,gesc3,gisc3)

                    else if (forma == 'ensight' ) then
                       write(6,*) 'NOT CODED'
                    end if
                    write(6,*) 'format=',forma
                    write(6,*) 'NOT CODED'
                 end if

              else
                 !
                 ! Int(:) without filter
                 !
                 allocate( gisca(npoin_2) )
                 ipois = 0
                 if( kfl_elimi == 0 ) then
                    do ipart = 1,npart_par
                       if( lsubd(ipart) == 1 ) then
                          read(ii) npoin
                          read(ii) (gisca(ipoin),ipoin=1+ipois,npoin+ipois)
                          ipois = ipois + npoin_par(ipart)
                       else
                          read(ii) npoin
                          read(ii) (dummi,ipoin=1,npoin)
                       end if
                    end do
                 else

                    if( wwwww(5) == 'NPOIN' ) then
                       do ipart = 1,npart_par
                          if( lsubd(ipart) == 1 ) then
                             read(ii) npoin
                             allocate( gisc2(npoin_par(ipart)) )
                             read(ii) ( gisc2(ipoin),ipoin=1,npoin)
                             do ipoin = 1,npoin
                                jpoin = lninv(ipoin+ipois)
                                gisca(jpoin) = real(gisc2(ipoin),rp)
                             end do
                             deallocate( gisc2 )
                             ipois = ipois + npoin_par(ipart)
                          else
                             read(ii) npoin
                             read(ii) ( dummi,ipoin=1,npoin)
                          end if
                       end do
                    else if( wwwww(5) == 'NELEM' ) then
                       call runend('INT NELEM NOT CODED')
                    else if( wwwww(5) == 'NBOUN' ) then
                       call runend('INT NBOUN NOT CODED')
                    end if

                 end if
                 if( forma == 'gid' ) then
                    do ipoin = 1,npoin_2
                       write(11,6) ipoin,gisca(ipoin)
                    end do
                 else if( forma == 'vu' ) then
                    call vu_res(&
                         ittim,npoin_2,nelem_2,dummr,gisca,102_ip,103_ip,lexis,&
                         lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                         wwwww(7),rttim,kfl_markm,kfl_multi,kfl_bound)
                 else if( forma == 'zfem' ) then
                    call zfemres(&
                         ittim,npoin_2,nelem_2,dummr,gisca,lexis,&
                         lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                         wwwww(7),wwwww(9),rttim,kfl_markm,kfl_multi)
                 else if( forma == 'ensight' ) then
                    if ( kfl_ensbi == 0 ) then
                       call ensres(&
                            ittim,npoin_2,nelem_2,dummr,gisca,lexis,&
                            lbxis,ltype,leinv,ndime,namda,wwwww(3),wwwww(4),wwwww(5),wwwww(6),&
                            wwwww(7),wwwww(9),rttim,kfl_markm,kfl_multi,kfl_field)
                    else
                       call ensres_bin(&
                            ittim,npoin_2,nelem_2,dummr,gisca,lexis,&
                            lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                            wwwww(7),rttim,kfl_markm,kfl_multi)

                    end if
                 else if( forma == 'txt' ) then
                    call txtres(&
                         ittim,npoin_2,nelem_2,dummr,gisca,lexis,&
                         lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                         wwwww(7),rttim,kfl_markm,kfl_multi)
                 end if
                 deallocate(gisca)
              end if
              !end if

           else
              !
              ! SCALAR REAL
              !

              !if( npoin /= 0 ) then

              if( wwwww(9) == 'FILTE' ) then
                 !
                 ! real(:) with filter
                 !
                 ipois = 0
                 if( kfl_elimi == 0 ) then
                    do ipart = 1,npart_par
                       if( lsubd(ipart) == 1 ) then
                          read(ii) npoin
                          if( npoin > 0 ) then
                             allocate( gisca(npoin) )
                             allocate( gesca(npoin) )
                             read(ii) ( gisca(ipoin), ipoin=1,npoin )
                             read(ii) ( gesca(ipoin), ipoin=1,npoin )
                             if( forma == 'gid' ) then
                                do ipoin = 1,npoin
                                   write(11,4) ipois+gisca(ipoin),gesca(ipoin)
                                end do
                             else if( forma == 'vu' .or. forma == 'filtre_vu' .or. forma == 'ensight') then
                                do ipoin = 1,npoin
                                   gisc3(ipois+gisca(ipoin))=ipois+gisca(ipoin)
                                   gesc3(ipois+gisca(ipoin))=gesca(ipoin)
                                   !write(*,*)'gisc3,gesc3',gisc3(ipoin),gesc3(ipoin)
                                end do
                             end if
                             deallocate( gisca )
                             deallocate( gesca )
                          end if
                          ipois = ipois + npoin_par(ipart)
                       else
                          read(ii) npoin
                          if( npoin > 0 ) then
                             read(ii) ( dummi,ipoin=1,npoin)
                             read(ii) ( dummr,ipoin=1,pdime*npoin)
                          end if
                       end if
                    end do
                    if ( forma == 'ensight') then
                       if(filt_msh_ensi==0)then
                          open(unit=113,file=trim(namda)//'-filter.ensi.case',form='formatted')
                          open(unit=114,file=trim(namda)//'-filter.ensi.geo',form='formatted')

                          call ensmsh_filter(&
                               kfl_bound,mnode,mnodb,npoin_total,nelem_total,nboun,113_ip,&
                               114_ip,lexis,ltype,lnods,lbxis,ltypb,lnodb,lelch,coord,&
                               ndime,namda,kfl_markm,npart_par,lsubd,nelem_par,&
                               gesc3,gisc3,dummr)

                          close(114)
                          !close(113)
                          filt_msh_ensi=1
                       end if
                       call ensres_filter(&
                            ittim,npoin_total,nelem_total,gesc3,dummi,lexis,&
                            lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                            wwwww(7),wwwww(9),rttim,kfl_markm,kfl_multi,&
                            gisc3,mnode,lnods)

                    end if
                    !
                    ! binaire
                    !
                    if ( forma == 'vu' .or. forma == 'filtre_vu' ) then
                       mpoin_2=0
                       do ipoin=1,npoin_total
                          if (gisc3(ipoin) /= 0) then
                             mpoin_2=mpoin_2+1
                             gesc3(ipoin)=1_rp
                          else
                             gesc3(ipoin)=0_rp
                          end if
                       end do
                    end if
                    if ( forma == 'vu' ) then
                       call vu_res(&
                            ittim,npoin_total,nelem_total,gesc3,dummi,102_ip,103_ip,lexis,&
                            lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                            wwwww(7),rttim,kfl_markm,kfl_multi,kfl_bound)
                    else if ( forma == 'filtre_vu' ) then
                       call vu_filter(&
                            kfl_bound,mnode,mnodb,npoin_total,nelem_total,nboun_total,104_ip,105_ip,&
                            lexis,ltype,lnods,dummi,ltypb,lnodb,dummi,coord,&
                            ndime,namda,kfl_markm,dummi,dummi,dummi,&
                            ittim,wwwww(3),wwwww(4),wwwww(6),wwwww(7),&
                            rttim,kfl_multi,gesc3,gisc3)
                    end if

                 else
                    print*,'format=',forma
                    write(6,*) 'NOT CODED'
                 end if

              else
                 !
                 ! real(:) without filter
                 !
                 if( wwwww(5) == 'NELEM' ) then

                    allocate( gesca(nelem_total) )
                    ieles = 0
                    do ipart = 1,npart_par
                       if( lsubd(ipart) == 1 ) then
                          read(ii) nelem
                          read(ii) (gesca(ielem),ielem=1+ieles,nelem+ieles)
                          ieles = ieles + nelem_par(ipart)
                       else
                          read(ii) nelem
                          read(ii) (dummr,ielem=1,nelem)
                       end if
                    end do

                    if( forma == 'gid' ) then
                       if( kfl_order == 1 ) call maths_heapsort_real(2_ip,1_ip,nelem_total,leinv_cpy,gesca,SAVING=.true.)
                       do ielem = 1,nelem_total
                          write(11,4) leinv(ielem),gesca(ielem)
                       end do

                    else if( forma == 'alya' .and. wwwww(3) == 'MATER' ) then
                       do ielem = 1,nelem_total
                          write(97,*) ielem,int(gesca(ielem))
                       end do
                       stop
                    else if( forma == 'ensight') then
                       if ( kfl_ensbi == 0 ) then
                          call ensres(&
                               ittim,npoin,nelem_total,gesca,dummi,lexis,&
                               lbxis,ltype,leinv,ndime,namda,wwwww(3),wwwww(4),wwwww(5),wwwww(6),&
                               wwwww(7),wwwww(9),rttim,kfl_markm,kfl_multi,kfl_field)
                       end if
                    end if

                    deallocate(gesca)

                  else if( wwwww(5) == 'NBOUN' ) then

                    allocate( gesca(nboun_total) )
                    ieles = 0
                    do ipart = 1,npart_par
                       if( lsubd(ipart) == 1 ) then
                          read(ii) nboun
                          read(ii) (gesca(ielem),ielem=1+ieles,nboun+ieles)
                          ieles = ieles + nboun_par(ipart)
                       else
                          read(ii) nboun
                          read(ii) (dummr,ielem=1,nboun)
                       end if
                    end do

                    if( forma == 'gid' ) then

                       !allocate(lbinv_cpy(nboun_total))
                       !lbinv_cpy = lbinv
                       !call maths_heapsort_real(2_ip,1_ip,nboun_total,lbinv_cpy,gesca)
                       if( kfl_order == 1 ) call maths_heapsort_real(2_ip,1_ip,nboun_total,lbinv_cpy,gesca,SAVING=.true.)                       
                       do iboun = 1,nboun_total 
                          write(11,4) lbinv(iboun),gesca(iboun)
                       end do
                       !deallocate(lbinv_cpy)

                    else
                       print*,'NBOUN REAL SCALAR MOT CODED'
                       stop
                    end if

                    deallocate(gesca)

                else

                    allocate( gesca(npoin_2) )
                    ipois = 0
                    if( kfl_elimi == 0 ) then
                       do ipart = 1,npart_par
                          if( lsubd(ipart) == 1 ) then
                             read(ii) npoin
                             read(ii) (gesca(ipoin),ipoin=1+ipois,npoin+ipois)
                             ipois = ipois + npoin_par(ipart)
                          else
                             read(ii) npoin
                             read(ii) (dummr,ipoin=1,npoin)
                          end if
                       end do
                    else
                       do ipoin = 1,npoin_2
                          gesca(ipoin) = 0.0_rp
                       end do
                       do ipart = 1,npart_par
                          if( lsubd(ipart) == 1 ) then
                             read(ii) npoin
                             allocate( gesc2(npoin_par(ipart)) )
                             read(ii) ( gesc2(ipoin),ipoin=1,npoin)
                             do ipoin = 1,npoin
                                jpoin = lninv(ipoin+ipois)
                                gesca(jpoin) = gesc2(ipoin)
                             end do
                             deallocate( gesc2 )
                             ipois = ipois + npoin_par(ipart)
                          else
                             read(ii) npoin
                             read(ii) ( dummr,ipoin=1,npoin)
                          end if
                       end do
                    end if

                    if( forma == 'gid' ) then
                       do ipoin = 1,npoin_2
                          write(11,4) ipoin,gesca(ipoin)
                       end do
                    else if( forma == 'vu' ) then
                       call vu_res(&
                            ittim,npoin_2,nelem_2,gesca,dummi,102_ip,103_ip,lexis,&
                            lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                            wwwww(7),rttim,kfl_markm,kfl_multi,kfl_bound)
                    else if( forma == 'zfem' ) then
                       call zfemres(&
                            ittim,npoin_2,nelem_2,gesca,dummi,lexis,&
                            lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                            wwwww(7),wwwww(9),rttim,kfl_markm,kfl_multi)
                    else if( forma == 'ensight' ) then
                       if ( kfl_ensbi == 0 ) then
                          call ensres(&
                               ittim,npoin_2,nelem_2,gesca,dummi,lexis,&
                               lbxis,ltype,leinv,ndime,namda,wwwww(3),wwwww(4),wwwww(5),wwwww(6),&
                               wwwww(7),wwwww(9),rttim,kfl_markm,kfl_multi,kfl_field)
                       else
                          call ensres_bin(&
                               ittim,npoin_2,nelem_2,gesca,dummi,lexis,&
                               lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                               wwwww(7),rttim,kfl_markm,kfl_multi)
                       end if
                    else if( forma == 'txt' ) then
                       call txtres(&
                            ittim,npoin_2,nelem_2,gesca,dummi,lexis,&
                            lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                            wwwww(7),rttim,kfl_markm,kfl_multi)
                    end if
                    deallocate(gesca)
                 end if
              end if
              !end if
           end if

           if( forma == 'gid' ) then
              write(11,1) 'End values'
              write(11,1) ''
           else if ( forma == 'vu' .or. forma == 'filtre_vu'.or.forma=='ensight') then
              deallocate(gisc3)
              deallocate(gesc3)
           end if

        else if( wwwww(4) == 'VECTO' ) then

           !-------------------------------------------------------------
           !
           ! Vector
           !
           !-------------------------------------------------------------

           if( forma == 'gid' ) then
              if( wwwww(5) == 'NELEM' .or. wwwww(5) == 'NBOUN' ) then
                 do ielty=iesta,iesto
                    if(lexis(ielty)/=0) then
                       if(ndime==2) then
                          if(     nnode(ielty)==2) then
                             elemt='Linear'
                          else if(nnode(ielty)==3.or.nnode(ielty)==6.or.nnode(ielty)==7) then
                             elemt='Triangle'
                          else
                             elemt='Quadrilateral' 
                          end if
                       else
                          if(nnode(ielty)==4.or.nnode(ielty)==10) then
                             elemt='Tetrahedra'
                          else if(nnode(ielty)==8.or.nnode(ielty)==20.or.nnode(ielty)==27) then
                             elemt='Hexahedra'
                          else if(nnode(ielty)==6.or.nnode(ielty)==15) then
                             elemt='Prism'
                          end if
                       end if

                       write(11,81) 'GaussPoints '//'GP'//' Elemtype '//trim(elemt)
                       write(11,85)  'Number of Gauss Points: ',1
                       write(11,81) 'Natural Coordinates: Internal'
                       write(11,81) 'End GaussPoints'
                    end if
                 end do
                 write(11,82) wwwww(3),'ALYA',rrrrr(1),'Vector','GP'
              else
                 write(11,2) wwwww(3),'ALYA',rrrrr(1),'Vector'
              end if
              write(11,3) wwwww(3)//'_X,'//wwwww(3)//'_Y,'//wwwww(3)//'_Z'
              write(11,1) 'Values'

           else if( forma == 'vu' .or. forma == 'filtre_vu'.or. forma == 'ensight') then

              mpoin=0
              allocate(gisc3(npoin_total))
              allocate(geve3(pdime,npoin_total))

              do ipoin = 1,npoin_total
                 gisc3(ipoin) = 0_ip
                 do idime=1,pdime
                    geve3(idime,ipoin) = 0_rp
                 end do
              end do

           end if

           !if( npoin /= 0 ) then

           if( wwwww(9) == 'FILTE' ) then
              !
              ! real(:,:) with filter
              !
              ipois = 0
              if( kfl_elimi == 0 ) then
                 do ipart = 1,npart_par
                    if( lsubd(ipart) == 1 ) then
                       read(ii) npoin
                       if( npoin > 0 ) then
                          allocate( gisca(npoin) )
                          allocate( gevec(pdime,npoin) )
                          read(ii) ( gisca(ipoin), ipoin=1,npoin )
                          read(ii) ( (gevec(idime,ipoin),idime=1,pdime),ipoin=1,npoin)
                          if( forma == 'gid' ) then
                             do ipoin = 1,npoin
                                write(11,4) ipois+gisca(ipoin),(gevec(idime,ipoin),idime=1,pdime)
                             end do
                          else if( forma == 'vu' .or. forma == 'filtre_vu'.or. forma == 'ensight') then
                             do ipoin = 1,npoin
                                gisc3(ipois+gisca(ipoin))=ipois+gisca(ipoin)
                                do idime = 1,pdime
                                   geve3(idime,ipois+gisca(ipoin))=gevec(idime,ipoin)
                                end do
                             end do
                          end if
                          deallocate( gisca )
                          deallocate( gevec )
                       end if
                       ipois = ipois + npoin_par(ipart)
                    else
                       read(ii) npoin
                       if( npoin > 0 ) then
                          read(ii) ( dummi,ipoin=1,npoin)
                          read(ii) ( dummr,ipoin=1,pdime*npoin)
                       end if
                    end if
                 end do

                 if( forma == 'vu' ) then
                    call vu_res(&
                         ittim,npoin_2,nelem_2,geve3,dummi,102_ip,103_ip,lexis,&
                         lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                         wwwww(7),rttim,kfl_markm,kfl_multi,kfl_bound)
                 else if( forma == 'filtre_vu' ) then
                    call vu_filter(&
                         kfl_bound,mnode,mnodb,npoin_total,nelem_total,nboun,104_ip,105_ip,&
                         lexis,ltype,lnods,dummi,dummi,dummi,dummi,coord,&
                         ndime,namda,kfl_markm,dummi,dummi,dummi,&
                         ittim,wwwww(3),wwwww(4),wwwww(6),wwwww(7),&
                         rttim,kfl_multi,geve3,gisc3)
                 else if( forma == 'ensight' ) then
                    if(filt_msh_ensi==0)then
                       open(unit=113,file=trim(namda)//'-filter.ensi.case',form='formatted')
                       open(unit=114,file=trim(namda)//'-filter.ensi.geo',form='formatted')

                       call ensmsh_filter(&
                            kfl_bound,mnode,mnodb,npoin_total,nelem_total,nboun,113_ip,&
                            114_ip,lexis,ltype,lnods,lbxis,ltypb,lnodb,lelch,coord,&
                            ndime,namda,kfl_markm,npart_par,lsubd,nelem_par,&
                            gesc3,gisc3,dummr)

                       close(114)
                       !close(113)
                       filt_msh_ensi=1
                    end if
                    call ensres_filter(&
                         ittim,npoin_total,nelem_total,geve3,dummi,lexis,&
                         lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                         wwwww(7),wwwww(9),rttim,kfl_markm,kfl_multi,&
                         gisc3,mnode,lnods)
                 end if
              else
                 print*,'format=',forma
                 write(6,*) 'NOT CODED'
              end if

           else
              !
              ! real(:,:) without filter
              !
              if( wwwww(5) == 'NELEM' ) then

                 allocate( gevec(pdime,nelem_total) )
                 gevec = 0.0_rp
                 ieles = 0
                 do ipart = 1,npart_par
                    if(lsubd(ipart) == 1 ) then
                       read(ii) nelem
                       read(ii) ( (gevec(idime,ielem),idime=1,pdime),ielem=ieles+1,nelem+ieles)
                       ieles = ieles + nelem_par(ipart)
                    else
                       read(ii) nelem
                       read(ii) ( dummr,ielem=1,pdime*nelem)
                    end if
                 end do

                 if( forma == 'gid' ) then

                    if( kfl_order == 1 ) call maths_heapsort_real(2_ip,pdime,nelem_total,leinv_cpy,gevec,SAVING=.true.)
                    do ielem = 1,nelem_total
                       write(11,4) leinv(ielem),(gevec(idime,ielem),idime=1,pdime)
                    end do

                 else if (forma == 'ensight') then

                    if ( kfl_ensbi == 0 ) then
                       call ensres(&
                            ittim,npoin_2,nelem_2,gevec,dummi,lexis,                    &
                            lbxis,ltype,leinv,ndime,namda,wwwww(3),wwwww(4),wwwww(5),wwwww(6),&
                            wwwww(7),wwwww(9),rttim,kfl_markm,kfl_multi,kfl_field)
                    else
                       print*,'NOT CODED'
                    end if
                 else
                    print*,'FORMAT NOT CODED'
                 end if

              else if( wwwww(5) == 'NBOUN' ) then

                 if( nboun_total == 0 ) print*,'ACTIVATE BOUNDARY OUTPUT'
                 allocate( gevec(pdime,nboun_total) )
                 gevec = 0.0_rp
                 ieles = 0
                 do ipart = 1,npart_par
                    if(lsubd(ipart) == 1 ) then
                       read(ii) nboun
                       read(ii) ( (gevec(idime,iboun),idime=1,pdime),iboun=ieles+1,nboun+ieles)
                       ieles = ieles + nboun_par(ipart)
                    else
                       read(ii) nboun
                       read(ii) ( dummr,iboun=1,pdime*nboun)
                    end if
                 end do

                 if( forma == 'gid' ) then

                    !allocate(lbinv_cpy(nboun_total))
                    !lbinv_cpy = lbinv
                    !call maths_heapsort_real(2_ip,pdime,nboun_total,lbinv_cpy,gevec)
                    if( kfl_order == 1 ) call maths_heapsort_real(2_ip,pdime,nboun_total,lbinv_cpy,gevec,SAVING=.true.)                    
                    do iboun = 1,nboun_total 
                       write(11,4) lbinv(iboun),(gevec(idime,iboun),idime=1,pdime)
                    end do
                    !deallocate(lbinv_cpy)
                    
                 else
                    
                    print*,'FORMAT NOT CODED'
                    
                 end if

              else
                 if( forma == 'alya' .and. wwwww(3) == 'CODNO' ) then
                    allocate( gevec(pdime,npoin_2) )
                    read(ii) npoin_2
                    read(ii) ( dummr,ipoin=1,pdime*npoin_2)
                    do ipoin = 1,npoin_2
                       if( abs(gevec(1,ipoin)-3.0_8) < 1.0e-3 ) then
                          write(97,*) ipoin,2,1
                       else if( abs(gevec(1,ipoin)-4.0_8) < 1.0e-3 ) then
                          write(97,*) ipoin,1,2
                       end if
                    end do
                    print*,'caca'
                    stop
                 end if
                 allocate( gevec(pdime,npoin_2) )
                 ipois = 0
                 if( kfl_elimi == 0 ) then
                    do ipart = 1,npart_par
                       if( lsubd(ipart) == 1 ) then
                          read(ii) npoin
                          read(ii) ( (gevec(idime,ipoin),idime=1,pdime),ipoin=1+ipois,npoin+ipois)
                          ipois = ipois + npoin_par(ipart)
                       else
                          read(ii) npoin
                          read(ii) ( dummr,ipoin=1,pdime*npoin)
                       end if
                    end do
                 else
                    do ipoin = 1,npoin_2
                       do idime = 1,ndime
                          gevec(idime,ipoin) = 0.0_rp
                       end do
                    end do
                    do ipart = 1,npart_par
                       if( lsubd(ipart) == 1 ) then
                          read(ii) npoin
                          allocate( geve2(ndime,npoin_par(ipart)) )
                          read(ii) ( (geve2(idime,ipoin),idime=1,pdime),ipoin=1,npoin)
                          do ipoin = 1,npoin
                             jpoin = lninv(ipoin+ipois)
                             do idime = 1,ndime
                                gevec(idime,jpoin) = geve2(idime,ipoin)
                             end do
                          end do
                          deallocate( geve2 )
                          ipois = ipois + npoin_par(ipart)
                       else
                          read(ii) npoin
                          read(ii) ( dummr,ipoin=1,pdime*npoin)
                       end if
                    end do
                 end if
                 if( forma == 'gid' ) then
                    do ipoin = 1,npoin_2
                       write(11,4) ipoin,(gevec(idime,ipoin),idime=1,ndime)
                    end do

                 else if( forma == 'zfem' ) then
                    call zfemres(&
                         ittim,npoin_2,nelem_2,gevec,dummi,lexis,&
                         lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                         wwwww(7),wwwww(9),rttim,kfl_markm,kfl_multi)

                 else if( forma == 'vu' ) then
                    call vu_res(&
                         ittim,npoin_2,nelem_2,gevec,dummi,102_ip,103_ip,lexis,&
                         lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                         wwwww(7),rttim,kfl_markm,kfl_multi,kfl_bound)
                 else if( forma == 'ensight' ) then
                    if ( kfl_ensbi == 0 ) then
                       call ensres(&
                            ittim,npoin_2,nelem_2,gevec,dummi,lexis,&
                            lbxis,ltype,leinv,ndime,namda,wwwww(3),wwwww(4),wwwww(5),wwwww(6),&
                            wwwww(7),wwwww(9),rttim,kfl_markm,kfl_multi,kfl_field)
                    else
                       call ensres_bin(&
                            ittim,npoin_2,nelem_2,gevec,dummi,lexis,&
                            lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                            wwwww(7),rttim,kfl_markm,kfl_multi)
                    end if
                 else if( forma == 'txt' ) then
                    call txtres(&
                         ittim,npoin_2,nelem_2,gevec,dummi,lexis,&
                         lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
                         wwwww(7),rttim,kfl_markm,kfl_multi)
                 end if

              end if
              !----------------------------------------------------------------------
              !
              ! STL with displacements
              !
              !----------------------------------------------------------------------
              if( wwwww(3)=='DISPL' .and. kfl_stlbo >= 1 .and. ndime == 3) then
                 call wristl(kfl_stlbo,namda,ittim,200,coord,1_ip,gevec,npoin_2,lnodb,mnodb,ltypb,nboun_total)
              end if

              deallocate(gevec)

           end if

           if( forma == 'gid' ) then
              write(11,1) 'End values'
              write(11,1) ''
           else if ( forma == 'vu' .or. forma == 'filtre_vu' .or. forma == 'ensight' ) then
              deallocate(gisc3)
              deallocate(geve3)
           end if

        else if( wwwww(4) == 'R3P  ' .or. wwwww(4) == 'R3PVE' ) then

           !-------------------------------------------------------------
           !
           ! R3P
           !
           !-------------------------------------------------------------
           !
           ! Read information and results
           !
           read(ii) ( ngaus(ielty),ielty=iesta,iesto)
           mgaus = 0
           do ielty = iesta,iesto
              mgaus = max(mgaus,ngaus(ielty))
           end do
           !
           ! Header for Gauss Points
           !
           if (ittim == 1_ip) call gidres_he(iesta,iesto,ndime,ngaus(:),lexis(:),flag_coh)

           ieles = 0
           if( kfl_elimi == 0 ) then
              allocate( geve4(ndime,mgaus,nelem_2) )
              ieles = 0
              do ipart = 1,npart_par
                 if( lsubd(ipart) == 1 ) then
                    read(ii) nelem
                    read(ii) dummi
                    read(ii) ( ( (geve4(idim1,idim2,ielem),idim1=1,pdime),&
                         idim2=1,ngaus(abs(ltype(ielem)))),ielem=ieles+1,ieles+nelem)
                 else
                    read(ii) nelem
                    read(ii) dummi
                    read(ii) ( ( (dummr,idim1=1,pdime),idim2=1,ngaus(abs(ltype(ielem)))),ielem=ieles+1,ieles+nelem)
                 end if
                 ieles = ieles + nelem_par(ipart)
              end do
              !
              ! Write results
              !
              call gidres_gp(iesta,iesto,wwwww(3),rrrrr(1),ndime,pdime, &
                   nelem_2,mgaus,ngaus(:),lexis(:),ltype(:),leinv(:),lelch(:),flag_coh,geve4(:,:,:))

              deallocate( geve4 )

           else
              print*, 'R3P ELIMINATING BOUNDARY NODES NOT CODED'
              stop
           end if

        end if

        close(10)

     end if
  end do

999 continue

  if( forma == 'ensight' ) then
     if ( kfl_ensbi == 0 ) then
        call ensres(&
             -1_ip,npoin_2,nelem_2,dummr,dummi,lexis,&
             lbxis,ltype,leinv,ndime,namda,wwwww(3),wwwww(4),wwwww(5),wwwww(6),&
             wwwww(7),wwwww(9),rttim,kfl_markm,kfl_multi,kfl_field)
     else
        call ensres_bin(&
             -1_ip,npoin_2,nelem_2,dummr,dummi,lexis,&
             lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
             wwwww(7),rttim,kfl_markm,kfl_multi)
     end if
  else if( forma == 'txt' ) then
     call txtres(&
          -1_ip,npoin_2,nelem_2,dummr,dummi,lexis,&
          lbxis,ltype,pdime,namda,wwwww(3),wwwww(4),wwwww(6),&
          wwwww(7),rttim,kfl_markm,kfl_multi)
  end if

!  call livinf(1_ip,'KEEP READING (YES=1)?',0_ip)
!  read(5,'(i1)') ipoin
!  if( ipoin == 1 ) goto 111

  call livinf(0_ip,' ',0_ip)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  close(113)  ! file de ensi-filter-case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! GiD formats
  !
1 format(a)
2 format('Result ',a,' ',a,' ',e15.8,' ',a,' OnNodes')
3 format('ComponentNames ',a)
4 format(i11, 3(1x,e16.8E3))
5 format('ComponentNames ',a,a,a)
6 format(i11, 3(1x,i8))
81 format(a)
82 format('Result ',a,' ',a,' ',e15.8,' ',a,' OnGaussPoints ',a)
83 format('ComponentNames ',a)
84 format(i11, 3(1x,e16.8E3))
85 format(a,1x,i2)
86 format(3(1x,e16.8E3))
  !
  ! Femview formats
  !
10 format(1x,i4,a1,a6,e12.5,32x,i2,i5)
20 format(1x,i2,2x,a5,3x,3i5)
30 format(1x,i2,2x,a5,3x,2i5)
40 format(1x,i2,i5,e12.5)
50 format(1x,i2)
  !
  ! Ensight Gold formats
  !
100 format(a)
110 format(i11)
120 format(e16.8E3)
  !
  ! VU format
  !
205 format('TEXTE Time(" ',e13.6,'");')
200 format('FIELD<double> ',a,'("',a,'",',i7,',',i12,');')
201 format('FIELD<float>  ',a,'("',a,'",',i7,',',i12,');')
202 format('FIELD<int>    ',a,'("',a,'",',i7,',',i12,');')
210 format('SOLUTION Solution( ) =',/,'{')
220 format('   VARIABLE ',a,'( ',a,',',a,',',a,',',a,');')
221 format('   VARIABLE ',a,'( ',a,',',a,',',a,',',a,',',a,');')
230 format('};')
  !
  ! Gnuplot format
  !
600 format(a)
610 format(i11)
620 format(10(1x,e16.8E3))
  !
  ! Alya ASCII format
  !
800 format(a)
810 format(i11)
820 format(i11,10(1x,e16.8E3))

end program alya2pos
