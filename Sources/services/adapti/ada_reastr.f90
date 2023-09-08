subroutine ada_reastr(itask)
!-----------------------------------------------------------------------
!****f* adapti/Adapti
! NAME 
!    ada_reastr
! DESCRIPTION
!    This routine reads adapti strategy
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_inpout
  use      def_domain
  use      mod_memchk

  use      def_adapti
  use mod_ecoute, only :  ecoute

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: &
       istat,iread,iimmo,ielem,kelem,jelem,inode,ipoin,kpoin,jpoin,&
       idime,iauxi
  integer(ip)             :: lnolo(10),coolo(3)

  real(rp)                :: cauxi(3,2)

  if (itask == 1) then
     !
     ! Reading just the parameters
     !
     
     !
     ! Initializations (defaults)
     !
     kfl_remst_ada = 0                                    ! No remesh strategy, do nothing
     kfl_redgr_ada = 0                                    ! No regdreen strategy
     kfl_inter_ada = 0                                    ! No interpolation
     kfl_chpoi_ada = 0                                    ! No chekpoint
     kfl_dista_ada = 0                                    ! Do not compute distance function
     kfl_prpro_ada = 1                                    ! Just pre-process
     kfl_ougid_ada = 1                                    ! Output gid format
     nimmo_ada     = 0                                    ! No immersed objects
     npoii_ada     = 0
     nelei_ada     = 0
     nelto_ada     = 0
     npoto_ada     = 0
     ndimi_ada     = 0
     nnodi_ada     = 0    

     vloca_ada = 0.0_rp                                  ! Origin in local reference
     vorig_ada = 0.0_rp                                  ! Origin in background reference
     vscal_ada = 1.0_rp                                  ! Scale factors in each direction
     vroti_ada = 0.0_rp                                  ! Rotations around X,Y,Z axes in sexa
     
     !
     ! Read and write files
     !
     lispa = 0
     lisda = lun_pdata_ada                                ! Reading file
     lisre = lun_outpu_ada                                ! Writing file

     if(kfl_paral<=0) then
        !
        ! Reach the section
        !
        rewind(lisda)
        call ecoute('ada_reastr')
        do while(words(1)/='ADAPT')
           call ecoute('ada_reastr')
        end do
        
        !
        ! Begin to read data
        !
        do while(words(1)/='ENDAD')     
           call ecoute('ada_reastr')        
           
           if(words(1) .eq.'IMMER') then
              iread= 0
              iimmo= 0
              do while(words(1)/='ENDIM')
                 call ecoute('ada_reastr')
                 if (words(1) .eq.'TOTAL') then
                    nimmo_ada=getint('TOTAL',1,'#Number of immersed objects')                 
                 else if (words(1) .eq.'DIMEN') then
                    iimmo = iimmo + 1
                    do while(words(1)/='ENDDI')
                       call ecoute('ada_reastr')
                       if (words(1) .eq.'NODAL') npoii_ada(iimmo)= getint('NODAL',1,'#Nodal points')
                       if (words(1) .eq.'ELEME') nelei_ada(iimmo)= getint('ELEME',1,'#Elements')
                       if (words(1) .eq.'SPACE') ndimi_ada       = getint('SPACE',1,'#Spaced dimension')
                       if (words(1) .eq.'NODES') nnodi_ada       = getint('NODES',1,'#Nodes per element')
                    end do
                 else if (words(1) .eq.'LOCAT') then
                    do while(words(1)/='ENDLO')
                       call ecoute('ada_reastr')
                       if (words(1) .eq.'LOCAL') vloca_ada(1:ndime,iimmo) = param(1:ndime)
                       if (words(1) .eq.'GLOBA') vorig_ada(1:ndime,iimmo) = param(1:ndime)
                       if (words(1) .eq.'XYZDE') vroti_ada(1:ndime,iimmo) = param(1:ndime)
                       if (words(1) .eq.'SCALE') vscal_ada(1:ndime,iimmo) = param(1:ndime)
                    end do
                 end if
              end do
           else if (words(1) .eq. 'REMES') then
              do while(words(1)/='ENDRE')
                 call ecoute('ada_reastr')
                 if (words(1) .eq.'EXTRA') then
                    kfl_remst_ada = 1
                 else if (words(1) .eq.'SPRIN') then
                    kfl_remst_ada = 2
                 else if (words(1) .eq.'REDGR') then                    
                    kfl_remst_ada = 11
                    do while(words(1)/='ENDRE')                       
                       call ecoute('ada_reastr')
                       if (words(1) .eq.'ERROR') then
                          if (words(2) .eq.'NONDE') then
                             !
                             ! non-defined error estimator
                             !
                             kfl_redgr_ada = 1
                          else if (words(2) .eq.'ELEME') then
                             !
                             ! by-element error estimator
                             !
                             kfl_redgr_ada = 11                             
                          end if
                       end if
                    end do
                 end if
              end do
           else if (words(1) .eq. 'PREPR') then
              do while(words(1)/='ENDPR')
                 call ecoute('ada_reastr')
                 if (words(1) .eq.'ALONE') kfl_prpro_ada = 1
                 if (words(1) .eq.'MESHO') then
                    kfl_oudom_ada = 1
                    if (exists('FILEN')) fil_oudom_ada = adjustl(trim(wname))
                 else if (words(1) .eq.'GIDOU') then
                    kfl_ougid_ada = 1
                    if (exists('FILEN')) fil_ougid_ada = adjustl(trim(wname))
                 end if
              end do
           else if (words(1) .eq.'CHECK') then
              kfl_chpoi_ada = 1
           else if (words(1) .eq.'INTER') then
              kfl_inter_ada= 1                             ! nested: simple linear interp.
              if (words(2) .eq.'ELSES') kfl_inter_ada= 2   ! use elsest
           end if
           
        end do
        
     end if
     
  else if (itask == 2 .and. nimmo_ada > 0) then
     !
     ! Reading the meshes for the immersed objects (when present)
     !
     !
     ! Read and write files
     !
     lispa = 0
     lisda = lun_pdata_ada                                ! Reading file
     lisre = lun_outpu_ada                                ! Writing file
     
     if(kfl_paral<=0) then
        !
        ! Reach the section
        !
        rewind(lisda)
        call ecoute('ada_reastr')
        do while(words(1)/='ADAPT')
           call ecoute('ada_reastr')
        end do

        do while(words(1)/='ENDAD')     
           call ecoute('ada_reastr')                   
           if(words(1) .eq.'IMMER') then
              iread= 0
              iimmo= 0
              do while(words(1)/='ENDIM')
                 call ecoute('ada_reastr')
                 if (words(1) .eq.'GEOME') then
                    iimmo = iimmo + 1
                    do while(words(1)/='ENDGE')
                       call ecoute('ada_reastr')
                       if (words(1) .eq.'ELEME') then
                          do ielem= 1,nelei_ada(iimmo)
                             call ecoute('ada_reastr')
                             lnodi_ada(1:nnodi_ada,ielem,iimmo)= int(param(2:nnpar-1))
                          end do
                       else if (words(1) .eq.'COORD') then
                          ! read 
                          do ipoin=1,npoii_ada(iimmo)
                             call ecoute('ada_reastr')
                             coori_ada(1:ndime,ipoin,iimmo) = vloca_ada(1:ndime,iimmo)
                             coori_ada(1:ndime,ipoin,iimmo) = param(2:nnpar-1)
                          end do
                          ! transform
                          do ipoin= 1, npoii_ada(iimmo)
                             cauxi(1:ndime,1) = coori_ada(1:ndime,ipoin,iimmo)
                             ! scale and set the local origin in 0,0,0
                             do idime=1,ndime
                                cauxi(idime,1) = &
                                     ( cauxi(idime,1) - vloca_ada(idime,iimmo)) &
                                     * vscal_ada(idime,iimmo )
                             end do
                             ! rotate around X
                             cauxi(1,2) = cauxi(1,1)
                             cauxi(2,2) = cauxi(2,1)*cos(vroti_ada(1,iimmo)) &
                                  -       cauxi(3,1)*sin(vroti_ada(1,iimmo)) 
                             cauxi(3,2) = cauxi(2,1)*sin(vroti_ada(1,iimmo)) &
                                  +       cauxi(3,1)*cos(vroti_ada(1,iimmo)) 
                             ! rotate around Y
                             cauxi(1,2) = cauxi(1,1)*cos(vroti_ada(2,iimmo)) &
                                  +       cauxi(3,1)*sin(vroti_ada(2,iimmo)) 
                             cauxi(2,2) = cauxi(2,1)
                             cauxi(3,2) = - cauxi(1,1)*sin(vroti_ada(2,iimmo)) &
                                  +         cauxi(3,1)*cos(vroti_ada(2,iimmo)) 
                             ! rotate around Z
                             cauxi(1,2) = cauxi(1,1)*cos(vroti_ada(3,iimmo)) &
                                  -       cauxi(2,1)*sin(vroti_ada(3,iimmo)) 
                             cauxi(2,2) = cauxi(1,1)*sin(vroti_ada(3,iimmo)) &
                                  +       cauxi(2,1)*cos(vroti_ada(3,iimmo)) 
                             cauxi(3,2) = cauxi(3,1)
                             ! refer to global origin and update
                             do idime=1,ndime
                                coori_ada(idime,ipoin,iimmo) = &
                                     cauxi(idime,2) + vorig_ada(idime,iimmo)
                             end do
                          end do                          
                       end if
                    end do
                 end if
              end do
           end if
        end do
     end if

  end if

end subroutine ada_reastr
