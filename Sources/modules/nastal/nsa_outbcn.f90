subroutine nsa_outbcn
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      mod_memchk
  
  use      def_nastal
  use mod_messages, only : livinf
  
  implicit none
  character(5)  :: chaux(mnodb)
  character(300)           :: messa, numstring
  
  integer(ip)   :: idime,ielem,idofn,ipoin,iboun,inode,inodb,ncouf(nbhie_nsa),ibhie,kcoun,nnodb

  !
  ! Some global checks on the physical setup
  !

  call nsa_evaluate_uinlet(0_ip,0_ip)


  !
  ! Master to output a summary on screen
  !
  if (INOTSLAVE) then
     
     messa = &
          '        MODULE SETUP SUMMARY:'
     call livinf(0_ip,messa,one)
     messa = '        REFERENCE VALUES:'
     call livinf(0_ip,messa,one)
     write(numstring,'(d17.10)') speed_nsa
     messa = '          SPEED= '//numstring
     call livinf(0_ip,messa,one)
     write(numstring,'(d17.10)') press_nsa
     messa = '          PRESSURE= '//numstring
     call livinf(0_ip,messa,one)
     write(numstring,'(d17.10)') densi_nsa
     messa = '          DENSITY= '//numstring
     call livinf(0_ip,messa,one)
     write(numstring,'(d17.10)') tempe_nsa
     messa = '          TEMPERATURE= '//numstring
     call livinf(0_ip,messa,one)
     write(numstring,'(d17.10)') press_dynamic_nsa
     messa = '          DYNAMIC PRESSURE= '//numstring
     call livinf(0_ip,messa,one)
     write(numstring,'(d17.10)') tstag_nsa
     messa = '          STAGNATION TEMPERATURE= '//numstring
     call livinf(0_ip,messa,one)
     write(numstring,'(d17.10)') rmach_nsa
     messa = '          MACH NUMBER= '//numstring
     call livinf(0_ip,messa,one)
     write(numstring,'(d17.10)') xsoun_nsa
     messa = '          SOUND SPEED= '//numstring
     call livinf(0_ip,messa,one)
     write(numstring,'(d17.10)') cpcoe_nsa
     messa = '          C_p= '//numstring
     call livinf(0_ip,messa,one)
     write(numstring,'(d17.10)') cvcoe_nsa
     messa = '          C_v= '//numstring
     call livinf(0_ip,messa,one)
     write(numstring,'(d17.10)') rgasc_nsa
     messa = '          UNIVERSAL GAS CONSTANT (R)= '//numstring
     call livinf(0_ip,messa,one)
     if (kfl_visco_nsa == 0) then
        messa = '          VISCOSITY TERMS SET TO OFF: EULER EQUATIONS '
        call livinf(0_ip,messa,one)
     else
        messa = '          VISCOSITY TERMS SET TO ON: NAVIER-STOKES EQUATIONS '
        call livinf(0_ip,messa,one)
     end if
     write(numstring,'(d17.10)') visco_nsa
     messa = '            DYNAMIC VISCOSITY (mu)= '//numstring
     call livinf(0_ip,messa,one)
     write(numstring,'(d17.10)') thdif_nsa
     messa = '            THERMAL DIFFUSION (mu*C_p / Pr)= '//numstring
     call livinf(0_ip,messa,one)
     messa = '        BOUNDARY CONDITIONS:'
     call livinf(0_ip,messa,one)
     if (nuinlet_nsa > 0) then
        write(numstring,'(i10)') nuinlet_nsa
        messa = '          TOTAL PRESSURE CONDITION DETECTED ON # NODES: '//numstring
        call livinf(0_ip,messa,one)
        if (kfl_mod_elmop_nsa > 0) then
           call runend('MOD_NSA_ELMOPERATIONS:TOTAL PRESSURE NOT YET IMPLEMENTED IN ELMOPERATIONS') 
        end if
     end if
  end if


  if (kfl_dumbo_nsa==0) return
  if (kfl_dumbo_nsa==1) then
     if (nboun==0) call runend('OUTBCN: NO NBOUN DEFINED')

!     write(6,*) 'uieee ',chhie_nsa(1:nbhie_nsa)

     write (lun_dumb1_nsa,*) '# CODES'
     do ibhie=1,nbhie_nsa
        write (lun_dumb1_nsa,*) '# ',ibhie,chhie_nsa(ibhie)(1:ndofn_nsa)
     end do
     write (lun_dumb1_nsa,*) '# END_CODES'
     write (lun_dumb1_nsa,'(a)') 'ON_BOUNDARIES, CODED'
     write (lun_dumb2_nsa,'(a)') 'BOUNDARIES'
     do iboun=1,nboun
        ncouf= 0
        chaux= ''
        nnodb=nnode(ltypb(iboun))
        do inodb=1,nnodb
           ipoin=lnodb(inodb,iboun)
           do idofn=1,ndofn_nsa
              chaux(inodb)(idofn:idofn)= intost(kfl_fixno_nsa(idofn,ipoin))
           end do
           do ibhie=1,nbhie_nsa
              if (chaux(inodb)(1:ndofn_nsa)==chhie_nsa(ibhie)(1:ndofn_nsa)) ncouf(ibhie)=ncouf(ibhie)+1
           end do
        end do
        kcoun=0
        do ibhie=1,nbhie_nsa
           if (ncouf(ibhie) > 0) then
              kcoun=ibhie
              exit
           end if
        end do
        
!             write(6,*) kcoun,ncouf(1:nbhie_nsa),chaux(1:nnodb)
        
        if (kcoun==0) then
           write (6,*) ' | '
           write (6,*) ' |------>   Missing fixity:  ',chaux(1:nnodb)
           write (6,*) ' | '
           call runend('OUTBCN: ALL FIXITIES MUST BE DECLARED IN THE HIERARCHY!')
        end if
        
        write (lun_dumb1_nsa,*) iboun,kcoun
        write (lun_dumb2_nsa,*) iboun,lnodb(1:nnodb,iboun),lelbo(iboun)

     end do
     
     write (lun_dumb1_nsa,'(a)') 'END_ON_BOUNDARIES'
     write (lun_dumb2_nsa,'(a)') 'END_BOUNDARIES'

  else if (kfl_dumbo_nsa==2) then

     call runend('OUTCON: BOUNDARY-TO-NODE NOT PROGRAMMED YET')

  end if


end subroutine nsa_outbcn
