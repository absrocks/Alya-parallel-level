subroutine elsest_statis(itask,imesh,ipara,ithre)
  !------------------------------------------------------------------------
  !****f* elsest/elsest_statis
  ! NAME 
  !    elsest_statis
  ! DESCRIPTION
  !    Output statistics
  ! USES
  ! USED BY
  !    elsest
  !***
  !------------------------------------------------------------------------
  use def_elsest
  implicit none
  integer(ip), intent(in) :: itask,imesh,ipara(*),ithre
  real(rp)                :: cputo,cputo_master
  integer(ip)             :: tomem,ii,dummi,jj
  real(rp)                :: rbyte,rbyt2,rbyt3,dummr
  character(2)            :: lbyte,lbyt2,lbyt3

  if(iunit(1)/=0) then

     select case(itask)

     case(1)
        !
        ! CPU time
        !
        if(ipara(8)==0) then
           !*OMP CRITICAL(write)
           write(iunit(1),1) imesh,'(BIN STRATEGY)'
           cputo=cputi(1,ithre)+cputi(2,ithre)+cputi(3,ithre)+cputi(4,ithre)
           !*OMP END CRITICAL(write)
           if(cputo==0.0_rp) cputo=1.0_rp
           !*OMP CRITICAL(write)
           write(iunit(1),100) cputo,&
                cputi(1,ithre),100.0_rp*cputi(1,ithre)/cputo,&
                cputi(2,ithre),100.0_rp*cputi(2,ithre)/cputo,&
                cputi(3,ithre),100.0_rp*cputi(3,ithre)/cputo,&
                cputi(4,ithre),100.0_rp*cputi(4,ithre)/cputo
           !*OMP END CRITICAL(write)
        else
           !*OMP CRITICAL(write)
           write(iunit(1),1) imesh,'(QUAD/OCT TREE STRATEGY)'
           !*OMP END CRITICAL(write)
           cputo=cputi(1,ithre)+cputi(2,ithre)+cputi(3,ithre)+cputi(4,ithre)
           if(cputo==0.0_rp) cputo=1.0_rp
           !*OMP CRITICAL(write)
           write(iunit(1),101) cputo,&
                cputi(1,ithre),100.0_rp*cputi(1,ithre)/cputo,&
                cputi(2,ithre),100.0_rp*cputi(2,ithre)/cputo,&
                cputi(3,ithre),100.0_rp*cputi(3,ithre)/cputo,&
                cputi(4,ithre),100.0_rp*cputi(4,ithre)/cputo
           !*OMP END CRITICAL(write)
        end if
        !
        ! Memory
        !
        tomem=0
        if(nthre>1) then
           do ii=2,5
              do jj=1,nthre
                 tomem=tomem+memor(ii,jj)
              end do
           end do
           do jj=2,nthre
              memor(2,1)=memor(2,1)+memor(2,jj)
              memor(3,1)=memor(3,1)+memor(3,jj)
              memor(5,1)=memor(5,1)+memor(5,jj)
           end do
        else
           do ii=2,5
              tomem=tomem+memor(ii,ithre)
           end do
        end if
        if(tomem==0) tomem=1

        call elsest_inbyte(tomem,rbyte,lbyte)
        call elsest_inbyte(memor(5,ithre),rbyt2,lbyt2)
        call elsest_inbyte(memax,rbyt3,lbyt3)

        if(ipara(8)==0) then
           !*OMP CRITICAL(write)
           write(iunit(1),200) real(tomem)/rbyte,lbyte,&
                real(memor(2,ithre))/rbyte,lbyte,100.0_rp*real(memor(2,ithre))/real(tomem),&
                real(memor(3,ithre))/rbyte,lbyte,100.0_rp*real(memor(3,ithre))/real(tomem),&
                real(memax)/rbyt3,lbyt3
           !*OMP END CRITICAL(write)
        else if(ipara(8)==1) then
           !*OMP CRITICAL(write)
           write(iunit(1),201) real(tomem)/rbyte,lbyte,&
                real(memor(2,ithre))/rbyte,lbyte,100.0_rp*real(memor(2,ithre))/real(tomem),&
                real(memor(3,ithre))/rbyte,lbyte,100.0_rp*real(memor(3,ithre))/real(tomem),&
                real(memax)/rbyt3,lbyt3
           !*OMP END CRITICAL(write)
        end if
        !
        ! Structure statistics
        !
        if(ipara(8)==0) then
           if(kstat(6,ithre)/=0) then
              dummr=real(kstat(6,ithre))
           else
              dummr=1.0_rp
           end if
           !*OMP CRITICAL(write)
           write(iunit(1),300)&
                kstat(6,ithre),kstat(5,ithre),100.0_rp*real(kstat(5,ithre))/dummr,kstat(1,ithre),kstat(2,ithre),&
                int(real(kstat(8,ithre))/dummr),kstat(3,ithre),kstat(4,ithre),int(real(kstat(7,ithre))/dummr)
           !*OMP END CRITICAL(write)
        else if(ipara(8)==1) then
           if(kstat(6,ithre)/=0) then
              dummr=real(kstat(6,ithre))
           else
              dummr=1.0_rp
           end if
           !*OMP CRITICAL(write)
           write(iunit(1),301)&
                kstat(6,ithre),kstat(1,ithre),kstat(2,ithre),int(real(kstat(8,ithre))/dummr),&
                kstat(3,ithre),kstat(4,ithre),int(real(kstat(7,ithre))/dummr)
           !*OMP END CRITICAL(write)
        end if

     case(2)
        !
        ! Process
        !
        if(ipara(8)==0) then
           !*OMP CRITICAL(write)
           write(iunit(1),2) imesh,'(BIN STRATEGY)',ithre
           !*OMP END CRITICAL(write)
        else if(ipara(8)==1) then
           !*OMP CRITICAL(write)
           write(iunit(1),2) imesh,'(QUAD/OCT TREE STRATEGY)'
           !*OMP END CRITICAL(write)
        end if

        !
        ! USING OMP
        !
        if(ksear(ithre)/=0) then
           if (nthre>1) then
              dummr=0.
              do ii=1,nthre
                 dummr=dummr+real(ksear(ii))
              end do
           else
              dummr=real(ksear(ithre))
           end if
        else
           dummr=1.0_rp
        end if
        if(nthre>1) then
           cputo=0.
           cputo_master=cputi(6,ithre)+cputi(7,ithre)
           do ii=1,nthre
              cputo=cputo+cputi(6,ii)+cputi(7,ii)
           end do
           do ii=2,nthre
              ksear(1)   = ksear(1)+ksear(ii)
              kfirs(1)   = kfirs(1)+kfirs(ii)
              kseco(1)   = kseco(1)+kseco(ii)
              cputi(1,1) = cputi(1,1)+cputi(1,ii)
              cputi(2,1) = cputi(2,1)+cputi(2,ii)
              cputi(3,1) = cputi(3,1)+cputi(3,ii)
              cputi(4,1) = cputi(4,1)+cputi(4,ii)
              cputi(5,1) = cputi(5,1)+cputi(5,ii)
              cputi(6,1) = cputi(6,1)+cputi(6,ii)
              cputi(7,1) = cputi(7,1)+cputi(7,ii)
           end do
        else
           cputo=cputi(6,ithre)+cputi(7,ithre)
           cputo_master=cputi(6,ithre)+cputi(7,ithre)
        end if
        if(cputo==0.0_rp) cputo=1.0_rp
        dummi=ksear(ithre)-kfirs(ithre)-kseco(ithre)
        !*OMP CRITICAL(write)
        write(iunit(1),402)&
             ksear(ithre),&
             nthre,&
             kfirs(ithre),&
             100.0_rp*real(kfirs(ithre))/real(ksear(ithre)),&
             kseco(ithre),100.0_rp*real(kseco(ithre))/real(ksear(ithre)),&
             dummi,100.0_rp*real(dummi)/real(ksear(ithre)),&
             cputo,&
             cputo_master,&
             cputi(6,ithre),100.0_rp*cputi(6,ithre)/cputo,&
             cputi(7,ithre),100.0_rp*cputi(7,ithre)/cputo
        !*OMP END CRITICAL(write)
     end select

  end if

1 format(/,&
       & 5x,'   ----------------','---',/&
       & 5x,'|- PREPROCESS MESH ',i3,' ',a,/,&
       & 5x,'   ----------------','---')
2 format(/,&
       & 5x,'   ----------------','---',/&
       & 5x,'|- PROCESS MESH    ',i3,' ',a,/,&
       & 5x,'   ----------------','---')
100 format(/,&
       & 5x,'   |- SUMMARY OF COMPUTING TIMES  ',//,&
       & 5x,'      TOTAL CPUT TIME:            ',f11.2,3x,/,&
       & 5x,'      INTIALIZATION:              ',f11.2,3x,' (',f6.2,' % )',/,&
       & 5x,'      FILL BIN WITH NODES:        ',f11.2,3x,' (',f6.2,' % )',/,&
       & 5x,'      FILL BIN WITH ELEMENTS:     ',f11.2,3x,' (',f6.2,' % )',/,&
       & 5x,'      SECOND TRY STRATEGY:        ',f11.2,3x,' (',f6.2,' % )')
101 format(/,&
       & 5x,'   |- SUMMARY OF COMPUTING TIMES  ',//,&
       & 5x,'      TOTAL CPUT TIME:            ',f11.2,3x,/,&
       & 5x,'      INITIALIZATION:             ',f11.2,3x,' (',f6.2,' % )',/,&
       & 5x,'      QUAD/OCT STRUCTURE:         ',f11.2,3x,' (',f6.2,' % )',/,&
       & 5x,'      FILL BIN WITH ELEMENTS:     ',f11.2,3x,' (',f6.2,' % )',/,&
       & 5x,'      SECOND TRY STRATEGY:        ',f11.2,3x,' (',f6.2,' % )')
200 format(/,&
       & 5x,'   |- SUMMARY OF MEMORY           ',//,&
       & 5x,'      TOTAL MEMORY:               ',f11.2,1x,a2,/,&
       & 5x,'      BIN STRUCTURE:              ',f11.2,1x,a2,' (',f6.2,' % )',/,&
       & 5x,'      SECOND TRY STRATEGY:        ',f11.2,1x,a2,' (',f6.2,' % )',/,&
       & 5x,/,&
       & 5x,'      MAXIMUM MEMORY REQUIRED:    ',f11.2,1x,a2)
201 format(/,&
       & 5x,'   |- SUMMARY OF MEMORY           ',//,&
       & 5x,'      TOTAL MEMORY:               ',f11.2,1x,a2,/,&
       & 5x,'      QUAD/OCT STRUCTURE:         ',f11.2,1x,a2,' (',f6.2,' % )',/,&
       & 5x,'      SECOND TRY STRATEGY:        ',f11.2,1x,a2,' (',f6.2,' % )',/,&
       & 5x,/,&
       & 5x,'      MAXIMUM MEMORY REQUIRED:    ',f11.2,1x,a2)
300 format(/,&
       & 5x,'   |- SUMMARY OF BIN STRUCTURE    ',//,&
       & 5x,'      # OF BINS:                  ',i11,/,&
       & 5x,'      # OF BINS WITH ELEMENTS:    ',i11,3x,' (',f6.2,' % )',/,&
       & 5x,'      MIN # NODES IN A BIN:       ',i11,/,&
       & 5x,'      MAX # NODES IN A BIN:       ',i11,/,&
       & 5x,'      AVERAGE # NODES IN BIN:     ',i11,/,&
       & 5x,'      MIN # ELEMENTS IN A BIN:    ',i11,/,&
       & 5x,'      MAX # ELEMENTS IN A BIN:    ',i11,/,&
       & 5x,'      AVERAGE # ELEMENTS IN BIN:  ',i11)
301 format(/,&
       & 5x,'   |- SUMMARY OF QUAD/OCT STRUCTURE',//,&
       & 5x,'      # OF BINS WITH ELEMENTS:    ',i11,/,&
       & 5x,'      MIN # NODES IN A BIN:       ',i11,/,&
       & 5x,'      MAX # NODES IN A BIN:       ',i11,/,&
       & 5x,'      AVERAGE # NODES IN BIN:     ',i11,/,&
       & 5x,'      MIN # ELEMENTS IN A BIN:    ',i11,/,&
       & 5x,'      MAX # ELEMENTS IN A BIN:    ',i11,/,&
       & 5x,'      AVERAGE # ELEMENTS IN BIN:  ',i11)
400 format(/,&
       & 5x,'   |- ELEMENT SEARCH:             ',//,&
       & 5x,'      # OF SEARCHES:              ',i11,/,&
       & 5x,'      # ELEMS FOUND 1ST STRATEGY: ',i11,  3x,   ' (',f6.2,' % )',/,&
       & 5x,'      # ELEMS FOUND 2ND STRATEGY: ',i11,  3x,   ' (',f6.2,' % )',/,&
       & 5x,'      # ELEMS NOT FOUND:          ',i11  ,3x,   ' (',f6.2,' % )',/,&
       & 5x,'      TOTAL CPU TIME :            ',f11.2,3x,                      /,&
       & 5x,'      CPU TIME 1ST STRATEGY:      ',f11.2,3x,   ' (',f6.2,' % )',/,&
       & 5x,'      CPU TIME 2ND STRATEGY:      ',f11.2,3x,   ' (',f6.2,' % )')
402 format(/,&
       & 5x,'   |- ELEMENT SEARCH:             ',//,&
       & 5x,'      # OF SEARCHES:              ',i11,/,&
       & 5x,'      # OF THREADS:               ',i11,/,&
       & 5x,'      # ELEMS FOUND 1ST STRATEGY: ',i11,  3x,   ' (',f6.2,' % )',/,&
       & 5x,'      # ELEMS FOUND 2ND STRATEGY: ',i11,  3x,   ' (',f6.2,' % )',/,&
       & 5x,'      # ELEMS NOT FOUND:          ',i11  ,3x,   ' (',f6.2,' % )',/,&
       & 5x,'      TOTAL CPU TIME :            ',f11.2,3x,                      /,&
       & 5x,'      CPU TIME MASTER THREAD:     ',f11.2,3x,                      /,&
       & 5x,'      CPU TIME 1ST STRATEGY:      ',f11.2,3x,   ' (',f6.2,' % )',/,&
       & 5x,'      CPU TIME 2ND STRATEGY:      ',f11.2,3x,   ' (',f6.2,' % )')

end subroutine elsest_statis
