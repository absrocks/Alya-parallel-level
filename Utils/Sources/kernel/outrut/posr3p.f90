subroutine posr3p(bridge,wopos,itste,ttime)

  !-----------------------------------------------------------------------
  !
  ! Write a matrix at Gauss points in postprocess file
  !
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_elmtyp
  implicit none
  character(*), intent(in)  :: wopos
  type(r3p),    intent(in)  :: bridge(:) 
  integer(ip) , intent(in)  :: itste
  real(rp)    , intent(in)  :: ttime
  integer(ip)               :: ielem,igaus
  character(8)              :: state
  character(4)              :: CNUL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
  integer(ip),  save        :: ipass=0
  integer(ip)               :: iesta,iesto,ielty,icomp,ncomp
  character(13)             :: elemt

  state='ANALYSIS'

  if(kfl_outfo==1) then
     if(ndime==2) then
        iesta=10
        iesto=29
     else if(ndime==2) then
        iesta=30
        iesto=50
     end if
     if(ipass==0) then
        ipass=1
        do ielty=iesta,iesto
           if(lexis(ielty)==1) then
              if(ndime==2) then
                 if(nnode(ielty)==3.or.nnode(ielty)==6.or.nnode(ielty)==7) then
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
              write(lun_postp,1) 'GaussPoints '//'GP_'//trim(cenam(ielty))&
                   &             //' Elemtype '//trim(elemt)
              write(lun_postp,5) 'Number of Gauss Points: ',ngaus(ielty)
              write(lun_postp,1) 'Natural Coordinates: Internal'
              write(lun_postp,1) 'End GaussPoints'
           end if
        end do
     end if
     do ielty=iesta,iesto
        if(lexis(ielty)==1) then 
           ncomp=size(bridge(1)%a,KIND=ip)
           if(ncomp==ndime) then
              write(lun_postp,2) wopos,state,ttime,'Vector','GP_'//trim(cenam(ielty))
              write(lun_postp,3) wopos//'_X,'//wopos//'_Y,'//wopos//'_Z'
              write(lun_postp,1) 'Values'
              do ielem=1,nelem
                 if(ltype(ielem)==ielty) then
                    write(lun_postp,7,advance='no') ielem
                    do igaus=1,ngaus(ielty)
                       write(lun_postp,6) (bridge(ielem)%a(icomp,igaus,1),icomp=1,ncomp)
                    end do
                 end if
              end do
           else
              write(lun_postp,2) wopos,state,ttime,'Matrix','GP_'//trim(cenam(ielty))
              if(ndime==2) then
                 write(lun_postp,3) wopos//'_XX,'//wopos//'_YY,'//wopos//'_ZZ,'//wopos//'_XY'
              else
                 write(lun_postp,3) wopos//'_XX,'//wopos//'_YY,'//wopos//'_ZZ,'//wopos//'_XY,'&
                      //wopos//'_YZ,'//wopos//'_XZ'
              end if
              write(lun_postp,1) 'Values'
              if(ndime==2) then
                 do ielem=1,nelem
                    if(ltype(ielem)==ielty) then
                       write(lun_postp,7,advance='no') ielem                    
                       do igaus=1,ngaus(ielty)
                          write(lun_postp,6) &
                               bridge(ielem)%a(1,igaus,1),&
                               bridge(ielem)%a(2,igaus,1),&
                               bridge(ielem)%a(4,igaus,1),&
                               bridge(ielem)%a(3,igaus,1)
                       end do
                    end if
                 end do
              else
                 do ielem=1,nelem
                    if(ltype(ielem)==ielty) then
                       write(lun_postp,7,advance='no') ielem                    
                       do igaus=1,ngaus(ielty)
                          write(lun_postp,6) &
                               bridge(ielem)%a(1,igaus,1),&
                               bridge(ielem)%a(2,igaus,1),&
                               bridge(ielem)%a(4,igaus,1),&
                               bridge(ielem)%a(3,igaus,1),&
                               bridge(ielem)%a(6,igaus,1),&
                               bridge(ielem)%a(5,igaus,1)
                       end do
                    end if
                 end do
              end if
           end if
           write(lun_postp,1) 'End Values'
        end if
     end do
     flush(lun_postp)
  else if(kfl_outfo==3.or.kfl_outfo==4) then
     call GiD_BeginScalarResult(wopos,state,ttime,1,CNUL,CNUL,CNUL)
     do ielem=1,nelem
        !call GiD_WriteScalar(ielem,bridge(ielem)%a(:,:,1))
     end do
     CALL GiD_EndResult
  end if
  !
  ! GiD formats
  !
1 format(a)
2 format('Result ',a,' ',a,' ',e12.6,' ',a,' OnGaussPoints ',a)
3 format('ComponentNames ',a)
4 format(i7, 3(1x,e16.8E3))
5 format(a,1x,i2)
6 format(6(1x,e16.8E3))
7 format(i7)
  !
  ! Femview formats
  !
10 format(1x,i4,a1,a6,e12.5,32x,i2,i5)
20 format(1x,i2,2x,a5,3x,3i5)
30 format(1x,i2,2x,a5,3x,2i5)
40 format(1x,i2,i5,e12.5)
50 format(1x,i2)


end subroutine posr3p
