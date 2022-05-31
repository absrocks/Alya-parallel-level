subroutine ibm_dsepar(iimbo,cooj1,cooj2,chkdi,dista,proji,projj,cooi1,cooi2)

  !-----------------------------------------------------------------------
  ! NAME
  !    nepoib
  ! DESCRIPTION
  !    Shortest signed distance from a given segment to the particle
  !    INPUT
  !       iimbo:       id of the particle
  !       nboib:       number of faces of rthe particle
  !       coor1,coor2: segment coordinates
  !       chkdi:       distance use to check the boundary box of the faces in the particle 
  !                    for find the shortest distance (use chkdi = 1.0e9_rp in general)
  !    OUTPUT
  !       dista: is less than 0 is point is inside the particle
  !              is more than 0 is point is outside the particle
  !       proji: point in iimbo nearest to the segment
  !       projj: point in the segment nearest to iimbo
  ! USED BY
  !    inouib
  !----------------------------------------------------------------------- 
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  imbou,kfl_paral,cutim,zeror
  use def_domain, only     :  ndime,mnoib,nnode
  implicit none
  integer(ip), intent(in)  :: iimbo
  real(rp),    intent(in)  :: cooj1(ndime),cooj2(ndime),chkdi
  real(rp),    intent(out) :: dista,proji(ndime),projj(ndime),cooi1(ndime),cooi2(ndime)

  integer(ip), pointer     :: blink(:),struc(:)
  real(rp),    pointer     :: sabox(:,:,:),ldist(:)

  integer(ip)              :: idime,indst,inoib,jnoib,curr,bound
  integer(ip)              :: pblty,pnodb,ipoib,iboib,jboib,kboib,ilist,jlist
  real(rp)                 :: temdi,temp,dismi1,dismi2,pladi,facto
  real(rp)                 :: tenor(ndime),norma(ndime),vect1(ndime),tepri(ndime),teprj(ndime)
  real(rp)                 :: bobox(2,ndime),bocod(ndime,mnoib),dummi
  real(rp)                 :: dsegi,dsegj,toler

  bound = 0_ip
  toler = 1.0e-2_rp

  sabox => imbou(iimbo) % sabox
  blink => imbou(iimbo) % blink
  struc => imbou(iimbo) % struc
  ldist => imbou(iimbo) % ldist

  bobox(1,1)     = min(cooj1(1)     , cooj2(1))
  bobox(2,1)     = max(cooj1(1)     , cooj2(1))
  bobox(1,2)     = min(cooj1(2)     , cooj2(2))
  bobox(2,2)     = max(cooj1(2)     , cooj2(2))  
  bobox(1,ndime) = min(cooj1(ndime) , cooj2(ndime))
  bobox(2,ndime) = max(cooj1(ndime) , cooj2(ndime))

  do idime = 1,ndime
     bobox(1,idime) = bobox(1,idime) - abs((bobox(2,idime) - bobox(1,idime))*1.0e-2)
     bobox(2,idime) = bobox(2,idime) + abs((bobox(2,idime) - bobox(1,idime))*1.0e-2)
  end do
  

  dista = 0.0_rp
  do idime = 1,ndime
     temp  = max(bobox(2,idime) - sabox(1,idime,1),sabox(2,idime,1) - bobox(1,idime))
     dista = dista + temp*temp
  end do
  
  if ( chkdi*chkdi < dista ) then
     dista = chkdi*chkdi
  end if

  indst        = 1_ip
  struc(indst) = 1_ip
  
  dismi1 = 0.0_rp
  do idime = 1,ndime
     temp  = max(0.0_rp,sabox(1,idime,1) - bobox(2,idime))  
     temp  = max(temp,  bobox(1,idime)   - sabox(2,idime,1))  
     dismi1 = dismi1 + temp * temp
  end do
  ldist(1) = dismi1
  !
  ! Assemble a list of candidate patches by traversing skd-tree
    !
  do while( indst > 0_ip )
     curr  = struc(indst)      
     indst = indst - 1_ip
     !
     ! Minimum distance
     !       
     if ( ldist(curr) < dista ) then
        !
        ! If currnode is a leaf in the tree structure
        !      
        if( blink(curr) < 0_ip ) then
           !
           ! Element properties and dimensions
           !
           iboib = -blink(curr)
           pblty = imbou(iimbo) % ltyib(iboib)
           pnodb = nnode(pblty)        
           do inoib = 1,pnodb
              ipoib = imbou(iimbo)%lnoib(inoib,iboib)
              do idime = 1,ndime
                 bocod(idime,inoib) = imbou(iimbo)%cooib(idime,ipoib)
              end do
           end do
           !
           ! Minimun distance between the segments of the differents particles
           !
           do inoib = 1,pnodb
              jnoib = inoib + 1
              if (jnoib > mnoib) jnoib = 1
              call ibm_dseseg(bocod(1,inoib),bocod(1,jnoib),cooj1,cooj2,temdi,tepri,teprj,dsegi,dsegj)
              if( temdi < dista ) then
                 if ( dsegi < 1.0_rp - toler .and. dsegi > toler .and. dsegj < 1.0_rp - toler .and. dsegj > toler ) then
                    bound = iboib
                    !
                    ! norma: Exterior normal           
                    !        
                    norma(1)     = 0.0_rp
                    norma(2)     = 0.0_rp
                    norma(ndime) = 0.0_rp
                    facto = -1.0_rp
                    
                    do ilist    = 1,imbou(iimbo)%lnele(  imbou(iimbo)%lnoib(inoib,iboib)  )%nelem
                       do jlist = 1,imbou(iimbo)%lnele(  imbou(iimbo)%lnoib(jnoib,iboib)  )%nelem
                          jboib =   imbou(iimbo)%lnele(  imbou(iimbo)%lnoib(inoib,iboib)  )%eleme(ilist)
                          kboib =   imbou(iimbo)%lnele(  imbou(iimbo)%lnoib(jnoib,iboib)  )%eleme(jlist)
                          if (jboib == kboib) then
                             call extbou(1_ip,nnode(imbou(iimbo)%ltyib(jboib)),imbou(iimbo) %lnoib(1,jboib),imbou(iimbo) %cooib,tenor)
                             do idime = 1,ndime  
                                norma(idime) = norma(idime) + tenor(idime)
                             end  do
                             call vecuni(ndime,norma,dummi)
                             
                             pladi = 0.0_rp
                             do idime = 1,ndime
                                pladi = pladi + tenor(idime) * ( teprj(idime) - bocod(idime,inoib) )
                             end do
                             if ( pladi >= 0.0_rp)  facto = 1.0_rp                 
                             
                          end if
                       end  do
                    end  do
                    dista = temdi
                    
                    proji(1)     = tepri(1)
                    proji(2)     = tepri(2)
                    proji(ndime) = tepri(ndime)
                    projj(1)     = teprj(1)
                    projj(2)     = teprj(2)
                    projj(ndime) = teprj(ndime)
                    
                    cooi1(1)     = bocod(1,inoib)
                    cooi1(2)     = bocod(2,inoib)
                    cooi1(ndime) = bocod(ndime,inoib)
                    cooi2(1)     = bocod(1,jnoib)
                    cooi2(2)     = bocod(2,jnoib)
                    cooi2(ndime) = bocod(ndime,jnoib)
                 
                 end if
              end if
           end do
        else                         
           dismi1 = 0.0_rp
           dismi2 = 0.0_rp
           do idime = 1,ndime
              temp = 0.0_rp
              temp  = max(0.0_rp,sabox(1,idime,blink(curr)) - bobox(2,idime))  
              temp  = max(temp,  bobox(1,idime)             - sabox(2,idime,blink(curr)))  
              dismi1 = dismi1 + temp * temp
              
              temp  = max(0.0_rp,sabox(1,idime,blink(curr)+1_ip) - bobox(2,idime))  
              temp  = max(temp,  bobox(1,idime)                  - sabox(2,idime,blink(curr)+1_ip))  
              dismi2 = dismi2 + temp * temp
           end do
           ldist(blink(curr))    = dismi1             
           ldist(blink(curr)+1)  = dismi2
           
           indst = indst + 2_ip
           if (dismi1 > dismi2) then               
              struc(indst-1) = blink(curr)
              struc(indst)   = blink(curr) + 1_ip
           else
              struc(indst-1) = blink(curr) + 1_ip
              struc(indst)   = blink(curr)                              
           end if
        end if
     end if
  end do

  if (bound /= 0_ip) then
     dista = sqrt(dista)*facto
  else
     dista = 1.0e10_rp
  end if

     
end subroutine ibm_dsepar
