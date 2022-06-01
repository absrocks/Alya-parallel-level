      program Changenod
c******************************************************************************
c
c**** This program changes the element type of a finite element mesh.
c**** The possibilities are the following:      
c
c**** For NDIME = 2
c      
c****     NNODE = 3     ---->    NNODE = 4, 6 & 7     (more nodes are created)
c****     NNODE = 4 P1+ ---->    NNODE = 3            (same number of nodes)
c****     NNODE = 4     ---->    NNODE = 4, 5, 9 & 16 (more nodes are created)
c****     NNODE = 6     ---->    NNODE = 7            (more nodes are created)
c****     NNODE = 6     ---->    NNODE = 3            (same number of nodes)
c****     NNODE = 9     ---->    NNODE = 3, 4 , 6, 7  (same number of nodes)
c****     NNODE = 16    ---->    NNODE = 10           (same number of nodes)
c
c**** For NDIME = 3
c
c****     NNODE = 10    ---->    NNODE = 4       (same number of nodes)     
c****     NNODE = 4     ---->    NNODE = 5 & 10  (more nodes are created)
c****     NNODE = 8     ---->    NNODE = 20      (more nodes are created)
c
c**** The data from the initial mesh is supposed to be read from the .geo
c**** file of FANTOM. The results are also presented in this format.
c
c******************************************************************************
      implicit     real*8 (a-h,o-z)
      common/contr/nin,nou,ndime,nelem,nnode,npoin
      parameter    (mxnod=20, mxelm=40000, mxpoi=40000, mxdim=3,
     .              mxnow=20, mxelw=100000, mxpow=80000)
      dimension    lnods(mxelm,mxnod), lnodw(mxelw,mxnow),
     .             coord(mxdim,mxpoi), coorw(mxdim,mxpow),
     .             tempo(mxdim,mxpow), lpoiw(      mxpow)
c
c***  Open files
c
      call opfile
c
c***  Identify control parameters
c
      call contro
c
c***  Read geometry (initial mesh)
c
      call reageo(lnods,coord)
c
c***  Undertakes the mesh change
c
      call mesdiv(mxnow,mxelw,mxpow,coord,coorw,lnods,lnodw,
     .  npoiw,nelew,nnodw)
c
c***  Checks if there are repeated nodes and output of results
c
      call mescek(mxnow,mxelw,mxpow,coord,coorw,lnods,lnodw,
     .  tempo,lpoiw,npoiw,nelew,nnodw)

      stop
      end
c*************************************************************************     
      subroutine opfile
      common/contr/nin,nou,ndime,nelem,nnode,npoin
      character*20 file_mesh, file_resu

      write(6,'(a,$)') ' >>> Original mesh file:'
      read(5,'(a)') file_mesh
      write(6,'(a,$)') ' >>> Final mesh file:'
      read(5,'(a)') file_resu
      nin=7
      nou=8
      open(nin,file=file_mesh,status='old')
      open(nou,file=file_resu,status='unknown')

      end
c*************************************************************************     
      subroutine contro
      implicit     real*8 (a-h,o-z)
      common/contr/nin,nou,ndime,nelem,nnode,npoin
      character*80 string
      character*5  wopos
      
      ndime=-1
      i_chk= 1
      nelem= 0
      npoin= 0
      icoun= 0
      nnode=-1
c      
      wopos=''
      do while(wopos.ne.'ELEME'.and.wopos.ne.' ELEM')
        read(nin,'(a5)')wopos
      end do
   20 jcoun= 1
      read(nin,'(a80)') string
   10 do while(string(jcoun:jcoun).eq.' '.and.jcoun.le.80)
        jcoun = jcoun+1
      enddo
      kcoun = 1
      do while(string(jcoun:jcoun).ne.' '.and.jcoun.le.80)
        if(string(jcoun:jcoun).eq.'/') then
          nelem = nelem+1
          i_chk = i_chk+1
          goto 20
        endif
        if(kcoun.eq.1) nnode = nnode+1
        kcoun = kcoun+1
        jcoun = jcoun+1
      enddo
      if(jcoun.lt.80) goto 10
      do while(string(1:7).ne.'END_ELE')
        read(nin,'(a80)') string
        nelem=nelem+1
      enddo
      rewind(nin)
      wopos=''
      do while(wopos.ne.'COORD'.and.wopos.ne.' COOR')
        read(nin,'(a5)')wopos
      end do     
      read(nin,'(a80)') string
      jcoun= 1
   30 do while(string(jcoun:jcoun).eq.' '.and.jcoun.le.80)
        jcoun = jcoun+1
      enddo
      kcoun = 1
      do while(string(jcoun:jcoun).ne.' '.and.jcoun.le.80)
        if(kcoun.eq.1) ndime = ndime+1
        kcoun = kcoun+1
        jcoun = jcoun+1
      enddo
      if(jcoun.lt.80) goto 30
      do while(string(1:7).ne.'END_COO')
        read(nin,'(a80)') string
        npoin=npoin+1
      enddo
      nelem=nelem/i_chk
      rewind(nin)
c
      write(6,'(a   )') '---------------------------------'
      write(6,'(a,i7)') 'NUMBER OF ELEMENTS FOUND: ',nelem
      write(6,'(a,i7)') 'NUMBER OF NODES FOUND:    ',npoin
      write(6,'(a,i7)') 'SPACE DIMENSION FOUND:    ',ndime
      write(6,'(a   )') '---------------------------------'
c      
      end
c*************************************************************************     
      subroutine reageo(lnods,coord)
      implicit     real*8 (a-h,o-z)
      common/contr/nin,nou,ndime,nelem,nnode,npoin
      dimension    coord(ndime,npoin), lnods(nelem,nnode)
      character*5  wopos

      wopos=''
      do while(wopos.ne.'ELEME'.and.wopos.ne.' ELEM')
        read(nin,'(a5)')wopos
      end do
      do ielem = 1,nelem
        read(nin,*) jelem,(lnods(jelem,inode),inode=1,nnode)
      enddo
      rewind(nin)
      wopos=''
      do while(wopos.ne.'COORD'.and.wopos.ne.' COOR')
        read(nin,'(a5)')wopos
      end do
      do ipoin = 1,npoin
        read(nin,*) jpoin,(coord(idime,jpoin),idime=1,ndime)
      enddo

      end
c*************************************************************************     
      subroutine mesdiv(mxnow,mxelw,mxpow,coord,coorw,lnods,lnodw,
     .  npoiw,nelew,nnodw)
      implicit     real*8 (a-h,o-z)
      common/contr/nin,nou,ndime,nelem,nnode,npoin
      dimension    coord(ndime,npoin), lnods(nelem,nnode),
     .  coorw(ndime,mxpow), lnodw(mxelw,mxnow)
c
c***  Initializations
c
      do idime=1,ndime
        do ipoin=1,npoin
          coorw(idime,ipoin)=coord(idime,ipoin)
        end do
      end do
      nelew=nelem
      npoiw=npoin
c
c***  Splits the elements
c
      write(6,'(a,$)') ' >>> Nodes of the final elements:'
      read(5,*) nnodw
      
      do ielem=1,nelem
        if(ndime.eq.2) then
c
c*** 2D: NNODE = 3 --> NNODE = 4 6 & 7
c          
          if(nnode.eq.3) then
            if (nnodw.ge.6) then
              do idime=1,2
                coorw(idime,npoiw+1)=0.5*(coord(idime,lnods(ielem,1))
     .            +coord(idime,lnods(ielem,2)))
                coorw(idime,npoiw+2)=0.5*(coord(idime,lnods(ielem,2))
     .            +coord(idime,lnods(ielem,3)))
                coorw(idime,npoiw+3)=0.5*(coord(idime,lnods(ielem,3))
     .            +coord(idime,lnods(ielem,1)))
              end do
              do inode=1,3
                lnodw(ielem,  inode)=lnods(ielem,inode)
                lnodw(ielem,3+inode)=npoiw+inode
              end do
              npoiw=npoiw+3
              if(nnodw.eq.7) then
                do idime=1,2
                  coorw(idime,npoiw+1)=(1.0/3.0)*(coord(idime,lnods(ielem,1))
     .              +coord(idime,lnods(ielem,2))+coord(idime,lnods(ielem,3)))
                end do
                lnodw(ielem,7)=npoiw+1
                npoiw=npoiw+1
              end if
            else if (nnodw.eq.4) then
              do idime=1,2
                coorw(idime,npoiw+1)=1./3.*(coord(idime,lnods(ielem,1))
     .            +coord(idime,lnods(ielem,2))
     .            +coord(idime,lnods(ielem,3)))
              end do
              npoiw=npoiw+1
              do inode=1,3
                lnodw(ielem,inode)=lnods(ielem,inode)
              end do
              lnodw(ielem,4)=npoiw
            end if
c
c*** 2D:  NNODE = 6 --> NNODE = 7
c
          else if (nnode.eq.6) then
            if (nnodw.eq.3) then
              lnodw(ielem,1)=lnods(ielem,1)
              lnodw(ielem,2)=lnods(ielem,4)
              lnodw(ielem,3)=lnods(ielem,6)
              nelew=nelew+1
              lnodw(nelew,1)=lnods(ielem,4)
              lnodw(nelew,2)=lnods(ielem,2)
              lnodw(nelew,3)=lnods(ielem,5)
              nelew=nelew+1
              lnodw(nelew,1)=lnods(ielem,4)
              lnodw(nelew,2)=lnods(ielem,5)
              lnodw(nelew,3)=lnods(ielem,6)
              nelew=nelew+1
              lnodw(nelew,1)=lnods(ielem,6)
              lnodw(nelew,2)=lnods(ielem,5)
              lnodw(nelew,3)=lnods(ielem,3)
            else if (nnodw.eq.7) then
              do idime=1,2
                coorw(idime,npoiw+1)=(1.0/3.0)*(coord(idime,lnods(ielem,1))
     .            +coord(idime,lnods(ielem,2))+coord(idime,lnods(ielem,3)))
              end do
              do inode=1,6
                lnodw(ielem,inode)=lnods(ielem,inode)
              end do
              lnodw(ielem,7)=npoiw+1
              npoiw=npoiw+1
            end if
c
c*** 2D:  NNODE = 4 --> NNODE = 4 (split), 5, 9 & 16
c
          else if(nnode.eq.4) then
            if (nnodw.eq.3) then
              continue
            else if (nnodw.eq.4.or.nnodw.eq.9) then
              do idime=1,2
                coorw(idime,npoiw+1)=0.5*(coord(idime,lnods(ielem,1))
     .            +coord(idime,lnods(ielem,2)))
                coorw(idime,npoiw+2)=0.5*(coord(idime,lnods(ielem,2))
     .            +coord(idime,lnods(ielem,3)))
                coorw(idime,npoiw+3)=0.5*(coord(idime,lnods(ielem,3))
     .            +coord(idime,lnods(ielem,4)))
                coorw(idime,npoiw+4)=0.5*(coord(idime,lnods(ielem,4))
     .            +coord(idime,lnods(ielem,1)))
                coorw(idime,npoiw+5)=0.25*(coord(idime,lnods(ielem,1))
     .            +coord(idime,lnods(ielem,2))
     .            +coord(idime,lnods(ielem,3))
     .            +coord(idime,lnods(ielem,4)))
              end do
            else if (nnodw.eq.5) then
              do idime=1,2
                coorw(idime,npoiw+1)=0.25*(coord(idime,lnods(ielem,1))
     .            +coord(idime,lnods(ielem,2))
     .            +coord(idime,lnods(ielem,3))
     .            +coord(idime,lnods(ielem,4)))
              end do
            else if(nnodw.eq.16) then
              do idime=1,2
                x1=coord(idime,lnods(ielem,1))
                x2=coord(idime,lnods(ielem,2))
                x3=coord(idime,lnods(ielem,3))
                x4=coord(idime,lnods(ielem,4))
                coorw(idime,npoiw+1) =(2.0*x1+    x2)/3.0
                coorw(idime,npoiw+2) =(    x1+2.0*x2)/3.0
                coorw(idime,npoiw+3) =(2.0*x2+    x3)/3.0
                coorw(idime,npoiw+4) =(    x2+2.0*x3)/3.0
                coorw(idime,npoiw+5) =(2.0*x3+    x4)/3.0
                coorw(idime,npoiw+6) =(    x3+2.0*x4)/3.0
                coorw(idime,npoiw+7) =(2.0*x4+    x1)/3.0
                coorw(idime,npoiw+8) =(    x4+2.0*x1)/3.0
                coorw(idime,npoiw+9) =(2.0*x1+    x3)/3.0
                coorw(idime,npoiw+10)=(    x4+2.0*x2)/3.0
                coorw(idime,npoiw+11)=(2.0*x3+    x1)/3.0
                coorw(idime,npoiw+12)=(    x2+2.0*x4)/3.0
              end do
            end if
            if (nnodw.eq.3) then
              lnodw(ielem  ,1)=lnods(ielem,1)
              lnodw(ielem  ,2)=lnods(ielem,2)
              lnodw(ielem  ,3)=lnods(ielem,4)
              nelew=nelew+1
              lnodw(nelew  ,1)=lnods(ielem,2)
              lnodw(nelew  ,2)=lnods(ielem,3)
              lnodw(nelew  ,3)=lnods(ielem,4)
            else if(nnodw.eq.4) then
              lnodw(ielem  ,1)=lnods(ielem,1)
              lnodw(ielem  ,2)=npoiw+1
              lnodw(ielem  ,3)=npoiw+5
              lnodw(ielem  ,4)=npoiw+4
              lnodw(nelew+1,1)=npoiw+1
              lnodw(nelew+1,2)=lnods(ielem,2)
              lnodw(nelew+1,3)=npoiw+2
              lnodw(nelew+1,4)=npoiw+5
              lnodw(nelew+2,1)=npoiw+5
              lnodw(nelew+2,2)=npoiw+2
              lnodw(nelew+2,3)=lnods(ielem,3)
              lnodw(nelew+2,4)=npoiw+3
              lnodw(nelew+3,1)=npoiw+4
              lnodw(nelew+3,2)=npoiw+5
              lnodw(nelew+3,3)=npoiw+3
              lnodw(nelew+3,4)=lnods(ielem,4)
              nelew=nelew+3
              npoiw=npoiw+5
            else if(nnodw.eq.9) then
              lnodw(ielem  ,1)=lnods(ielem,1)
              lnodw(ielem  ,2)=lnods(ielem,2)
              lnodw(ielem  ,3)=lnods(ielem,3)
              lnodw(ielem  ,4)=lnods(ielem,4)
              lnodw(ielem  ,5)=npoiw+1
              lnodw(ielem  ,6)=npoiw+2
              lnodw(ielem  ,7)=npoiw+3
              lnodw(ielem  ,8)=npoiw+4
              lnodw(ielem  ,9)=npoiw+5
              npoiw=npoiw+5
            else if(nnodw.eq.5) then
              lnodw(ielem  ,1)=lnods(ielem,1)
              lnodw(ielem  ,2)=lnods(ielem,2)
              lnodw(ielem  ,3)=lnods(ielem,3)
              lnodw(ielem  ,4)=lnods(ielem,4)
              lnodw(ielem  ,5)=npoiw+1
              npoiw=npoiw+1
            else if(nnodw.eq.16) then
              lnodw(ielem  ,1) =lnods(ielem,1)
              lnodw(ielem  ,2) =lnods(ielem,2)
              lnodw(ielem  ,3) =lnods(ielem,3)
              lnodw(ielem  ,4) =lnods(ielem,4)
              lnodw(ielem  ,5) =npoiw+1
              lnodw(ielem  ,6) =npoiw+2
              lnodw(ielem  ,7) =npoiw+3
              lnodw(ielem  ,8) =npoiw+4
              lnodw(ielem  ,9) =npoiw+5
              lnodw(ielem  ,10)=npoiw+6
              lnodw(ielem  ,11)=npoiw+7
              lnodw(ielem  ,12)=npoiw+8
              lnodw(ielem  ,13)=npoiw+9
              lnodw(ielem  ,14)=npoiw+10
              lnodw(ielem  ,15)=npoiw+11
              lnodw(ielem  ,16)=npoiw+12
              npoiw=npoiw+12
            end if
c
c***  2D: NNODE = 9 --> NNODE = 3
c
          else if(nnode.eq.9) then
            if(nnodw.eq.4) then
              lnodw(ielem  ,1)=lnods(ielem,1)
              lnodw(ielem  ,2)=lnods(ielem,5)
              lnodw(ielem  ,3)=lnods(ielem,9)
              lnodw(ielem  ,4)=lnods(ielem,8)
              nelew=nelew+1
              lnodw(nelew  ,1)=lnods(ielem,5)
              lnodw(nelew  ,2)=lnods(ielem,2)
              lnodw(nelew  ,3)=lnods(ielem,6)
              lnodw(nelew  ,4)=lnods(ielem,9)
              nelew=nelew+1
              lnodw(nelew  ,1)=lnods(ielem,8)
              lnodw(nelew  ,2)=lnods(ielem,9)
              lnodw(nelew  ,3)=lnods(ielem,7)
              lnodw(nelew  ,4)=lnods(ielem,4)
              nelew=nelew+1
              lnodw(nelew  ,1)=lnods(ielem,9)
              lnodw(nelew  ,2)=lnods(ielem,6)
              lnodw(nelew  ,3)=lnods(ielem,3)
              lnodw(nelew  ,4)=lnods(ielem,7)
            else if(nnodw.eq.3) then
              lnodw(ielem  ,1)=lnods(ielem,1)
              lnodw(ielem  ,2)=lnods(ielem,5)
              lnodw(ielem  ,3)=lnods(ielem,9)
              lnodw(nelew+1,1)=lnods(ielem,1)
              lnodw(nelew+1,2)=lnods(ielem,9)
              lnodw(nelew+1,3)=lnods(ielem,8)
              lnodw(nelew+2,1)=lnods(ielem,8)
              lnodw(nelew+2,2)=lnods(ielem,9)
              lnodw(nelew+2,3)=lnods(ielem,7)
              lnodw(nelew+3,1)=lnods(ielem,8)
              lnodw(nelew+3,2)=lnods(ielem,7)
              lnodw(nelew+3,3)=lnods(ielem,4)
              lnodw(nelew+4,1)=lnods(ielem,5)
              lnodw(nelew+4,2)=lnods(ielem,2)
              lnodw(nelew+4,3)=lnods(ielem,6)
              lnodw(nelew+5,1)=lnods(ielem,5)
              lnodw(nelew+5,2)=lnods(ielem,6)
              lnodw(nelew+5,3)=lnods(ielem,9)
              lnodw(nelew+6,1)=lnods(ielem,9)
              lnodw(nelew+6,2)=lnods(ielem,6)
              lnodw(nelew+6,3)=lnods(ielem,3)
              lnodw(nelew+7,1)=lnods(ielem,9)
              lnodw(nelew+7,2)=lnods(ielem,3)
              lnodw(nelew+7,3)=lnods(ielem,7)
              nelew=nelew+7
            else if (nnodw.eq.6.or.nnodw.eq.7) then
              lnodw(ielem  ,1)=lnods(ielem,1)
              lnodw(ielem  ,2)=lnods(ielem,2)
              lnodw(ielem  ,3)=lnods(ielem,3)
              lnodw(ielem  ,4)=lnods(ielem,5)
              lnodw(ielem  ,5)=lnods(ielem,6)
              lnodw(ielem  ,6)=lnods(ielem,9)
              lnodw(nelew+1,1)=lnods(ielem,1)
              lnodw(nelew+1,2)=lnods(ielem,3)
              lnodw(nelew+1,3)=lnods(ielem,4)
              lnodw(nelew+1,4)=lnods(ielem,9)
              lnodw(nelew+1,5)=lnods(ielem,7)
              lnodw(nelew+1,6)=lnods(ielem,8)
              nelew=nelew+1
              if (nnodw.eq.7) then
                do idime=1,2
                  coorw(idime,npoiw+1)=-1.0/9.0*(coord(idime,lnodw(ielem,1))
     .              +coord(idime,lnodw(ielem,2))+coord(idime,lnodw(ielem,3)))
     .              +4.0/9.0*(coord(idime,lnodw(ielem,4))
     .              +coord(idime,lnodw(ielem,5))+coord(idime,lnodw(ielem,6)))
                  coorw(idime,npoiw+2)=-1.0/9.0*(coord(idime,lnodw(nelew,1))
     .              +coord(idime,lnodw(nelew,2))+coord(idime,lnodw(nelew,3)))
     .              +4.0/9.0*(coord(idime,lnodw(nelew,4))
     .              +coord(idime,lnodw(nelew,5))+coord(idime,lnodw(nelew,6)))
                end do
                lnodw(ielem,7)=npoiw+1
                lnodw(nelew,7)=npoiw+2
                npoiw=npoiw+2
              end if
            end if
c
c***  2D: NNODE = 16 --> NNODE = 10
c
          else if(nnode.eq.16) then
            if(nnodw.eq.10) then
              lnodw(ielem  , 1)=lnods(ielem, 1)
              lnodw(ielem  , 2)=lnods(ielem, 2)
              lnodw(ielem  , 3)=lnods(ielem, 3)
              lnodw(ielem  , 4)=lnods(ielem, 5)
              lnodw(ielem  , 5)=lnods(ielem, 6)
              lnodw(ielem  , 6)=lnods(ielem, 7)
              lnodw(ielem  , 7)=lnods(ielem, 8)
              lnodw(ielem  , 8)=lnods(ielem,15)
              lnodw(ielem  , 9)=lnods(ielem,13)
              lnodw(ielem  ,10)=lnods(ielem,14)
              lnodw(nelew+1, 1)=lnods(ielem, 1)
              lnodw(nelew+1, 2)=lnods(ielem, 3)
              lnodw(nelew+1, 3)=lnods(ielem, 4)
              lnodw(nelew+1, 4)=lnods(ielem,13)
              lnodw(nelew+1, 5)=lnods(ielem,15)
              lnodw(nelew+1, 6)=lnods(ielem, 9)
              lnodw(nelew+1, 7)=lnods(ielem,10)
              lnodw(nelew+1, 8)=lnods(ielem,11)
              lnodw(nelew+1, 9)=lnods(ielem,12)
              lnodw(nelew+1,10)=lnods(ielem,16)
              nelew=nelew+1
            end if
          end if                                            ! nnode
        else if(ndime.eq.3) then
c
c*** 3D: NNODE = 8 --> NNODE = 4 (same number of nodes)
c
          if(nnode.eq.8) then
            if(nnodw.eq.4) then
              lnodw(ielem,1) =lnods(ielem, 1)
              lnodw(ielem,2) =lnods(ielem, 3)
              lnodw(ielem,3) =lnods(ielem, 6)
              lnodw(ielem,4) =lnods(ielem, 2)
              nelew = nelew+1  
              lnodw(nelew,1) =lnods(ielem, 6)
              lnodw(nelew,2) =lnods(ielem, 3)
              lnodw(nelew,3) =lnods(ielem, 8)
              lnodw(nelew,4) =lnods(ielem, 7)
              nelew = nelew+1  
              lnodw(nelew,1) =lnods(ielem, 1)
              lnodw(nelew,2) =lnods(ielem, 6)
              lnodw(nelew,3) =lnods(ielem, 8)
              lnodw(nelew,4) =lnods(ielem, 5)
              nelew = nelew+1  
              lnodw(nelew,1) =lnods(ielem, 4)
              lnodw(nelew,2) =lnods(ielem, 8)
              lnodw(nelew,3) =lnods(ielem, 1)
              lnodw(nelew,4) =lnods(ielem, 3)
              nelew = nelew+1  
              lnodw(nelew,1) =lnods(ielem, 1)
              lnodw(nelew,2) =lnods(ielem, 3)
              lnodw(nelew,3) =lnods(ielem, 8)
              lnodw(nelew,4) =lnods(ielem, 6)
c
c*** 3D: NNODE = 8 --> NNODE = 20 (create new nodes)
c
            else if (nnodw.eq.20) then
              do inode=1,8
                lnodw(ielem,inode)=lnods(ielem,inode)
              end do
              do inode=9,20
                lnodw(ielem,inode) = npoiw+inode-8
              end do
              do inode=1,12
                do idime=1,ndime
                  coorw(idime,npoiw+inode)=0.0d0
                end do
              end do
              do idime=1,ndime                              ! 9
                coorw(idime,npoiw+1)=coorw(idime,npoiw+1)+0.5d0*
     .            (coord(idime,lnods(ielem,1))
     .            +coord(idime,lnods(ielem,2)))
              end do
              do idime=1,ndime                              ! 10
                coorw(idime,npoiw+2)=coorw(idime,npoiw+2)+0.5d0*
     .            (coord(idime,lnods(ielem,2))
     .            +coord(idime,lnods(ielem,3)))
              end do
              do idime=1,ndime                              ! 11
                coorw(idime,npoiw+3)=coorw(idime,npoiw+3)+0.5d0*
     .            (coord(idime,lnods(ielem,3))
     .            +coord(idime,lnods(ielem,4)))
              end do
              do idime=1,ndime                              ! 12
                coorw(idime,npoiw+4)=coorw(idime,npoiw+4)+0.5d0*
     .            (coord(idime,lnods(ielem,4))
     .            +coord(idime,lnods(ielem,1)))
              end do
              do idime=1,ndime                              ! 13
                coorw(idime,npoiw+5)=coorw(idime,npoiw+5)+0.5d0*
     .            (coord(idime,lnods(ielem,1))
     .            +coord(idime,lnods(ielem,5)))
              end do
              do idime=1,ndime                              ! 14
                coorw(idime,npoiw+6)=coorw(idime,npoiw+6)+0.5d0*
     .            (coord(idime,lnods(ielem,2))
     .            +coord(idime,lnods(ielem,6)))
              end do
              do idime=1,ndime                              ! 15
                coorw(idime,npoiw+7)=coorw(idime,npoiw+7)+0.5d0*
     .            (coord(idime,lnods(ielem,3))
     .            +coord(idime,lnods(ielem,7)))
              end do
              do idime=1,ndime                              ! 16
                coorw(idime,npoiw+8)=coorw(idime,npoiw+8)+0.5d0*
     .            (coord(idime,lnods(ielem,4))
     .            +coord(idime,lnods(ielem,8)))
              end do
              do idime=1,ndime                              ! 17
                coorw(idime,npoiw+9)=coorw(idime,npoiw+9)+0.5d0*
     .            (coord(idime,lnods(ielem,5))
     .            +coord(idime,lnods(ielem,6)))
              end do
              do idime=1,ndime                              ! 18
                coorw(idime,npoiw+10)=coorw(idime,npoiw+10)+0.5d0*
     .            (coord(idime,lnods(ielem,6))
     .            +coord(idime,lnods(ielem,7)))
              end do
              do idime=1,ndime                              ! 19
                coorw(idime,npoiw+11)=coorw(idime,npoiw+11)+0.5d0*
     .            (coord(idime,lnods(ielem,7))
     .            +coord(idime,lnods(ielem,8)))
              end do
              do idime=1,ndime                              ! 20
                coorw(idime,npoiw+12)=coorw(idime,npoiw+12)+0.5d0*
     .            (coord(idime,lnods(ielem,5))
     .            +coord(idime,lnods(ielem,8)))
              end do
              npoiw=npoiw+12
            else if(nnodw.eq.9) then
              npoiw=npoiw+1
              lnodw(ielem,1) =lnods(ielem, 1)
              lnodw(ielem,2) =lnods(ielem, 2)
              lnodw(ielem,3) =lnods(ielem, 3)
              lnodw(ielem,4) =lnods(ielem, 4)
              lnodw(ielem,5) =lnods(ielem, 5)
              lnodw(ielem,6) =lnods(ielem, 6)
              lnodw(ielem,7) =lnods(ielem, 7)
              lnodw(ielem,8) =lnods(ielem, 8)
              lnodw(ielem,9) =npoiw

              do idime=1,3
                coorw(idime,npoiw)=0.0d0
                do inode=1,4
                  coorw(idime,npoiw)=coorw(idime,npoiw)
     .              -0.25d0*coord(idime,lnods(ielem,inode))
     .              +0.50d0*coord(idime,lnods(ielem,inode+4))
                end do
              end do
            end if
c
c*** 3D: NNODE = 10 --> NNODE = 4 (same number of nodes)
c
          else if(nnode.eq.10) then
            if(nnodw.eq.4) then
              lnodw(ielem,1) = lnods(ielem,1)
              lnodw(ielem,2) = lnods(ielem,5)
              lnodw(ielem,3) = lnods(ielem,7)
              lnodw(ielem,4) = lnods(ielem,8)
              nelew = nelew+1
              lnodw(nelew,1) = lnods(ielem,2)
              lnodw(nelew,2) = lnods(ielem,6)
              lnodw(nelew,3) = lnods(ielem,5)
              lnodw(nelew,4) = lnods(ielem,9)
              nelew = nelew+1
              lnodw(nelew,1) = lnods(ielem,3)
              lnodw(nelew,2) = lnods(ielem,7)
              lnodw(nelew,3) = lnods(ielem,6)
              lnodw(nelew,4) = lnods(ielem,10)
              nelew = nelew+1
              lnodw(nelew,1) = lnods(ielem,4)
              lnodw(nelew,2) = lnods(ielem,10)
              lnodw(nelew,3) = lnods(ielem,9)
              lnodw(nelew,4) = lnods(ielem,8)
              nelew = nelew+1
              lnodw(nelew,1) = lnods(ielem,8)
              lnodw(nelew,2) = lnods(ielem,10)
              lnodw(nelew,3) = lnods(ielem,9)
              lnodw(nelew,4) = lnods(ielem,7)
              nelew = nelew+1
              lnodw(nelew,1) = lnods(ielem,8)
              lnodw(nelew,2) = lnods(ielem,9)
              lnodw(nelew,3) = lnods(ielem,5)
              lnodw(nelew,4) = lnods(ielem,7)
              nelew = nelew+1
              lnodw(nelew,1) = lnods(ielem,9)
              lnodw(nelew,2) = lnods(ielem,10)
              lnodw(nelew,3) = lnods(ielem,6)
              lnodw(nelew,4) = lnods(ielem,7)
              nelew = nelew+1
              lnodw(nelew,1) = lnods(ielem,9)
              lnodw(nelew,2) = lnods(ielem,6)
              lnodw(nelew,3) = lnods(ielem,5)
              lnodw(nelew,4) = lnods(ielem,7)
            end if
c
c*** 3D: NNODE = 4 --> NNODE = 5 (mini-element)
c
          else if(nnode.eq.4) then
            if(nnodw.eq.5) then
              npoiw=npoiw+1
              do idime=1,ndime
                coorw(idime,npoiw)=0.0
              end do
              do inode=1,4
                ipoin=lnods(ielem,inode)
                lnodw(ielem,inode)=ipoin
                do idime=1,ndime
                  coorw(idime,npoiw)=coorw(idime,npoiw)
     .              +coord(idime,ipoin)/4.0
                end do
              end do
              lnodw(ielem,5)=npoiw
c
c*** 3D: NNODE = 4 --> NNODE = 9
c
            else if (nnodw.eq.9) then
              do idime=1,3
                coorw(idime,npoiw+1)=0.5*(coord(idime,lnods(ielem,1))
     .            +coord(idime,lnods(ielem,2)))
                coorw(idime,npoiw+2)=0.5*(coord(idime,lnods(ielem,2))
     .            +coord(idime,lnods(ielem,3)))
                coorw(idime,npoiw+3)=0.5*(coord(idime,lnods(ielem,3))
     .            +coord(idime,lnods(ielem,4)))
                coorw(idime,npoiw+4)=0.5*(coord(idime,lnods(ielem,4))
     .            +coord(idime,lnods(ielem,1)))
                coorw(idime,npoiw+5)=0.25*(coord(idime,lnods(ielem,1))
     .            +coord(idime,lnods(ielem,2))
     .            +coord(idime,lnods(ielem,3))
     .            +coord(idime,lnods(ielem,4)))
              end do
              lnodw(ielem  ,1)=lnods(ielem,1)
              lnodw(ielem  ,2)=lnods(ielem,2)
              lnodw(ielem  ,3)=lnods(ielem,3)
              lnodw(ielem  ,4)=lnods(ielem,4)
              lnodw(ielem  ,5)=npoiw+1
              lnodw(ielem  ,6)=npoiw+2
              lnodw(ielem  ,7)=npoiw+3
              lnodw(ielem  ,8)=npoiw+4
              lnodw(ielem  ,9)=npoiw+5
              npoiw=npoiw+5
c
c*** 3D: NNODE = 4 --> NNODE = 10
c
            else if (nnodw.eq.10) then
              do idime=1,3
                coorw(idime,npoiw+1)=0.5*(coord(idime,lnods(ielem,1))
     .            +coord(idime,lnods(ielem,2)))
                coorw(idime,npoiw+2)=0.5*(coord(idime,lnods(ielem,2))
     .            +coord(idime,lnods(ielem,3)))
                coorw(idime,npoiw+3)=0.5*(coord(idime,lnods(ielem,3))
     .            +coord(idime,lnods(ielem,1)))
                coorw(idime,npoiw+4)=0.5*(coord(idime,lnods(ielem,1))
     .            +coord(idime,lnods(ielem,4)))
                coorw(idime,npoiw+5)=0.5*(coord(idime,lnods(ielem,2))
     .            +coord(idime,lnods(ielem,4)))
                coorw(idime,npoiw+6)=0.5*(coord(idime,lnods(ielem,3))
     .            +coord(idime,lnods(ielem,4)))
              end do
              do inode=1,4
                lnodw(ielem,  inode)=lnods(ielem,inode)
                lnodw(ielem,4+inode)=npoiw+inode
              end do
              lnodw(ielem, 9)=npoiw+5
              lnodw(ielem,10)=npoiw+6
              npoiw=npoiw+6
            end if                                          ! nnodw
          end if                                            ! nnode
        end if                                              ! ndime
      end do                                                ! ielem=1,nelem

      if(nelew.gt.mxelw) then
        write(6,'(a,i5)') 'Increase MXELW to ' , nelew
        stop
      end if
      if(npoiw.gt.mxpow) then
        write(6,'(a,i5)') 'Increase MXPOW to ' , npoiw
        stop
      end if

      end
c*************************************************************************     
      subroutine mescek(mxnow,mxelw,mxpow,coord,coorw,lnods,lnodw,
     .                  tempo,lpoiw,npoiw,nelew,nnodw)
      implicit     real*8 (a-h,o-z)
      common/contr/nin,nou,ndime,nelem,nnode,npoin
      dimension    coord(ndime,npoin), lnods(nelem,nnode),
     .             coorw(ndime,mxpow), lnodw(mxelw,mxnow),
     .             tempo(ndime,npoiw), lpoiw(npoiw)
       
      do ipoiw=1,npoiw
        lpoiw(ipoiw)=ipoiw
      end do

      z1=0.0
      z0=0.0
      do ipoiw=npoin+2,npoiw
        x1=coorw(1,ipoiw)
        y1=coorw(2,ipoiw)
        if(ndime.eq.3) z1=coorw(3,ipoiw)
        do jpoiw=npoin+1,ipoiw-1
          x0=coorw(1,jpoiw)
          y0=coorw(2,jpoiw)
          if(ndime.eq.3) z0=coorw(3,jpoiw)
          dist=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
          if(dist.lt.1.0e-7) then
            lpoiw(ipoiw)=-abs(lpoiw(jpoiw))
            do kpoiw=ipoiw+1,npoiw
              lpoiw(kpoiw)=lpoiw(kpoiw)-1
            end do
          end if
        end do
      end do

      npoif=0
      do ipoiw=1,npoiw
        if(lpoiw(ipoiw).gt.0) then
          npoif=npoif+1
          do idime=1,ndime
            tempo(idime,npoif)=coorw(idime,ipoiw)
          end do
        end if
      end do

      do ielew=1,nelew
        do inodw=1,nnodw
          lnodw(ielew,inodw)=abs(lpoiw(lnodw(ielew,inodw)))
        end do
      end do

      write(nou,'(a)') 'ELEMENTS'
      do ielem=1,nelew
        if(nnodw.eq.16) then
          write(nou,11) ielem,(lnodw(ielem,inode),inode=1,nnodw)
        else if (nnodw.gt.16) then
          write(nou,12) ielem,(lnodw(ielem,inode),inode=1,nnodw)
        else
          write(nou,10) ielem,(lnodw(ielem,inode),inode=1,nnodw)
        end if
      end do
      write(nou,'(a)') 'END_ELEMENTS'

      write(nou,'(a)') 'COORDINATES'
      do ipoin=1,npoif
        write(nou,20) ipoin,(tempo(idime,ipoin),idime=1,ndime)
      end do
      write(nou,'(a)') 'END_COORDINATES'
      
   10 format(1x,i5,10(1x,i6))
   11 format(1x,i3,16(1x,i6))      
   12 format(1x,i5,10(1x,i6)," /",/,6x,10(1x,i5))
   20 format(5x,i6,3(2x,e16.9))
      
      end
