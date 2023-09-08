subroutine ibm_blende(itask) 
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  imbou,cutim,intost,dtime,ittim,zeror
  use def_domain, only     :  ndime,nimbo
  use def_immbou, only     :  twall_ibm,nwaib

  implicit none

  integer(ip),intent(in)   :: itask
  integer(ip)              :: iimbo,iwaib,iboib,ipoib,idime,frame
  real(rp)                 :: angle
  character(20)            :: messa,mess1,mess2,mess3


  
  !open(unit=77,file='/Users/csamaniego/svn/blender/particles.py',status='old')

  if ( itask == 1 ) then
     open(unit=13,file="/Users/csamaniego/svn/Alya/blender/particles.py")
 
     write(13,'(A)') "import  bpy"
     write(13,'(A)') "bpy.data.scenes[0].frame_set(0)"
     do iimbo=1,nimbo        
        messa = intost(iimbo)

        write(13,'(A)')"mesh"//trim(messa)//" = bpy.data.meshes.new('mesh"//trim(messa)//"')"
        write(13,'(A)')"meshObj"//trim(messa)//" = bpy.data.objects.new('meshObj"//trim(messa)//"', mesh"//trim(messa)//" )"

        iboib = 1
        mess1 = intost(imbou(iimbo)%lnoib(1,iboib)-1_ip)
        mess2 = intost(imbou(iimbo)%lnoib(2,iboib)-1_ip)
        mess3 = intost(imbou(iimbo)%lnoib(3,iboib)-1_ip)        
        write(13,'(A)')"faces = (("//trim(mess1)//","//trim(mess2)//","//trim(mess3)//"), "//trim("\ ")

        do iboib = 2,imbou(iimbo) % nboib-1
           mess1 = intost(imbou(iimbo)%lnoib(1,iboib)-1_ip)
           mess2 = intost(imbou(iimbo)%lnoib(2,iboib)-1_ip)
           mess3 = intost(imbou(iimbo)%lnoib(3,iboib)-1_ip)
           write(13,'(A)') "("//trim(mess1)//","//trim(mess2)//","//trim(mess3)//"), "//trim("\ ")
        end do

        iboib = imbou(iimbo) % nboib
        mess1 = intost(imbou(iimbo)%lnoib(1,iboib)-1_ip)
        mess2 = intost(imbou(iimbo)%lnoib(2,iboib)-1_ip)
        mess3 = intost(imbou(iimbo)%lnoib(3,iboib)-1_ip)        
        write(13,'(A)') "("//trim(mess1)//","//trim(mess2)//","//trim(mess3)//"))"
        

        ipoib = 1
        write(13,'(A,e16.8E3,A,e16.8E3,A,e16.8E3,A)') "verts = ((", &
             imbou(iimbo)%cooin(1,ipoib),",",&
             imbou(iimbo)%cooin(2,ipoib),",",&
             imbou(iimbo)%cooin(3,ipoib),"), "//trim("\ ")

        do ipoib = 2,imbou(iimbo) % npoib-1
           write(13,'(A,e16.8E3,A,e16.8E3,A,e16.8E3,A)') "(",&
                imbou(iimbo)%cooin(1,ipoib),",",&
                imbou(iimbo)%cooin(2,ipoib),",",&
                imbou(iimbo)%cooin(3,ipoib),"), "//trim("\ ")
        end do

        ipoib = imbou(iimbo) % npoib
        write(13,'(A,e16.8E3,A,e16.8E3,A,e16.8E3,A)') "(",&
             imbou(iimbo)%cooin(1,ipoib),",",&
             imbou(iimbo)%cooin(2,ipoib),",",&
             imbou(iimbo)%cooin(3,ipoib),"))"
        
        write(13,'(A)') "mesh"//trim(messa)//".from_pydata( verts, (), faces )"
        write(13,'(A)') "bpy.context.scene.objects.link( meshObj"//trim(messa)//" )"

        write(13,'(A)') "bpy.data.objects['meshObj"//trim(messa)//"'].rotation_mode='QUATERNION'"        
        do idime=1,ndime+1
           mess1 = intost(idime-1_ip)
           write(13,'(A,e16.8E3)') "bpy.data.objects['meshObj"//trim(messa)//"'].rotation_quaternion["//trim(mess1)//"]=",&                
                imbou(iimbo)%quate(idime,1)
        end do
        write(13,'(A)') "bpy.data.objects['meshObj"//trim(messa)//"'].keyframe_insert('rotation_quaternion')"                
        do idime=1,ndime
           mess1 = intost(idime-1_ip)
           write(13,'(A,e16.8E3)') "bpy.data.objects['meshObj"//trim(messa)//"'].location["//trim(mess1)//"]=",&
                imbou(iimbo)%posil(idime,1)
        end do
        write(13,'(A)') "bpy.data.objects['meshObj"//trim(messa)//"'].keyframe_insert('location')"                
     end do

     do iwaib=1,nwaib
        messa = intost(iwaib + nimbo)

        write(13,'(A)')"mesh"//trim(messa)//" = bpy.data.meshes.new('mesh"//trim(messa)//"')"
        write(13,'(A)')"meshObj"//trim(messa)//" = bpy.data.objects.new('meshObj"//trim(messa)//"', mesh"//trim(messa)//" )"

        iboib = 1
        mess1 = intost(twall_ibm(iwaib) % lnodb(1,iboib)-1_ip)
        mess2 = intost(twall_ibm(iwaib) % lnodb(2,iboib)-1_ip)
        mess3 = intost(twall_ibm(iwaib) % lnodb(3,iboib)-1_ip)        
        write(13,'(A)')"faces = (("//trim(mess1)//","//trim(mess2)//","//trim(mess3)//"), "//trim("\ ")

        do iboib = 2,twall_ibm(iwaib) % nboun-1
           mess1 = intost(twall_ibm(iwaib) % lnodb(1,iboib)-1_ip)
           mess2 = intost(twall_ibm(iwaib) % lnodb(2,iboib)-1_ip)
           mess3 = intost(twall_ibm(iwaib) % lnodb(3,iboib)-1_ip)
           write(13,'(A)') "("//trim(mess1)//","//trim(mess2)//","//trim(mess3)//"), "//trim("\ ")
        end do

        iboib = twall_ibm(iwaib) % nboun
        mess1 = intost(twall_ibm(iwaib) % lnodb(1,iboib)-1_ip)
        mess2 = intost(twall_ibm(iwaib) % lnodb(2,iboib)-1_ip)
        mess3 = intost(twall_ibm(iwaib) % lnodb(3,iboib)-1_ip)        
        write(13,'(A)') "("//trim(mess1)//","//trim(mess2)//","//trim(mess3)//"))"
        

        ipoib = 1
        write(13,'(A,e16.8E3,A,e16.8E3,A,e16.8E3,A)') "verts = ((", &
             twall_ibm(iwaib) % coord(1,ipoib),",",&
             twall_ibm(iwaib) % coord(2,ipoib),",",&
             twall_ibm(iwaib) % coord(3,ipoib),"), "//trim("\ ")

        do ipoib = 2,twall_ibm(iwaib) % npoin-1
           write(13,'(A,e16.8E3,A,e16.8E3,A,e16.8E3,A)') "(",&
                twall_ibm(iwaib) % coord(1,ipoib),",",&
                twall_ibm(iwaib) % coord(2,ipoib),",",&
                twall_ibm(iwaib) % coord(3,ipoib),"), "//trim("\ ")
        end do

        ipoib = twall_ibm(iwaib) % npoin
        write(13,'(A,e16.8E3,A,e16.8E3,A,e16.8E3,A)') "(",&
             twall_ibm(iwaib) % coord(1,ipoib),",",&
             twall_ibm(iwaib) % coord(2,ipoib),",",&
             twall_ibm(iwaib) % coord(3,ipoib),"))"
        
        write(13,'(A)') "mesh"//trim(messa)//".from_pydata( verts, (), faces )"
        write(13,'(A)') "bpy.context.scene.objects.link( meshObj"//trim(messa)//" )"     
     end do

     close(13)
  elseif ( itask == 2 ) then

     messa = intost(int(real(ittim,rp)/10.0_rp,ip))
     if ( mod(ittim,10_ip) == 1) then
        open(unit=14,file="/Users/csamaniego/svn/Alya/blender/movement"//trim(messa)//".py")
        write(14,'(A)') "import bpy"
        close(14)
     end if

     ! Previous time step
      call runend('ibm_blende. comentarie open que daba error')
!     open(unit=14,file="/Users/csamaniego/svn/Alya/blender/movement"//trim(messa)//".py",access="append")

     frame = int(cutim/ 1.0e-3_rp,ip) + 1_ip      
     messa = intost(frame)

     write(14,'(A)')  "bpy.data.scenes[0].frame_set("//trim(messa)//")"
     do iimbo=1,nimbo      
        messa = intost(iimbo)

        angle = 2.0_rp*acos(imbou(iimbo)%quate(1,1))
        if ( abs(cos(angle/2.0_rp)) > zeror ) then
           write(14,'(A,e16.8E3)') "bpy.data.objects['meshObj"//trim(messa)//"'].rotation_quaternion[0]=",&                
                imbou(iimbo)%quate(1,1)*(1.0_rp/cos(angle/2.0_rp))*cos(angle/4.0_rp)
        else
           write(14,'(A,e16.8E3)') "bpy.data.objects['meshObj"//trim(messa)//"'].rotation_quaternion[0]=",&                
                imbou(iimbo)%quate(1,1)
        end if
        do idime=2,ndime+1
              
           mess1 = intost(idime-1_ip)
           if ( abs(sin(angle/2.0_rp)) > zeror ) then
              write(14,'(A,e16.8E3)') "bpy.data.objects['meshObj"//trim(messa)//"'].rotation_quaternion["//trim(mess1)//"]=",&                
                   imbou(iimbo)%quate(idime,1)*(1.0_rp/sin(angle/2.0_rp))*sin(angle/4.0_rp)                
           else
              write(14,'(A,e16.8E3)') "bpy.data.objects['meshObj"//trim(messa)//"'].rotation_quaternion["//trim(mess1)//"]=",&                
                   imbou(iimbo)%quate(idime,1)
           end if
        end do
        write(14,'(A)') "bpy.data.objects['meshObj"//trim(messa)//"'].keyframe_insert('rotation_quaternion')"

        do idime=1,ndime
           mess1 = intost(idime-1_ip)
           write(14,'(A,e16.8E3)') "bpy.data.objects['meshObj"//trim(messa)//"'].location["//trim(mess1)//"]=",&
                imbou(iimbo)%posil(idime,1)
        end do
        write(14,'(A)') "bpy.data.objects['meshObj"//trim(messa)//"'].keyframe_insert('location')"

     end do
     close(14)
  end  if
              
  

end subroutine ibm_blende
