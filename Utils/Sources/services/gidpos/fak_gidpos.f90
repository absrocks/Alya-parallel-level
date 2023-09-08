subroutine GID_OPENPOSTRESULTFILE(fil_postp,itask)
  use def_kintyp
  character(150) :: fil_postp
  integer(ip)    :: itask
  call runend('SERVICE GIDPOS WAS NOT COMPILED')
  return
end subroutine GID_OPENPOSTRESULTFILE

subroutine GID_CLOSEPOSTRESULTFILE
  return
end subroutine GID_CLOSEPOSTRESULTFILE

subroutine GID_BEGINSCALARRESULT(wopos,state,ttime,ittim,char1,char2,char3)
  use def_kintyp
  character(*) :: wopos
  integer(ip)  :: ittim
  real(rp)     :: ttime
  character(8) :: state
  character(4) :: char1,char2,char3
  return
end subroutine GID_BEGINSCALARRESULT

subroutine GID_BEGINVECTORRESULT(wopos,state,ttime,ittim,char1,char2,char3,char4)
  use def_kintyp
  character(*)  :: wopos
  integer(ip)   :: ittim
  real(rp)      :: ttime
  character(8)  :: state
  character(4)  :: char1,char2,char3,char4
  return
end subroutine GID_BEGINVECTORRESULT


subroutine GID_WRITESCALAR(ipoin,bridge)
  use def_kintyp
  real(rp)    :: bridge
  integer(ip) :: ipoin
  return
end subroutine GID_WRITESCALAR

subroutine GID_WRITEVECTOR(ipoin,bridg1,bridg2,bridg3)
  use def_kintyp
  real(rp)    :: bridg1,bridg2,bridg3
  integer(ip) :: ipoin
  return
end subroutine GID_WRITEVECTOR

subroutine GID_ENDRESULT()
end subroutine GID_ENDRESULT

subroutine GiD_BeginCoordinates()
end subroutine GiD_BeginCoordinates

subroutine GiD_EndCoordinates()
end subroutine GiD_EndCoordinates

subroutine GiD_BeginElements()
end subroutine GiD_BeginElements

subroutine GiD_EndElements()
end subroutine GiD_EndElements

subroutine GiD_BeginMesh()
end subroutine GiD_BeginMesh

subroutine GiD_EndMesh()
end subroutine GiD_EndMesh

subroutine GiD_WriteCoordinates2D(ipoin,bridg1,bridg2)
  use def_kintyp
  real(rp)     :: bridg1(*),bridg2(*)
  integer(ip)  :: ipoin
  ipoin = ipoin
  bridg1(1)  = bridg1(1) 
  bridg2(1)  = bridg2(1) 
end subroutine GiD_WriteCoordinates2D

subroutine GiD_WriteCoordinates(ipoin,bridg1,bridg2,bridg3)
  use def_kintyp
  real(rp) :: bridg1(*),bridg2(*),bridg3(*)
  integer(ip)  :: ipoin
end subroutine GiD_WriteCoordinates

subroutine GiD_WriteElementMat(ipoin,bridg1)
  use def_kintyp
  integer(ip) :: ipoin,bridg1(*)
end subroutine GiD_WriteElementMat

subroutine GiD_ClosePostMeshFile()
end subroutine GiD_ClosePostMeshFile

subroutine GiD_OpenPostMeshFile(fil_postp,itask)
  use def_kintyp
  character(150) :: fil_postp
  integer(ip)    :: itask
  call runend('SERVICE GIDPOS WAS NOT COMPILED')
end subroutine GiD_OpenPostMeshFile
