module mod_usenetcdf

  use netCDF
  use KindType

  implicit none
  

  interface readv_netcdf
     module procedure readv_netcdf_scalar, &
          readv_netcdf_array, &
          readv_netcdf_matrix, &
          readv_netcdf_3d
     
          
  end interface readv_netcdf
  
  public  :: readv_netcdf
  public  :: read_dim_netcdf
  private :: readv_netcdf_array, &
       readv_netcdf_matrix, &
       readv_netcdf_3d
  
contains
  
  subroutine read_dim_netcdf(fileID, dim_name, nsize)
    ! reads the content in dimension name dim_name, and loads it in  nsize
    ! 
    !-------------------------------------------------------------------------
    
    implicit none
    integer(ip), intent(in)                     :: fileID
    character(len=s_name), intent(inout)        :: dim_name
    integer(ip), intent(out)                    :: nsize
    
    !    local variables
    integer(ip)               ::DimID
    
    if( nf90_inq_dimid(fileID,dim_name,DimID) /= 0 ) &
         call runend('Error in nf90_inq_dimid '//TRIM(dim_name))
    if( nf90_inquire_dimension(fileID,DimID,dim_name,nsize) /= 0 ) &
         call runend('Error in nf90_inquire_dimension '//TRIM(dim_name))
    
  end subroutine read_dim_netcdf
  
  
  subroutine readv_netcdf_3d(fileID, vname, varia, size1, size2, size3)
    ! reads the content in variable name vname, and loads it in varia(i,j,k)
    ! 
    !-------------------------------------------------------------------------
    
    implicit none
    integer(ip), intent(in)                     :: fileID
    character(len=s_name), intent(in)           :: vname
    integer(ip), intent(in)                     :: size1,size2,size3
    real(rp), intent(out)                       :: varia(size1,size2,size3)

    !    local variables
    integer(ip)               ::VarID
    
    if( nf90_inq_varid(fileID,vname,VarID) /= 0) &
         call runend('Error getting VarID for '//TRIM(vname))
    if( nf90_get_var  (fileID,VarID,varia,start=(/1,1,1/),count=(/size1,size2,size3/)) /= 0) &
         call runend('Error reading variable'//TRIM(vname))
    
  end subroutine readv_netcdf_3d

  subroutine readv_netcdf_matrix(fileID, vname, varia, size1, size2)
    ! reads the content in variable name vname, and loads it in varia(i,j,k)
    ! 
    !-------------------------------------------------------------------------
    
    implicit none
    integer(ip), intent(in)                     :: fileID
    character(len=s_name), intent(in)           :: vname
    integer(ip), intent(in)                     :: size1,size2
    real(rp), intent(out)                       :: varia(size1,size2)

    !    local variables
    integer(ip)               ::VarID
    
    if( nf90_inq_varid(fileID,vname,VarID) /= 0) &
         call runend('Error getting VarID for '//TRIM(vname))
    if( nf90_get_var  (fileID,VarID,varia,start=(/1,1/),count=(/size1,size2/)) /= 0) &
         call runend('Error reading variable'//TRIM(vname))
    
  end subroutine readv_netcdf_matrix

  subroutine readv_netcdf_array(fileID, vname, varia, size1)
    ! reads the content in variable name vname, and loads it in varia(i,j,k)
    ! 
    !-------------------------------------------------------------------------
    
    implicit none
    integer(ip), intent(in)                     :: fileID
    character(len=s_name), intent(in)           :: vname
    integer(ip), intent(in)                     :: size1
    real(rp), intent(out)                       :: varia(size1)

    !    local variables
    integer(ip)               ::VarID
    
    if( nf90_inq_varid(fileID,vname,VarID) /= 0) &
         call runend('Error getting VarID for '//TRIM(vname))
    if( nf90_get_var  (fileID,VarID,varia,start=(/1/),count=(/size1/)) /= 0) &
         call runend('Error reading variable'//TRIM(vname))
    
  end subroutine readv_netcdf_array
  subroutine readv_netcdf_scalar(fileID, vname, varia)
    ! reads the content in  variable name vname, and loads it in varia(i,j,k)
    ! 
    !-------------------------------------------------------------------------
    
    implicit none
    integer(ip), intent(in)                     :: fileID
    character(len=s_name), intent(in)           :: vname
    real(rp), intent(out)                       :: varia

    !    local variables
    integer(ip)               ::VarID
    
    if( nf90_inq_varid(fileID,vname,VarID) /= 0) &
         call runend('Error getting VarID for '//TRIM(vname))
    if( nf90_get_var  (fileID,VarID,varia) /= 0) &
         call runend('Error reading variable'//TRIM(vname))
    
  end subroutine readv_netcdf_scalar
  
end module mod_usenetcdf

