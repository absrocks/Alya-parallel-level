!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    def_mpio.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   MPI-IO format definitions
!> @details This module defines several fields necessary to deal with the MPI-IO file format
!>          \verbatim
!>          - MPI-IO file extension
!>          - geometry file suffixes
!>          - some constants (header size, option size)
!>          - some header fixed values (magic number, version, format name)
!>          - structure mpio_header
!>          \endverbatim
!> @}
!-----------------------------------------------------------------------

module def_mpio

    use def_kintyp,                     only : ip,rp,lg,r1p


!-----------------------------------------------------------------------------------------------------------------------------
!                              Definitions
!-----------------------------------------------------------------------------------------------------------------------------

    character(9), parameter                 ::  mpio_ext    = '.mpio.bin'           ! MPIO file format extension

    character(9), parameter                 ::  post_ext    = '.post'               ! Post format extension

    character(6), parameter                 ::  coord_ext   = '-COORD',&
                                                ltype_ext   = '-LTYPE',&
                                                lnods_ext   = '-LNODS',&
                                                ltypb_ext   = '-LTYPB',&
                                                lelbo_ext   = '-LELBO',&
                                                lnodb_ext   = '-LNODB',&
                                                leinv_ext   = '-LEINV',&
                                                lninv_ext   = '-LNINV',&
                                                lbinv_ext   = '-LBINV',&
                                                lnoch_ext   = '-LNOCH',&
                                                lelch_ext   = '-LELCH',&
                                                lboch_ext   = '-LBOCH',&
                                                lmate_ext   = '-LMATE',&
                                                lmast_ext   = '-LMAST',&
                                                lesub_ext   = '-LESUB',&
                                                codbo_ext   = '-CODBO',&
                                                codno_ext   = '-CODNO',&
                                                leset_ext   = '-LESET',&
                                                lbset_ext   = '-LBSET',&
                                                lnset_ext   = '-LNSET',&
                                                field_ext   = '-XFIEL'

    integer(ip), parameter                  ::  header_size = 256                   ! Header size of MPI-IO file

    integer(ip), parameter                  ::  option_size = 10                    ! Option size

    integer(ip), parameter                  ::  value_count = 1                     ! Count of an integer or real value

    integer(ip), parameter                  ::  string_count = 8                    ! Count of a string

    integer, parameter                      ::  header_ip = 8

    integer, parameter                      ::  header_rp = 8

    integer(8), parameter                   ::  header_magic_number = 27093

    character(3), parameter                 ::  c5_f                = '00'//char(0)
    character(1), parameter                 ::  c7_f                = char(0)

    character(8),parameter                  ::  header_format       = 'MPIAL'//c5_f,&
                                                header_version      = 'V0004'//c5_f,&
                                                header_align_chars  = '00000'//c5_f,&
                                                header_option_init  = 'NONE0'//c5_f,&
                                                header_no_filter    = 'NOFIL'//c5_f,&
                                                header_asc_sorting  = 'ASCEN'//c5_f,&
                                                header_no_id        = 'NOID0'//c5_f

    integer(ip),        parameter           ::  PAR_MPIO_OFF = 0_ip
    integer(ip),        parameter           ::  PAR_MPIO_ON = 1_ip
    integer(ip),        parameter           ::  PAR_MPIO_FORCE_SEQ = 2_ip
    integer(ip),        parameter           ::  PAR_MPIO_NOT_IMPL = 3_ip
    integer(ip),        parameter           ::  PAR_MPIO_COMM_OFF = 0_ip
    integer(ip),        parameter           ::  PAR_MPIO_COMM_SFC = 1_ip
    integer(ip),        parameter           ::  PAR_MPIO_COMM_ALL = 2_ip

    integer(ip),        parameter           ::  IO_DISABLED = -1_ip
    integer(ip),        parameter           ::  IO_CLASSIC = 0_ip
    integer(ip),        parameter           ::  IO_MPIO_SEQ = 1_ip
    integer(ip),        parameter           ::  IO_MPIO_PAR = 2_ip
    integer(ip),        parameter           ::  IO_MPIO_POST = 3_ip


    TYPE mpio_header_options
        character(string_count), dimension(option_size) :: opt = (/header_option_init,&
                                                                   header_option_init,&
                                                                   header_option_init,&
                                                                   header_option_init,&
                                                                   header_option_init,&
                                                                   header_option_init,&
                                                                   header_option_init,&
                                                                   header_option_init,&
                                                                   header_option_init,&
                                                                   header_option_init/)
    END TYPE

    TYPE mpio_header
        ! fields contained in the binary file header
        integer(header_ip)                          :: magic_number = header_magic_number
        character(string_count)                     :: format       = header_format
        character(string_count)                     :: version      = header_version
        character(string_count)                     :: object
        character(string_count)                     :: dimension
        character(string_count)                     :: resultson
        character(string_count)                     :: type
        character(string_count)                     :: size
        character(string_count)                     :: par
        character(string_count)                     :: filter       = header_no_filter
        character(string_count)                     :: sorting      = header_asc_sorting
        character(string_count)                     :: id           = header_no_id
        character(string_count)                     :: align_chars  = header_align_chars
        integer(header_ip)                          :: columns
        integer(header_ip)                          :: lines
        integer(header_ip)                          :: ittim
        integer(header_ip)                          :: nsubd
        integer(header_ip)                          :: divi
        integer(header_ip)                          :: tag1
        integer(header_ip)                          :: tag2
        real(header_rp)                             :: time
        type(mpio_header_options)                   :: options
        ! fields that are not exported
        integer(ip)                                 :: item_size
        integer(8)                                  :: file_size
    END TYPE

    integer(8)                                      :: mpio_memor(2)

    integer(ip)                                     :: mpio_flag_enabled,                   &    ! Parallel IO using MPI-IO
                                                       mpio_flag_geometry,                  &    ! Parallel reading of geometry/mesh
                                                       mpio_flag_geometry_export,           &    ! Export mesh with MPI-IO format
                                                       mpio_flag_geometry_read_post,        &    ! Parallel reading of geometry/mesh
                                                       mpio_flag_post,                      &    ! Parallel IO for posts
                                                       mpio_flag_post_light,                &    ! Post-process only will export the mesh main files
                                                       mpio_flag_rst,                       &    ! Parallel IO for restarts
                                                       mpio_flag_autoromio,                 &    ! Parallel IO ROMIO auto configuration
                                                       mpio_flag_synchro,                   &    ! Parallel IO synchronous calls
                                                       mpio_val_asyncbuffer,                &    ! Parallel IO communicator strategy
                                                       mpio_flag_collective,                &    ! Parallel IO collective calls
                                                       mpio_flag_communicator,              &    ! Parallel IO communicator strategy
                                                       mpio_flag_all_par,                   &    ! Do not gather after reading
                                                       mpio_flag_post_merge,                &    ! Merge subdomains
                                                       mpio_val_merge_block                      ! Merging minimal block size


    integer(ip)                                     :: kfl_mpio_input = IO_CLASSIC
    integer(ip)                                     :: kfl_mpio_post = IO_CLASSIC
    integer(ip)                                     :: kfl_mpio_rst = IO_CLASSIC
    integer(ip)                                     :: kfl_mpio_export = PAR_MPIO_OFF

    public                                          ::  mpio_ext,&
                                                        coord_ext,&
                                                        ltype_ext,&
                                                        lnods_ext,&
                                                        ltypb_ext,&
                                                        lelbo_ext,&
                                                        lnodb_ext,&
                                                        leinv_ext,&
                                                        lninv_ext,&
                                                        lbinv_ext,&
                                                        lnoch_ext,&
                                                        lelch_ext,&
                                                        lboch_ext,&
                                                        lmate_ext,&
                                                        codbo_ext,&
                                                        codno_ext,&
                                                        leset_ext,&
                                                        lbset_ext,&
                                                        lnset_ext,&
                                                        field_ext,&
                                                        header_size,&
                                                        option_size,&
                                                        value_count,&
                                                        string_count,&
                                                        header_magic_number,&
                                                        header_format,&
                                                        header_version,&
                                                        header_align_chars,&
                                                        mpio_header,&
                                                        mpio_header_options,&
                                                        c5_f,&
                                                        c7_f,&
                                                        header_ip,&
                                                        header_rp,&
                                                        mpio_memor,&
                                                        PAR_MPIO_COMM_ALL,&
                                                        PAR_MPIO_COMM_OFF,&
                                                        PAR_MPIO_COMM_SFC,&
                                                        PAR_MPIO_FORCE_SEQ,&
                                                        PAR_MPIO_NOT_IMPL,&
                                                        PAR_MPIO_OFF,&
                                                        PAR_MPIO_ON,&
                                                        IO_CLASSIC,&
                                                        IO_DISABLED,&
                                                        IO_MPIO_PAR,&
                                                        IO_MPIO_POST,&
                                                        IO_MPIO_SEQ,&
                                                        mpio_flag_all_par,&
                                                        mpio_flag_autoromio,&
                                                        mpio_flag_collective,&
                                                        mpio_flag_communicator,&
                                                        mpio_flag_enabled,&
                                                        mpio_flag_geometry,&
                                                        mpio_flag_geometry_export,&
                                                        mpio_flag_geometry_read_post,&
                                                        mpio_flag_post,&
                                                        mpio_flag_rst,&
                                                        mpio_flag_synchro,&
                                                        mpio_flag_post_merge,&
                                                        mpio_val_asyncbuffer,&
                                                        mpio_val_merge_block,&
                                                        kfl_mpio_export,&
                                                        kfl_mpio_input,&
                                                        kfl_mpio_post,&
                                                        kfl_mpio_rst


end module
