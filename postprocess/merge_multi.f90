! calculate the mean and standard deviation of NEMO-MEDUSA output
program merge
use netcdf
use mpi
implicit none
integer, parameter :: Ne = 30
! integer :: nx = 362, ny = 332, nz = 75, Ne = 30
! dimension of the variables
integer :: nt, nx, ny, nz = 1
! MPI related variables
integer :: MPIerr
logical :: iniflag
integer :: rank
integer :: npes
integer :: chunksize, remainder, displs
integer :: nfiles, i_file
CHARACTER(len=4) :: nfiles_c, i_file_index
! input arguments
CHARACTER(len=1) :: doLogC
CHARACTER(len=1024) :: BaseDir
CHARACTER(len=1024) :: filename
! switch for log10 transformation
logical             :: do_logspace

! netCDF4 variables
integer :: ncid
integer :: ierr
integer :: dimid
integer :: nVariables, nAttributes
CHARACTER(LEN=20), dimension(8) :: exclude_varnames
! added variables; at the moment, it only applies to total nitrogen and chlorophyll data
integer :: n_new_vars
logical :: has_nitrogen

! initialise number of added variables
n_new_vars = 0
has_nitrogen = .false.
! no need to calculate mean and std of these variables
exclude_varnames = [CHARACTER(LEN=20) :: 'nav_lat', 'nav_lon', 'deptht', 'deptht_bounds', &
                     'time_centered', 'time_centered_bounds', &
                     'time_counter', 'time_counter_bounds']

! initialise the MPI processors
CALL MPI_Initialized(iniflag, MPIerr)
IF (.not.iniflag) CALL MPI_Init(MPIerr)
CALL MPI_Comm_size(MPI_COMM_WORLD, npes, MPIerr)
CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, MPIerr)
! input argument
! BaseDir: the output file directory
CALL get_command_argument(1, BaseDir)
! filename: base filename of the output file 
CALL get_command_argument(2, filename)
! number of subdomains being processed
CALL get_command_argument(3, nfiles_c)
! switch for converting ensemble to log10 space
CALL get_command_argument(4, doLogC)
! convert string to boolean type
do_logspace = .false.
if (doLogC == 'T') do_logspace = .true.
! get nfiles
read(nfiles_c, *) nfiles

! decompose subdomains into processors
chunksize = nfiles / npes + 1
remainder = mod(nfiles, npes)
displs = rank * chunksize
if (rank >= remainder) then
   chunksize = nfiles / npes
   displs = rank * chunksize + remainder
end if


if (rank == 0) print *, 'compute mean and std in the log space', do_logspace

print *, 'rank', rank, 'start', displs, 'end', displs+chunksize - 1
do i_file = displs, displs + chunksize - 1
   write(i_file_index, '(I4.4)') i_file
   ! open the output file which is created by copying file from one ensmeble member
   if (do_logspace) then
      print *, 'Opening ',  trim(BaseDir)//'/'//trim(filename)//'-log_'//i_file_index//'.nc'
      call check( nf90_open(trim(BaseDir)//'/'//trim(filename)//'-log_'//i_file_index//'.nc', nf90_write, ncid) )
   else
      print *, 'Opening ',  trim(BaseDir)//'/'//trim(filename)//'_'//i_file_index//'.nc'
      call check( nf90_open(trim(BaseDir)//'/'//trim(filename)//'_'//i_file_index//'.nc', nf90_write, ncid) )
   endif
   ! get the total number of variables and attributes
   call check( nf90_inquire(ncid, nVariables=nVariables, nAttributes=nAttributes) )
   ! get the time and space dimension in each file
   call check( nf90_inq_dimid(ncid, 'time_counter', dimid) )
   call check( nf90_inquire_dimension(ncid, dimid, len=nt) )
   call check( nf90_inq_dimid(ncid, 'x', dimid) )
   call check( nf90_inquire_dimension(ncid, dimid, len=nx) )
   call check( nf90_inq_dimid(ncid, 'y', dimid) )
   call check( nf90_inquire_dimension(ncid, dimid, len=ny) )
   ierr = nf90_inq_dimid(ncid, 'deptht', dimid)
   if (ierr == nf90_noerr) &
      call check( nf90_inquire_dimension(ncid, dimid, len=nz) )

   call add_total_nitrogen_and_chlorophyll_variable_to_output
   if (has_nitrogen .and. i_file == displs) n_new_vars = n_new_vars + 2
   ! add ensemble standard deviation and mean information
   ! delete uuid and add new time stamp
   call redefine_metadata()
   ! put calculated ensemble standard deviation and mean to the netCDF file
   call calculate_mean_std()

   print *, rank, 'Closing ', trim(BaseDir)//'/'//trim(filename)//'_'//i_file_index//'.nc'
   call check( nf90_close(ncid) )
end do

! end MPI program gracefully
call MPI_Finalize(MPIerr)

contains
   !> Check status of NC operation
   !!   
   subroutine check(status)
      ! *** Aruments ***
      integer, intent ( in) :: status   ! Reading status
      if(status /= nf90_noerr) then 
         print *, rank, trim(nf90_strerror(status))
         call MPI_ABORT(MPI_COMM_WORLD, -1, MPIerr)
      end if
   end subroutine check

   function iso8601() result(datetime)
      !! Returns current date and time in ISO 8601 format.
      !! from https://cyber.dabamos.de/programming/modernfortran/date-and-time.html
      character(len=*), parameter :: ISO8601_FMT = '(i4, 2("-", i2.2), "T", 2(i0.2, ":"), i0.2, ".", i0.3, a, ":", a)'
      character(len=29) :: datetime
      character(len=5)  :: zone
      integer           :: dt(8)

      call date_and_time(values=dt, zone=zone)

      write (datetime, ISO8601_FMT) dt(1), dt(2), dt(3), dt(5), dt(6), &
          dt(7), dt(8), zone(1:3), zone(4:5)
   end function iso8601

   subroutine add_total_nitrogen_and_chlorophyll_variable_to_output()
      integer :: varid     ! varid of size classes
      integer :: varid_new ! varid of nitrogen/chlorophyll variable
      integer :: xtype     ! datatype
      integer :: ndims     ! number of dimensions of the variable
      integer :: nAtts     ! number of attributes of each variable
      integer :: i_att     ! variable attribute counter
      integer, allocatable :: dimids(:)

      CHARACTER(LEN=1024) :: attname
      CHARACTER(LEN=1024) :: attval

      ! check if the file contains required variables and get the varid of 'PHD'
      ierr = nf90_inq_varid(ncid, 'PHD', varid)
      ! if ierr != nf90_noerr, this means variable PHD is not in this file
      ! we do not need to add variables
      if (ierr /= nf90_noerr) return
      print *, rank, 'add new variable nitrogen and chlorophyll data'
      ! open the mode to redefine metadata
      call check( nf90_redef(ncid) )
      ! now create nitrogen variable
      ! get variable dimensions and datatype
      call check( nf90_inquire_variable(ncid, varid, xtype=xtype, ndims=ndims, nAtts=nAtts) )
      if (allocated(dimids)) deallocate(dimids)
      allocate(dimids(ndims))
      call check( nf90_inquire_variable(ncid, varid, dimids=dimids) )
      ! define variable
      call check( nf90_def_var(ncid, 'nitrogen', xtype, dimids, varid_new) )
      do i_att = 1, nAtts
         call check( nf90_inq_attname(ncid, varid, i_att, attname) )
         ! copy attribute to newly created variable
         call check( nf90_copy_att(ncid, varid, trim(attname), ncid, varid_new) )
         ! a new long name should be given
         if (trim(attname) == 'long_name') then
            attval = 'pelagic phytoplankton nitrogen'
            if (do_logspace) attval = trim(attval)//'_in log10 space'
            call check( nf90_put_att(ncid, varid_new, trim(attname), trim(attval)) )
         end if
      end do
      ! as chlorophyll has the same dimension and datatype as PHD
      ! we can define chlorophyll as well
      call check( nf90_def_var(ncid, 'chlorophyll', xtype, dimids, varid_new) )
      do i_att = 1, nAtts
         call check( nf90_inq_attname(ncid, varid, i_att, attname) )
         ! copy attribute to newly created variable
         call check( nf90_copy_att(ncid, varid, trim(attname), ncid, varid_new) )
         ! a new long name should be given
         if (trim(attname) == 'long_name') then
            attval = 'pelagic phytoplankton chlorophyll'
            if (do_logspace) attval = trim(attval)//'_in log10 space'
            call check( nf90_put_att(ncid, varid_new, trim(attname), trim(attval)) )
         end if
      end do
      ! stop metadata redefinition
      call check( nf90_enddef(ncid) )
      ! clean allocated array
      if (allocated(dimids)) deallocate(dimids)
      ! two new vars are added
      has_nitrogen = .true.
   end subroutine add_total_nitrogen_and_chlorophyll_variable_to_output

   subroutine calculate_total_phytoplankton_mean(varname, xtype, &
                                                 missing_value_f, missing_value_d, &
                                                 val_4d_f, val_4d_d, intmd_4d_f, intmd_4d_d, &
                                                 mean_4d_f, mean_4d_d)
      CHARACTER(LEN=*), intent(in) :: varname
      integer, intent(in)          :: xtype
      real(4), intent(in)          :: missing_value_f
      real(8), intent(in)          :: missing_value_d
      real(4), intent(inout)       :: val_4d_f(:, :, :, :)
      real(8), intent(inout)       :: val_4d_d(:, :, :, :)
      real(4), intent(inout)       :: intmd_4d_f(:, :, :, :)
      real(8), intent(inout)       :: intmd_4d_d(:, :, :, :)
      real(4), intent(inout)       :: mean_4d_f(:, :, :, :)
      real(8), intent(inout)       :: mean_4d_d(:, :, :, :)

      integer :: i, j
      integer :: ncid_ens
      integer :: varid_ens
      CHARACTER(LEN=3), dimension(2) :: sub_varnames
      ! define varnames based on varname
      if (trim(varname) == 'chlorophyll') then
         sub_varnames(1) = 'CHD'
         sub_varnames(2) = 'CHN'
      else if (trim(varname) == 'nitrogen') then
         sub_varnames(1) = 'PHD'
         sub_varnames(2) = 'PHN'
      else
         print *, rank, 'Incorrect varname', trim(varname),' is used to calculate total phytoplankton values'
      end if
      ! 
      mean_4d_f = 0.
      mean_4d_d = 0.
      ! 
      do i = 1, Ne
         ! open ensemble member file
         call check( nf90_open(trim(BaseDir)//'/ensemble_'//trim(str(i))//'/'//trim(filename)//'_'//i_file_index//'.nc', &
                               nf90_nowrite, ncid_ens) )
         intmd_4d_f = 0.
         intmd_4d_d = 0.
         do j  = 1, 2
            val_4d_f = 0.
            val_4d_d = 0.
            call check( nf90_inq_varid(ncid_ens, trim(sub_varnames(j)), varid_ens) )
            if (xtype == NF90_float) &
               call check( nf90_get_var(ncid_ens, varid_ens, val_4d_f(:, :, :, :), &
                           start=[1, 1, 1, 1], count = [nx, ny, nz, nt]) )
            if (xtype == NF90_double) &
               call check( nf90_get_var(ncid_ens, varid_ens, val_4d_d(:, :, :, :), &
                           start=[1, 1, 1, 1], count = [nx, ny, nz, nt]) )
            intmd_4d_f = intmd_4d_f + val_4d_f
            intmd_4d_d = intmd_4d_d + val_4d_d
         end do
         ! compute the average value
         if (do_logspace) then
            where (intmd_4d_f > 0. .and. abs(intmd_4d_f - missing_value_f) > 1e-6)
               intmd_4d_f = log10(intmd_4d_f)
            else where (abs(intmd_4d_f - missing_value_f) < 1e-6)
               intmd_4d_f = missing_value_f
            else where
               intmd_4d_f = -14
            end where

            where (intmd_4d_d > 0..and. abs(intmd_4d_d - missing_value_d) > 1e-6)
               intmd_4d_d = log10(intmd_4d_d)
            else where (abs(intmd_4d_d - missing_value_d) < 1e-6)
               intmd_4d_d = missing_value_d
            else where
               intmd_4d_d = -14
            end where
         end if
         ! calculate mean value
         where (abs(intmd_4d_f - missing_value_f) < 1e-6)
            mean_4d_f = missing_value_f
         else where
            mean_4d_f = mean_4d_f + intmd_4d_f/Ne
         end where

         where (abs(intmd_4d_d - missing_value_d) < 1e-6)
            mean_4d_d = missing_value_d
         else where
            mean_4d_d = mean_4d_d + intmd_4d_d/Ne
         end where
         call check ( nf90_close (ncid_ens) )
      end do
      where (mean_4d_f > missing_value_f)
         mean_4d_f = missing_value_f
      end where
      where (mean_4d_d > missing_value_d)
         mean_4d_d = missing_value_d
      end where
   end subroutine calculate_total_phytoplankton_mean

   subroutine calculate_total_phytoplankton_std(varname, xtype, &
                                                 missing_value_f, missing_value_d, &
                                                 mean_4d_f, mean_4d_d, &
                                                 val_4d_f, val_4d_d, intmd_4d_f, intmd_4d_d, &
                                                 std_4d_f, std_4d_d)
      CHARACTER(LEN=*), intent(in) :: varname
      integer, intent(in)          :: xtype
      real(4), intent(in)          :: missing_value_f
      real(8), intent(in)          :: missing_value_d
      real(4), intent(in)          :: mean_4d_f(:, :, :, :)
      real(8), intent(in)          :: mean_4d_d(:, :, :, :)
      real(4), intent(inout)       :: val_4d_f(:, :, :, :)
      real(8), intent(inout)       :: val_4d_d(:, :, :, :)
      real(4), intent(inout)       :: intmd_4d_f(:, :, :, :)
      real(8), intent(inout)       :: intmd_4d_d(:, :, :, :)
      real(4), intent(out)         :: std_4d_f(:, :, :, :)
      real(8), intent(out)         :: std_4d_d(:, :, :, :)

      integer :: i, j
      integer :: ncid_ens
      integer :: varid_ens
      CHARACTER(LEN=3), dimension(2) :: sub_varnames

      ! define varnames based on varname
      if (trim(varname) == 'chlorophyll') then
         sub_varnames(1) = 'CHD'
         sub_varnames(2) = 'CHN'
      else if (trim(varname) == 'nitrogen') then
         sub_varnames(1) = 'PHD'
         sub_varnames(2) = 'PHN'
      else
         print *, rank, 'Incorrect varname', trim(varname),' is used to calculate total phytoplankton standard deviation'
      end if
      ! 
      ! get ensemble std
      std_4d_f = 0.
      std_4d_d = 0.

      do i = 1, Ne
         call check( nf90_open(trim(BaseDir)//'/ensemble_'//trim(str(i))//'/'//trim(filename)//'_'//i_file_index//'.nc', &
                               nf90_nowrite, ncid_ens) )
         intmd_4d_f = 0.
         intmd_4d_d = 0.
         do j = 1, 2
            val_4d_f = 0.
            val_4d_d = 0.
            call check( nf90_inq_varid(ncid_ens, trim(sub_varnames(j)), varid_ens) )
            if (xtype == NF90_float) &
               call check( nf90_get_var(ncid_ens, varid_ens, val_4d_f(:, :, :, :), &
                           start=[1, 1, 1, 1], count = [nx, ny, nz, nt]) )
            if (xtype == NF90_double) &
               call check( nf90_get_var(ncid_ens, varid_ens, val_4d_d(:, :, :, :), &
                           start=[1, 1, 1, 1], count = [nx, ny, nz, nt]) )
            intmd_4d_f = intmd_4d_f + val_4d_f
            intmd_4d_d = intmd_4d_d + val_4d_d
         end do
         ! compute the value in log-space
         if (do_logspace) then
            where (intmd_4d_f > 0. .and. abs(intmd_4d_f - missing_value_f) > 1e-6)
               intmd_4d_f = log10(intmd_4d_f)
            else where (abs(intmd_4d_f - missing_value_f) < 1e-6)
               intmd_4d_f = missing_value_f
            else where
               intmd_4d_f = -14
            end where

            where (intmd_4d_d > 0..and. abs(intmd_4d_d - missing_value_d) > 1e-6)
               intmd_4d_d = log10(intmd_4d_d)
            else where (abs(intmd_4d_d - missing_value_d) < 1e-6)
               intmd_4d_d = missing_value_d
            else where
               intmd_4d_d = -14
            end where
         end if
         ! compute standard deviation
         where (abs(intmd_4d_f - missing_value_f) < 1e-6)
            std_4d_f = missing_value_f
         else where
            std_4d_f = std_4d_f + (intmd_4d_f - mean_4d_f)*(intmd_4d_f - mean_4d_f)/Ne
         end where

         where (abs(intmd_4d_d - missing_value_d) < 1e-6)
            std_4d_d = missing_value_d
         else where
            std_4d_d = std_4d_d + (intmd_4d_d - mean_4d_d)*(intmd_4d_d - mean_4d_d)/Ne
         end where
         call check( nf90_close(ncid_ens) )
      end do
      where (std_4d_f > missing_value_f)
         std_4d_f = missing_value_f
      else where
         std_4d_f = sqrt(std_4d_f)
      end where
      where (std_4d_d > missing_value_d)
         std_4d_d = missing_value_d
      else where
         std_4d_d = sqrt(std_4d_d)
      end where
   end subroutine calculate_total_phytoplankton_std

   subroutine redefine_metadata()
      integer :: attnum
      integer :: varid, varid_std
      integer :: ndims
      integer :: xtype
      integer :: nAtts
      integer, allocatable :: dimids(:)

      CHARACTER(LEN=1024) :: varname
      CHARACTER(LEN=1024) :: attname
      CHARACTER(LEN=1024) :: attval
      CHARACTER(LEN=1024) :: pwd

      print *, rank, 'file definition'
      ! get current working directory
      call get_environment_variable('PWD',pwd)
      ! open the mode to redefine metadata
      call check( nf90_redef(ncid) )
      ! write the creating path of the data
      attval = trim(pwd)//'/'//trim(filename)
      call check( nf90_put_att(ncid, NF90_GLOBAL, 'name', trim(attval)) )
      ! write the creating time of the data
      call check( nf90_put_att(ncid, NF90_GLOBAL, 'timeStamp', trim(iso8601())) )
      ! remove uuid of the data  -- some unique strings
      if (nf90_inquire_attribute(ncid, NF90_GLOBAL, 'uuid') == NF90_NOERR) &
         call check( nf90_del_att(ncid, NF90_GLOBAL, 'uuid') )

      ! check if the current file has '_std' defined. 
      ! If so, we don't need to define std variables
      do varid = 1, nVariables + n_new_vars
         ! get the varname of the variable in this file
         call check( nf90_inquire_variable(ncid, varid, name=varname) )
         if (index(trim(varname), '_std') /= 0) then
           call check( nf90_enddef(ncid) )
           return
         end if
      end do

      do varid = 1, nVariables + n_new_vars
         call check( nf90_inquire_variable(ncid, varid, name=varname, xtype=xtype, ndims=ndims, nAtts=nAtts) )

         if (any(varname == exclude_varnames)) cycle

         ! change name for variable attributes
         if (nf90_inquire_attribute(ncid, varid, 'long_name') == NF90_NOERR) then
            call check( nf90_get_att(ncid, varid, 'long_name', attval) )
            attval = 'ensemble mean of '//trim(attval)
            if (do_logspace) attval = trim(attval)//'_in log10 space'
            call check( nf90_put_att(ncid, varid, 'long_name', trim(attval)) )
         end if

         if (nf90_inquire_attribute(ncid, varid, 'standard_name') == NF90_NOERR) then
            call check( nf90_get_att(ncid, varid, 'standard_name', attval) )
            attval = 'ens_mean_'//trim(attval)
            if (do_logspace) attval = trim(attval)//'_in log10 space'
            call check( nf90_put_att(ncid, varid, 'standard_name', trim(attval)) )
         end if   

         ! define standard deviation of the variables
         if (allocated(dimids)) deallocate(dimids)
         allocate(dimids(ndims))
         call check( nf90_inquire_variable(ncid, varid, dimids=dimids) )
         call check( nf90_def_var(ncid, trim(varname)//'_std', xtype, dimids, varid_std) )
         do attnum = 1, nAtts
            call check( nf90_inq_attname(ncid, varid, attnum, attname) )
            call check( nf90_copy_att(ncid, varid, trim(attname), ncid, varid_std) )
            if (trim(attname) == 'long_name') then
               call check( nf90_get_att(ncid, varid, trim(attname), attval) )
               attval = 'ensemble std of '//trim(attval)
               call check( nf90_put_att(ncid, varid_std, trim(attname), trim(attval)) )
            end if
            if (trim(attname) == 'standard_name') then
               call check( nf90_get_att(ncid, varid, trim(attname), attval) )
               attval = 'ens_std_'//trim(attval)
               call check( nf90_put_att(ncid, varid_std, trim(attname), trim(attval) ) )
            end if
         end do
      end do
      call check( nf90_enddef(ncid) )
      if (allocated(dimids)) deallocate(dimids)
   end subroutine redefine_metadata

   subroutine calculate_mean_std()
      integer :: varid, varid_std
      integer :: xtype
      integer :: ndims
      CHARACTER(LEN=1024) :: varname

      real(4), ALLOCATABLE :: mean_4d_f(:, :, :, :)
      real(8), ALLOCATABLE :: mean_4d_d(:, :, :, :)

      real(4), ALLOCATABLE :: val_4d_f(:, :, :, :)
      real(8), ALLOCATABLE :: val_4d_d(:, :, :, :)

      real(4), ALLOCATABLE :: intmd_4d_f(:, :, :, :)
      real(8), ALLOCATABLE :: intmd_4d_d(:, :, :, :)

      real(4), ALLOCATABLE :: std_4d_f(:, :, :, :)
      real(8), ALLOCATABLE :: std_4d_d(:, :, :, :)

      integer :: nc_ierr
      real(8) :: missing_value_d
      real(4) :: missing_value_f
      real(8) :: filled_value_d
      real(4) :: filled_value_f

      if (.not. allocated(val_4d_f)) allocate(val_4d_f(nx, ny, nz, nt))
      if (.not. allocated(val_4d_d)) allocate(val_4d_d(nx, ny, nz, nt))
      if (.not. allocated(intmd_4d_f)) allocate(intmd_4d_f(nx, ny, nz, nt))
      if (.not. allocated(intmd_4d_d)) allocate(intmd_4d_d(nx, ny, nz, nt))
      if (.not. allocated(mean_4d_f)) allocate(mean_4d_f(nx, ny, nz, nt))
      if (.not. allocated(mean_4d_d)) allocate(mean_4d_d(nx, ny, nz, nt))
      if (.not. allocated(std_4d_f)) allocate(std_4d_f(nx, ny, nz, nt))
      if (.not. allocated(std_4d_d)) allocate(std_4d_d(nx, ny, nz, nt))

      do varid = 1, nVariables + n_new_vars
         ! inquire the variable name, datatype and dimensions
         call check( nf90_inquire_variable(ncid, varid, name=varname, xtype=xtype, ndims=ndims) )
         ! no need to calculate mean values for excluded variables
         if (any(trim(varname) == exclude_varnames)) cycle
         ! std variables are calculated together with the mean variable
         if (index(trim(varname), '_std') /= 0) cycle
         ! get the missing value of the variable
         missing_value_f = 0.0
         missing_value_d = 0.0
         if (xtype == NF90_float) then
            nc_ierr = nf90_get_att(ncid, varid, 'missing_value', missing_value_f)
            nc_ierr = nf90_get_att(ncid, varid, '_FillValue', filled_value_f)
         else if (xtype == NF90_double) then
            nc_ierr = nf90_get_att(ncid, varid, 'missing_value', missing_value_d)
            nc_ierr = nf90_get_att(ncid, varid, '_FillValue', filled_value_d)
         end if

         print *, rank, trim(varname),  'writing mean value'
         ! calculate mean value of the ensemble
         call calculate_mean(ndims, xtype, trim(varname), missing_value_f, missing_value_d, &
                             val_4d_f, val_4d_d, intmd_4d_f, intmd_4d_d, mean_4d_f, mean_4d_d)
         ! write the mean value of the ensemble to netCDF file
         if (ndims == 3) then
            if (xtype == NF90_float) &
               call check( nf90_put_var(ncid, varid, mean_4d_f(:, :, 1, :), &
                           start=[1, 1, 1], count = [nx, ny, nt]) )
            if (xtype == NF90_double) &
               call check( nf90_put_var(ncid, varid, mean_4d_d(:, :, 1, :), &
                           start=[1, 1, 1], count = [nx, ny, nt]) )
         else if (ndims == 4) then
            if (xtype == NF90_float) &
               call check( nf90_put_var(ncid, varid, mean_4d_f(:, :, :, :), &
                           start=[1, 1, 1, 1], count = [nx, ny, nz, nt]) )
            if (xtype == NF90_double) &
               call check( nf90_put_var(ncid, varid, mean_4d_d(:, :, :, :), &
                           start=[1, 1, 1, 1], count = [nx, ny, nz, nt]) )
         end if

         ! calculate standard deviation (spread) of the ensemble
         call calculate_std(ndims, xtype, varname, missing_value_f, missing_value_d, &
                            val_4d_f, val_4d_d, intmd_4d_f, intmd_4d_d, mean_4d_f, mean_4d_d, std_4d_f, std_4d_d)
         print *, rank, 'writing std value'
         ! write the standard deviation to the netCDF file
         call check( nf90_inq_varid(ncid, trim(varname)//'_std', varid_std) )
         if (ndims == 3) then
            if (xtype == NF90_float) &
               call check( nf90_put_var(ncid, varid_std, std_4d_f(:, :, 1, :), &
                           start=[1, 1, 1], count = [nx, ny, nt]) )
            if (xtype == NF90_double) &
               call check( nf90_put_var(ncid, varid_std, std_4d_d(:, :, 1, :), &
                           start=[1, 1, 1], count = [nx, ny, nt]) )
         else if (ndims == 4) then
            if (xtype == NF90_float) &
               call check( nf90_put_var(ncid, varid_std, std_4d_f(:, :, :, :), &
                           start=[1, 1, 1, 1], count = [nx, ny, nz, nt]) )
            if (xtype == NF90_double) &
               call check( nf90_put_var(ncid, varid_std, std_4d_d(:, :, :, :), &
                           start=[1, 1, 1, 1], count = [nx, ny, nz, nt]) )
         end if
      end do

      if (allocated(val_4d_f)) deallocate(val_4d_f)
      if (allocated(val_4d_d)) deallocate(val_4d_d)
      if (allocated(intmd_4d_f)) deallocate(intmd_4d_f)
      if (allocated(intmd_4d_d)) deallocate(intmd_4d_d)
      if (allocated(mean_4d_f)) deallocate(mean_4d_f)
      if (allocated(mean_4d_d)) deallocate(mean_4d_d)
      if (allocated(std_4d_f)) deallocate(std_4d_f)
      if (allocated(std_4d_d)) deallocate(std_4d_d)
   end subroutine calculate_mean_std

   subroutine calculate_mean(ndims, xtype, varname, &
                             missing_value_f, missing_value_d, &
                             val_4d_f, val_4d_d, intmd_4d_f, intmd_4d_d, &
                             mean_4d_f, mean_4d_d)
      CHARACTER(LEN=*), intent(in) :: varname
      integer,          intent(in) :: ndims
      integer,          intent(in) :: xtype

      real(4), intent(inout) :: val_4d_f(:, :, :, :)
      real(8), intent(inout) :: val_4d_d(:, :, :, :)
      real(4), intent(inout) :: intmd_4d_f(:, :, :, :)
      real(8), intent(inout) :: intmd_4d_d(:, :, :, :)
      real(8), intent(in)    :: missing_value_d
      real(4), intent(in)    :: missing_value_f
      real(4), intent(out)   :: mean_4d_f(:, :, :, :)
      real(8), intent(out)   :: mean_4d_d(:, :, :, :)

      integer :: ncid_ens
      integer :: varid_ens
      integer :: i

      ! 
      if (trim(varname) == 'chlorophyll' .or. trim(varname) == 'nitrogen') then
         call calculate_total_phytoplankton_mean(trim(varname), xtype, &
                                              missing_value_f, missing_value_d, &
                                              val_4d_f, val_4d_d, intmd_4d_f, intmd_4d_d, &
                                              mean_4d_f, mean_4d_d)
      else
         mean_4d_f = 0.
         mean_4d_d = 0.
         ! 
         do i = 1, Ne
            ! open ensemble member file
            call check( nf90_open(trim(BaseDir)//'/ensemble_'//trim(str(i))//'/'//trim(filename)//'_'//i_file_index//'.nc', &
                                  nf90_nowrite, ncid_ens) )
            ! get the variable id
            call check( nf90_inq_varid(ncid_ens, trim(varname), varid_ens) )
            ! read variable values
            val_4d_f = 0.
            val_4d_d = 0.
            if (ndims == 3) then
               if (xtype == NF90_float) &
                  call check( nf90_get_var(ncid_ens, varid_ens, val_4d_f(:, :, 1, :), &
                              start=[1, 1, 1], count = [nx, ny, nt]) )
               if (xtype == NF90_double) &
                  call check( nf90_get_var(ncid_ens, varid_ens, val_4d_d(:, :, 1, :), &
                              start=[1, 1, 1], count = [nx, ny, nt]) )
            else if (ndims == 4) then
               if (xtype == NF90_float) &
                  call check( nf90_get_var(ncid_ens, varid_ens, val_4d_f(:, :, :, :), &
                              start=[1, 1, 1, 1], count = [nx, ny, nz, nt]) )
               if (xtype == NF90_double) &
                  call check( nf90_get_var(ncid_ens, varid_ens, val_4d_d(:, :, :, :), &
                              start=[1, 1, 1, 1], count = [nx, ny, nz, nt]) )
            end if
            ! compute the average value
            if (do_logspace) then
               where (val_4d_f > 0. .and. abs(val_4d_f - missing_value_f) > 1e-6)
                  val_4d_f = log10(val_4d_f)
               else where (abs(val_4d_f - missing_value_f) < 1e-6)
                  val_4d_f = missing_value_f
               else where
                  val_4d_f = -14
               end where
               where (val_4d_d > 0..and. abs(val_4d_d - missing_value_d) > 1e-6)
                  val_4d_d = log10(val_4d_d)
               else where (abs(val_4d_d - missing_value_d) < 1e-6)
                  val_4d_d = missing_value_d
               else where
                  val_4d_d = -14
               end where
            end if
            ! calculate mean value
            where (abs(val_4d_f - missing_value_f) < 1e-6)
               mean_4d_f = missing_value_f
            else where
               mean_4d_f = mean_4d_f + val_4d_f/Ne
            end where

            where (abs(val_4d_d - missing_value_d) < 1e-6)
               mean_4d_d = missing_value_d
            else where
               mean_4d_d = mean_4d_d + val_4d_d/Ne
            end where
            call check( nf90_close(ncid_ens) )
         end do

         where (mean_4d_f > missing_value_f)
            mean_4d_f = missing_value_f
         end where
         where (mean_4d_d > missing_value_d)
            mean_4d_d = missing_value_d
         end where
      end if
   end subroutine calculate_mean

   subroutine calculate_std(ndims, xtype, varname, missing_value_f, missing_value_d, &
                            val_4d_f, val_4d_d, intmd_4d_f, intmd_4d_d, &
                            mean_4d_f, mean_4d_d, std_4d_f, std_4d_d)
      CHARACTER(LEN=*), intent(in) :: varname
      integer,          intent(in) :: ndims
      integer,          intent(in) :: xtype

      real(4), intent(in) :: mean_4d_f(:, :, :, :)
      real(8), intent(in) :: mean_4d_d(:, :, :, :)
      real(8), intent(in) :: missing_value_d
      real(4), intent(in) :: missing_value_f

      real(4), intent(inout) :: val_4d_f(:, :, :, :)
      real(8), intent(inout) :: val_4d_d(:, :, :, :)
      real(4), intent(inout) :: intmd_4d_f(:, :, :, :)
      real(8), intent(inout) :: intmd_4d_d(:, :, :, :)
      real(4), intent(out) :: std_4d_f(:, :, :, :)
      real(8), intent(out) :: std_4d_d(:, :, :, :)

      integer :: ncid_ens
      integer :: varid_ens
      integer :: i

      if (trim(varname) == 'chlorophyll' .or. trim(varname) == 'nitrogen') then
         call calculate_total_phytoplankton_std(varname, xtype, &
                                                 missing_value_f, missing_value_d, &
                                                 mean_4d_f, mean_4d_d, &
                                                 val_4d_f, val_4d_d, intmd_4d_f, intmd_4d_d, &
                                                 std_4d_f, std_4d_d)
      else
         ! get ensemble std
         std_4d_f = 0.
         std_4d_d = 0.

         do i = 1, Ne
            call check( nf90_open(trim(BaseDir)//'/ensemble_'//trim(str(i))//'/'//trim(filename)//'_'//i_file_index//'.nc', &
                                  nf90_nowrite, ncid_ens) )
            call check( nf90_inq_varid(ncid_ens, trim(varname), varid_ens) )
            val_4d_f = 0.
            val_4d_d = 0.

            if (ndims == 3) then
               if (xtype == NF90_float) &
                  call check( nf90_get_var(ncid_ens, varid_ens, val_4d_f(:, :, 1, :), &
                              start=[1, 1, 1], count = [nx, ny, nt]) )
               if (xtype == NF90_double) &
                  call check( nf90_get_var(ncid_ens, varid_ens, val_4d_d(:, :, 1, :), &
                              start=[1, 1, 1], count = [nx, ny, nt]) )
            else if (ndims == 4) then
               if (xtype == NF90_float) &
                  call check( nf90_get_var(ncid_ens, varid_ens, val_4d_f(:, :, :, :), &
                              start=[1, 1, 1, 1], count = [nx, ny, nz, nt]) )
               if (xtype == NF90_double) &
                  call check( nf90_get_var(ncid_ens, varid_ens, val_4d_d(:, :, :, :), &
                              start=[1, 1, 1, 1], count = [nx, ny, nz, nt]) )
            end if
            if (do_logspace) then
               where (val_4d_f > 0. .and. abs(val_4d_f - missing_value_f) > 1e-6)
                  val_4d_f = log10(val_4d_f)
               else where (abs(val_4d_f - missing_value_f) < 1e-6)
                  val_4d_f = missing_value_f
               else where
                  val_4d_f = -14
               end where
               where (val_4d_d > 0..and. abs(val_4d_d - missing_value_d) > 1e-6)
                  val_4d_d = log10(val_4d_d)
               else where (abs(val_4d_d - missing_value_d) < 1e-6)
                  val_4d_d = missing_value_d
               else where
                  val_4d_d = -14
               end where
            end if
            ! calculate mean value
            where (abs(val_4d_f - missing_value_f) < 1e-6)
               std_4d_f = missing_value_f
            else where
               std_4d_f = std_4d_f + (val_4d_f - mean_4d_f)*(val_4d_f - mean_4d_f)/Ne
            end where

            where (abs(val_4d_d - missing_value_d) < 1e-6)
               std_4d_d = missing_value_d
            else where
               std_4d_d = std_4d_d + (val_4d_d - mean_4d_d)*(val_4d_d - mean_4d_d)/Ne
            end where
            call check( nf90_close(ncid_ens) )
         end do
         where (std_4d_f > missing_value_f)
            std_4d_f = missing_value_f
         else where
            std_4d_f = sqrt(std_4d_f)
         end where
         where (std_4d_d > missing_value_d)
            std_4d_d = missing_value_d
         else where
            std_4d_d = sqrt(std_4d_d)
         end where
      end if
   end subroutine calculate_std

   character(len=4) function str(k)

      implicit none

      integer, intent(in) :: k   !< number

      CHARACTER(LEN=1) :: kc
      write (kc, '(I1)') FLOOR(LOG10(REAL(k))+1)
      write (str, '(i'//kc//')') k
   end function str
end program merge
