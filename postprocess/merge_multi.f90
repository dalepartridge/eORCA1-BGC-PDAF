program merge
use netcdf
use mpi
implicit none
integer, parameter :: Ne = 30
! integer :: nx = 362, ny = 332, nz = 75, Ne = 30
integer :: nx, ny, nz = 1
integer :: MPIerr
logical :: iniflag
integer :: rank
integer :: npes
integer :: chunksize, remainder, displs
integer :: nfiles, i_file
CHARACTER(len=4) :: nfiles_c, i_file_index
CHARACTER(len=1) :: doLogC
CHARACTER(len=1024) :: BaseDir
CHARACTER(len=1024) :: filename
logical             :: do_logspace

integer :: ncid
integer :: ierr
integer :: dimid, nt
integer :: nVariables, nAttributes
CHARACTER(LEN=20), dimension(8) :: exclude_varnames


exclude_varnames = [CHARACTER(LEN=20) :: 'nav_lat', 'nav_lon', 'deptht', 'deptht_bounds', &
                     'time_centered', 'time_centered_bounds', &
                     'time_counter', 'time_counter_bounds']

CALL MPI_Initialized(iniflag, MPIerr)
IF (.not.iniflag) CALL MPI_Init(MPIerr)
CALL MPI_Comm_size(MPI_COMM_WORLD, npes, MPIerr)
CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, MPIerr)
CALL get_command_argument(1, BaseDir)
CALL get_command_argument(2, filename)
CALL get_command_argument(3, nfiles_c)
CALL get_command_argument(4, doLogC)
do_logspace = .false.
if (doLogC == 'T') do_logspace = .true.

! get nfiles
read(nfiles_c, *) nfiles

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
   if (do_logspace) then
      print *, 'Opening ',  trim(BaseDir)//'/'//trim(filename)//'-log_'//i_file_index//'.nc'
      call check( nf90_open(trim(BaseDir)//'/'//trim(filename)//'-log_'//i_file_index//'.nc', nf90_write, ncid) )
   else
      print *, 'Opening ',  trim(BaseDir)//'/'//trim(filename)//'_'//i_file_index//'.nc'
      call check( nf90_open(trim(BaseDir)//'/'//trim(filename)//'_'//i_file_index//'.nc', nf90_write, ncid) )
   endif
   call check( nf90_inquire(ncid, nVariables=nVariables, nAttributes=nAttributes) )

   call check( nf90_inq_dimid(ncid, 'time_counter', dimid) )
   call check( nf90_inquire_dimension(ncid, dimid, len=nt) )
   ! if (npes /= nt) then
   !    if (rank == 0) &
   !       print *, npes, 'number of PEs with ', nt, 'steps'
   !    CALL MPI_Abort(MPI_COMM_WORLD, 1, MPIerr)
   ! endif

   call check( nf90_inq_dimid(ncid, 'x', dimid) )
   call check( nf90_inquire_dimension(ncid, dimid, len=nx) )
   call check( nf90_inq_dimid(ncid, 'y', dimid) )
   call check( nf90_inquire_dimension(ncid, dimid, len=ny) )
   ierr = nf90_inq_dimid(ncid, 'deptht', dimid)
   if (ierr == nf90_noerr) &
      call check( nf90_inquire_dimension(ncid, dimid, len=nz) )

   ! add ensemble standard deviation and mean information
   ! delete uuid and add new time stamp
   call redefine_metadata()
   ! put calculated ensemble standard deviation and mean to the netCDF file
   call calculate_mean_std()

   print *, rank, 'Closing ', trim(BaseDir)//'/'//trim(filename)//'_'//i_file_index//'.nc'
   call check( nf90_close(ncid) )
end do
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
      call get_environment_variable('PWD',pwd)
      call check( nf90_redef(ncid) )
      attval = trim(pwd)//'/'//trim(filename)
      call check( nf90_put_att(ncid, NF90_GLOBAL, 'name', trim(attval)) )
      call check( nf90_put_att(ncid, NF90_GLOBAL, 'timeStamp', trim(iso8601())) )
      if (nf90_inquire_attribute(ncid, NF90_GLOBAL, 'uuid') == NF90_NOERR) &
         call check( nf90_del_att(ncid, NF90_GLOBAL, 'uuid') )

      do varid = 1, nVariables
         call check( nf90_inquire_variable(ncid, varid, name=varname, xtype=xtype, ndims=ndims, nAtts=nAtts) )
         if (index(trim(varname), '_std') /= 0) then
           call check( nf90_enddef(ncid) )
           return
         end if
      end do

      do varid = 1, nVariables
         call check( nf90_inquire_variable(ncid, varid, name=varname, xtype=xtype, ndims=ndims, nAtts=nAtts) )

         if (any(varname == exclude_varnames)) cycle
         print *, rank, varid, nVariables, trim(varname)
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
               if (do_logspace) attval = trim(attval)//'_in log10 space'
               call check( nf90_put_att(ncid, varid_std, trim(attname), trim(attval)) )
            end if
            if (trim(attname) == 'standard_name') then
               call check( nf90_get_att(ncid, varid, trim(attname), attval) )
               attval = 'ens_std_'//trim(attval)
               if (do_logspace) attval = trim(attval)//'_in log10 space'
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

      real(4), ALLOCATABLE :: std_4d_f(:, :, :, :)
      real(8), ALLOCATABLE :: std_4d_d(:, :, :, :)

      integer :: nc_ierr
      real(8) :: missing_value_d
      real(4) :: missing_value_f

      if (.not. allocated(val_4d_f)) allocate(val_4d_f(nx, ny, nz, nt))
      if (.not. allocated(val_4d_d)) allocate(val_4d_d(nx, ny, nz, nt))
      if (.not. allocated(mean_4d_f)) allocate(mean_4d_f(nx, ny, nz, nt))
      if (.not. allocated(mean_4d_d)) allocate(mean_4d_d(nx, ny, nz, nt))
      if (.not. allocated(std_4d_f)) allocate(std_4d_f(nx, ny, nz, nt))
      if (.not. allocated(std_4d_d)) allocate(std_4d_d(nx, ny, nz, nt))

      do varid = 1, nVariables
         call check( nf90_inquire_variable(ncid, varid, name=varname, xtype=xtype, ndims=ndims) )
         if (any(trim(varname) == exclude_varnames)) cycle
         if (index(trim(varname), '_std') /= 0) cycle
 
         missing_value_f = 0.0
         missing_value_d = 0.
         if (xtype == NF90_float) then
            nc_ierr = nf90_get_att(ncid, varid, 'missing_value', missing_value_f)
         else if (xtype == NF90_double) then
            nc_ierr = nf90_get_att(ncid, varid, 'missing_value', missing_value_d)
         end if

         print *, rank, trim(varname),  'writing mean value'
         call calculate_mean(ndims, xtype, trim(varname), missing_value_f, missing_value_d, &
                             val_4d_f, val_4d_d, mean_4d_f, mean_4d_d)
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

         call calculate_std(ndims, xtype, varname, missing_value_f, missing_value_d, &
                            val_4d_f, val_4d_d, mean_4d_f, mean_4d_d, std_4d_f, std_4d_d)
         print *, rank, 'writing std value'
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
      if (allocated(mean_4d_f)) deallocate(mean_4d_f)
      if (allocated(mean_4d_d)) deallocate(mean_4d_d)
      if (allocated(std_4d_f)) deallocate(std_4d_f)
      if (allocated(std_4d_d)) deallocate(std_4d_d)
   end subroutine calculate_mean_std

   subroutine calculate_mean(ndims, xtype, varname, missing_value_f, missing_value_d,&
                             val_4d_f, val_4d_d, mean_4d_f, mean_4d_d)
      CHARACTER(LEN=*), intent(in) :: varname
      integer,          intent(in) :: ndims
      integer,          intent(in) :: xtype

      real(4), intent(inout) :: val_4d_f(:, :, :, :)
      real(8), intent(inout) :: val_4d_d(:, :, :, :)
      real(8), intent(in) :: missing_value_d
      real(4), intent(in) :: missing_value_f
      real(4), intent(out) :: mean_4d_f(:, :, :, :)
      real(8), intent(out) :: mean_4d_d(:, :, :, :)

      integer :: ncid_ens
      integer :: varid_ens
      integer :: i


      mean_4d_f = 0.
      mean_4d_d = 0.

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
            where (val_4d_f > 0.)
               val_4d_f = log10(val_4d_f)
            else where
               val_4d_f = -14
            end where
            where (val_4d_d > 0.)
               val_4d_d = log10(val_4d_d)
            else where
               val_4d_d = -14
            end where
         end if
         mean_4d_f = mean_4d_f + val_4d_f/Ne
         mean_4d_d = mean_4d_d + val_4d_d/Ne
         call check( nf90_close(ncid_ens) )
      end do
      where (abs(val_4d_f - missing_value_f) < 1e-6)
         mean_4d_f = missing_value_f
      end where

      where (abs(val_4d_d - missing_value_d) < 1e-6)
         mean_4d_d = missing_value_d
      end where
   end subroutine calculate_mean

   subroutine calculate_std(ndims, xtype, varname, missing_value_f, missing_value_d, &
                            val_4d_f, val_4d_d, mean_4d_f, mean_4d_d, std_4d_f, std_4d_d)
      CHARACTER(LEN=*), intent(in) :: varname
      integer,          intent(in) :: ndims
      integer,          intent(in) :: xtype

      real(4), intent(in) :: mean_4d_f(:, :, :, :)
      real(8), intent(in) :: mean_4d_d(:, :, :, :)
      real(8), intent(in) :: missing_value_d
      real(4), intent(in) :: missing_value_f

      real(4), intent(inout) :: val_4d_f(:, :, :, :)
      real(8), intent(inout) :: val_4d_d(:, :, :, :)

      real(4), intent(out) :: std_4d_f(:, :, :, :)
      real(8), intent(out) :: std_4d_d(:, :, :, :)

      integer :: ncid_ens
      integer :: varid_ens
      integer :: i


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
            where (val_4d_f > 0.)
               val_4d_f = log10(val_4d_f)
            else where
               val_4d_f = -14
            end where
            where (val_4d_d > 0.)
               val_4d_d = log10(val_4d_d)
            else where
               val_4d_d = -14
            end where
         end if
         std_4d_f = std_4d_f + (val_4d_f - mean_4d_f)*(val_4d_f - mean_4d_f)/Ne
         std_4d_d = std_4d_d + (val_4d_d - mean_4d_d)*(val_4d_d - mean_4d_d)/Ne
         call check( nf90_close(ncid_ens) )
      end do
      std_4d_f = sqrt(std_4d_f)
      std_4d_d = sqrt(std_4d_d)
      where (abs(val_4d_f - missing_value_f) < 1e-6)
         std_4d_f = missing_value_f
      end where

      where (abs(val_4d_d - missing_value_d) < 1e-6)
         std_4d_d = missing_value_d
      end where
   end subroutine calculate_std

   character(len=4) function str(k)

      implicit none

      integer, intent(in) :: k   !< number

      CHARACTER(LEN=1) :: kc
      write (kc, '(I1)') FLOOR(LOG10(REAL(k))+1)
      write (str, '(i'//kc//')') k
   end function str
end program merge
