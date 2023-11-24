module mod_obs
implicit none
contains
   subroutine getObs(obstype, obsfilename, indexFilename, &
                     ind_i, ind_j, ngridpoints, model_obs, model_obs_unc)
      character(*), intent(in) :: obstype
      character(*), intent(in) :: obsfilename
      character(*), intent(in) :: indexFilename
      integer, intent(in) :: ngridpoints
      integer, intent(in) :: ind_i(ngridpoints), ind_j(ngridpoints)
      real(kind=8), intent(inout) :: model_obs(ngridpoints)
      real(kind=8), intent(inout) :: model_obs_unc(ngridpoints)

      integer :: dim_olat, dim_olon
      integer(kind=8), allocatable :: mobsi(:), mobsj(:)
      real(kind=8), allocatable :: phyto_aux(:, :)

      print *, 'Observation filename: ', obsfilename
      print *, 'obstype: ', obstype
      print *, 'indexFilename: ', indexFilename
      print *, 'ngridpoints', ngridpoints

      call getObsDim(indexFilename, dim_olat, dim_olon)
      allocate(mobsi(dim_olon*dim_olat), mobsj(dim_olon*dim_olat))
      allocate(phyto_aux(dim_olon*dim_olat, 4))
      call getModelIndexOfObs(indexFilename, mobsi, mobsj)
      print *, 'getting observation data'
      call getModelVarFromObs(obsfilename, obstype, phyto_aux)
      print *, 'processing data'
      call getObsOnModel(mobsi, mobsj, ind_i, ind_j, ngridpoints, &
                        phyto_aux, model_obs, model_obs_unc)
   end subroutine getObs

   !< get the rank of observation in the model foreacst
   !< ensembles
   subroutine getRank(nlength, Ne, nsize, x, rank)
      integer, intent(in) :: nlength, Ne, nsize
      logical, intent(in) :: x(nsize)
      integer, intent(inout) :: rank(nlength)
      integer :: i, j
      do i = 1, nlength
         do j = 1, Ne+1
            if (x(j + (i-1)*(Ne+1))) exit
         end do
         rank(i) = j
      end do
   end subroutine getRank

   !< gather observations for each model grid point
   !< The median of all observations are chosen
   subroutine getObsOnModel(mobsi, mobsj, ind_i, ind_j, ngridpoints, &
                            phyto_aux, model_obs, model_obs_unc)
      integer(kind=8), intent(in) :: mobsi(:), mobsj(:)
      integer, intent(in) :: ind_i(:), ind_j(:)
      integer, intent(in) :: ngridpoints
      real(kind=8), intent(in) :: phyto_aux(:, :)
      real(kind=8), intent(inout) :: model_obs(:), model_obs_unc(:)

      integer :: nobs
      integer :: i, j, k
      integer :: cntobs
      real(kind=8), allocatable :: obslist(:), obsunc(:)

      nobs = size(mobsi, dim=1)
      allocate(obslist(nobs), obsunc(nobs))

      model_obs = 2e20
      do k = 1, ngridpoints
         i = ind_i(k) + 1
         j = ind_j(k) + 1

         if (mod(k-1, ngridpoints/10) == 0) &
            print *, 'computed at: ', (k - 1)*100./ngridpoints,'%'

         cntobs = count ( (mobsi == i) .and. (mobsj == j) .and. &
                           (phyto_aux(:, 4) < 1e14) &
                        )
         if (cntobs > 0) then
            obslist(:cntobs) = pack(phyto_aux(:, 4), &
                                       (mobsi == i) .and. &
                                       (mobsj == j) .and. &
                                       (phyto_aux(:, 4) < 1e14))
            obsunc(:cntobs) = pack(phyto_aux(:, 3), &
                                       (mobsi == i) .and. &
                                       (mobsj == j) .and. &
                                       (phyto_aux(:, 4) < 1e14))
            call quicksort(cntobs, obslist(:cntobs), obsunc(:cntobs))
            if ( mod(cntobs, 2) == 0 ) then
               model_obs(k) =(obslist(cntobs/2+1) + obslist(cntobs/2))*0.5
               model_obs_unc(k) =(obsunc(cntobs/2+1) + obsunc(cntobs/2))*0.5
            else
               model_obs(k) = obslist(cntobs/2+1)
               model_obs_unc(k) = obsunc(cntobs/2+1)
            end if

         endif
      enddo
      deallocate(obslist, obsunc)
   end subroutine getObsOnModel

   !< get observation dimensions (number of lat/lon)
   subroutine getObsDim(indexFilename, dim_olat, dim_olon)
      use netcdf
      character(*), intent(in) :: indexFilename
      integer, intent(out) :: dim_olat, dim_olon

      integer :: ncid
      integer :: dimid

      call check( nf90_open(indexFilename, nf90_nowrite, ncid) )
      call check( nf90_inq_dimid(ncid, "ilat", dimid) )
      call check( nf90_Inquire_dimension(ncid, dimid, len=dim_olat) )
      call check( nf90_inq_dimid(ncid, "ilon", dimid) )
      call check( nf90_Inquire_dimension(ncid, dimid, len=dim_olon) )
      call check( nf90_close(ncid) )
   end subroutine getObsDim

   !< The corresponding grid point index of the observation.
   !< This is calculated when one observation is located within
   !< a box surrounding the model grid point. The length of the
   !< box is determined by the grid point interval.
   subroutine getModelIndexOfObs(indexFilename, mobsi, mobsj)
      use netcdf
      character(*), intent(in) :: indexFilename
      integer(kind=8), intent(inout) :: mobsi(:), mobsj(:)

      integer :: ncid
      integer :: dimid
      integer :: varid
      integer :: dim_olat, dim_olon

      call check( nf90_open(indexFilename, nf90_nowrite, ncid) )
      call check( nf90_inq_dimid(ncid, "ilat", dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=dim_olat) )
      call check( nf90_inq_dimid(ncid, "ilon", dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=dim_olon) )
      call check( nf90_inq_varid(ncid, "mobsi", varid) )
      call check( nf90_get_var(ncid, varid, mobsi, &
                  start=[1, 1], count=[dim_olon, dim_olat]) )
      call check( nf90_inq_varid(ncid, "mobsj", varid) )
      call check( nf90_get_var(ncid, varid, mobsj, &
                  start=[1, 1], count=[dim_olon, dim_olat]) )
      call check( nf90_close(ncid) )
      ! associated the high resolution grid index to the actual model grid index
      mobsi = mobsi + (mod(mobsi, 2) - 1)
      mobsi = (mobsi + 1)/2
      mobsj = mobsj + (mod(mobsj, 2) - 1)
      mobsj = (mobsj + 1)/2
   end subroutine getModelIndexOfObs

   !< convert from observation carbon data to non-/diatom phytoplankton
   subroutine getModelVarFromObs(obsfilename, obstype, phyto_aux)
      use netcdf
      character(*), intent(in) :: obsfilename
      character(*), intent(in) :: obstype
      real(kind=8), intent(inout) :: phyto_aux(:, :)
      integer :: ncid
      integer :: dimid
      integer :: varid
      integer :: dim_olat, dim_olon

      call check( nf90_open(obsfilename, nf90_nowrite, ncid) )
      if (obstype == 'chlo') then
         print *, 'reading chlorophyll data'
         call check( nf90_inq_dimid(ncid, "lat", dimid) )
         call check( nf90_Inquire_dimension(ncid, dimid, len=dim_olat) )
         call check( nf90_inq_dimid(ncid, "lon", dimid) )
         call check( nf90_Inquire_dimension(ncid, dimid, len=dim_olon) )
         call check( nf90_inq_varid(ncid, 'chlor_a', varid) )
         call check( nf90_get_var(ncid, varid, phyto_aux(:, 1), &
                                  start=[1,1,1], &
                                  count=[dim_olon,dim_olat,1]) )
         call check( nf90_inq_varid(ncid, 'chlor_a_log10_bias', varid) )
         call check( nf90_get_var(ncid, varid, phyto_aux(:, 2), &
                                  start=[1,1,1], &
                                  count=[dim_olon,dim_olat,1]) )

         call check( nf90_inq_varid(ncid, 'chlor_a_log10_rmsd', varid) )
         call check( nf90_get_var(ncid, varid, phyto_aux(:, 3), &
                                 start=[1,1,1], &
                                 count=[dim_olon,dim_olat,1]) )
         call check( nf90_close(ncid) )

         ! bias correction for chlor_a:
         phyto_aux(:, 3) = sqrt(abs(phyto_aux(:, 3)*phyto_aux(:, 3) - phyto_aux(:, 2)*phyto_aux(:, 2)))
         phyto_aux(:, 4) = log10(phyto_aux(:, 1)) + phyto_aux(:, 2) - 0.5*log(10.)*phyto_aux(:, 3)
         !phyto_aux(:, 2) = 1.057*(1-exp(-0.851*phyto_aux(:, 4)))
         !phyto_aux(:, 1) = phyto_aux(:, 4) - phyto_aux(:, 2)
      else if (obstype == 'pc') then
         call check( nf90_inq_dimid(ncid, "lat", dimid) )
         call check( nf90_Inquire_dimension(ncid, dimid, len=dim_olat) )
         call check( nf90_inq_dimid(ncid, "lon", dimid) )
         call check( nf90_Inquire_dimension(ncid, dimid, len=dim_olon) )
         call check( nf90_inq_varid(ncid, 'C_microphyto', varid) )
         call check( nf90_get_var(ncid, varid, phyto_aux(:, 1), &
                                  start=[1,1,1], &
                                  count=[dim_olon,dim_olat, 1]) )

         call check( nf90_inq_varid(ncid, 'C_nanophyto', varid) )
         call check( nf90_get_var(ncid, varid, phyto_aux(:, 2), &
                                  start=[1,1,1], &
                                  count=[dim_olon,dim_olat, 1]) )

         call check( nf90_inq_varid(ncid, 'C_picophyto', varid) )
         call check( nf90_get_var(ncid, varid, phyto_aux(:, 3), &
                                  start=[1,1,1], &
                                  count=[dim_olon,dim_olat, 1]) )
         call check( nf90_close(ncid) )

         ! bias correction for chlor_a:
         phyto_aux(:, 4) = phyto_aux(:, 1) + phyto_aux(:, 2) + phyto_aux(:, 3)
         phyto_aux(:, 4) = phyto_aux(:, 4)/6.625
         phyto_aux(:, 3) = 0.3*phyto_aux(:, 4)
      else
         print *, 'unrecoginised obstype:', obstype
      end if

   end subroutine getModelVarFromObs

   !< A quick sort subroutine
   recursive subroutine quicksort(n, obs, obsunc)
      integer  , intent(in)    :: n
      real(kind=8), intent(inout) :: obs(n)
      real(kind=8), intent(inout) :: obsunc(n)

      ! local variables
      integer      :: left, right
      real(kind=8)    :: pivot
      real(kind=8)    :: temp

      pivot = obs((1 + n)/2)   ! choice a random pivot (not best performance, but avoids worst-case)
      left = 1
      right = n
      ! partition loop
      do
         do while (obs(right) > pivot)
            right = right - 1
         end do
         do while (obs(left) < pivot)
            left = left + 1
         end do
         if (left >= right) exit

         temp = obs(left)
         obs(left) = obs(right)
         obs(right) = temp

         temp = obsunc(left)
         obsunc(left) = obsunc(right)
         obsunc(right) = temp

         left = left + 1
         right = right - 1
      end do

      if (left - 1 > 1) call quicksort(left-1, obs(:left-1), obsunc(:left-1))
      if (right + 1 < n) call quicksort(n - right, obs(right+1:), obsunc(right+1:))
   end subroutine quicksort

   !< check netCDF error
   subroutine check(status)
      use netcdf
      ! *** Aruments ***
      integer, intent ( in) :: status   ! Reading status

      if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         stop
      end if
   end subroutine check
end module mod_obs
