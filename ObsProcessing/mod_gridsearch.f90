module mod_gridsearch
implicit none

! INTEGER, PUBLIC, PARAMETER ::   sp = SELECTED_REAL_KIND( 6, 37)   !: single precision (real 4)
! INTEGER, PUBLIC, PARAMETER ::   dp = SELECTED_REAL_KIND(12,307)   !: double precision (real 8)
! INTEGER, PUBLIC, PARAMETER ::   wp = dp                              !: working precision
! logical :: ln_grid_global = .false.
contains
   SUBROUTINE obs_grd_bruteforce( kpi, kpj, kpiglo, kpjglo,       &
      &                            kldi, klei, kldj, klej,         &
      &                            kmyproc, ktotproc,              &
      &                            pglam, pgphi, pmask,            &
      &                            kobs, plam, pphi, kobsi, kobsj, &
      &                                   kproc)
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE obs_grd_bruteforce ***
      !!
      !! ** Purpose : Search gridpoints to find the grid box containing
      !!              the observations
      !!
      !! ** Method  : Call to linquad
      !!
      !! ** Action  : Return kproc holding the observation and kiobsi,kobsj
      !!              valid on kproc=kmyproc processor only.
      !!   
      !! History :
      !!        !  2001-11  (N. Daget, A. Weaver)
      !!        !  2006-03  (A. Weaver) NEMOVAR migration.
      !!        !  2006-05  (K. Mogensen) Moved to to separate routine.
      !!        !  2007-10  (A. Vidard) Bug fix in wrap around checks; cleanup
      !!----------------------------------------------------------------------

      !! * Arguments
      INTEGER, INTENT(IN) :: kpi                ! Number of local longitudes
      INTEGER, INTENT(IN) :: kpj                ! Number of local latitudes
      INTEGER, INTENT(IN) :: kpiglo             ! Number of global longitudes
      INTEGER, INTENT(IN) :: kpjglo             ! Number of global latitudes
      INTEGER, INTENT(IN) :: kldi               ! Start of inner domain in i
      INTEGER, INTENT(IN) :: klei               ! End of inner domain in i
      INTEGER, INTENT(IN) :: kldj               ! Start of inner domain in j
      INTEGER, INTENT(IN) :: klej               ! End of inner domain in j
      INTEGER, INTENT(IN) :: kmyproc            ! Processor number for MPP
      INTEGER, INTENT(IN) :: ktotproc           ! Total number of processors
      REAL(KIND=8), DIMENSION(kpi,kpj), INTENT(IN) :: &
         & pglam,   &               ! Grid point longitude
         & pgphi,   &               ! Grid point latitude
         & pmask                    ! Grid point mask
      INTEGER,INTENT(IN) :: kobs                ! Size of the observation arrays
      REAL(KIND=8), DIMENSION(kobs), INTENT(IN) :: &
         & plam, &                  ! Longitude of obsrvations 
         & pphi                     ! Latitude of observations
      INTEGER, DIMENSION(kobs), INTENT(OUT) :: &
         & kobsi, &                 ! I-index of observations 
         & kobsj, &                 ! J-index of observations 
         & kproc                    ! Processor number of observations
  
      !! * Local declarations
      REAL(8), DIMENSION(:), ALLOCATABLE :: &
         & zplam, zpphi
      REAL(8) :: zlammax
      REAL(8) :: zlam
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jk
      INTEGER :: jo
      INTEGER :: jlon
      INTEGER :: jlat
      INTEGER :: joffset
      INTEGER :: jostride
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: &
         & zlamg, &
         & zphig, &
         & zmskg, &
         & zphitmax,&
         & zphitmin,&
         & zlamtmax,&
         & zlamtmin
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: &
         & llinvalidcell
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zlamtm,  &
         & zphitm

      !-----------------------------------------------------------------------
      ! Define grid setup for grid search
      !-----------------------------------------------------------------------
      jlon     = kpi
      jlat     = kpj
      joffset  = 0
      jostride = 1
      !-----------------------------------------------------------------------
      ! Set up data for grid search
      !-----------------------------------------------------------------------
      ALLOCATE( &
         & zlamg(jlon,jlat),             &
         & zphig(jlon,jlat),             &
         & zmskg(jlon,jlat),             &
         & zphitmax(jlon-1,jlat-1),      &
         & zphitmin(jlon-1,jlat-1),      &
         & zlamtmax(jlon-1,jlat-1),      &
         & zlamtmin(jlon-1,jlat-1),      &
         & llinvalidcell(jlon-1,jlat-1), &
         & zlamtm(4,jlon-1,jlat-1),      &
         & zphitm(4,jlon-1,jlat-1)       &
         & )
      !-----------------------------------------------------------------------
      ! Copy data to local arrays
      !-----------------------------------------------------------------------

         DO jj = 1, jlat
            DO ji = 1, jlon
               zlamg(ji,jj) = pglam(ji,jj)
               zphig(ji,jj) = pgphi(ji,jj)
               zmskg(ji,jj) = pmask(ji,jj)
            END DO
         END DO

      !-----------------------------------------------------------------------
      ! Copy longitudes and latitudes
      !-----------------------------------------------------------------------
      ALLOCATE( &
         & zplam(kobs), &
         & zpphi(kobs)  &
         & )
      DO jo = 1, kobs
         zplam(jo) = plam(jo)
         zpphi(jo) = pphi(jo)
      END DO
      !-----------------------------------------------------------------------
      ! Set default values for output
      !-----------------------------------------------------------------------
      kproc(:) = -1
      kobsi(:) = -1
      kobsj(:) = -1
      !-----------------------------------------------------------------------
      ! Copy grid positions to temporary arrays and renormalize to 0 to 360.
      !-----------------------------------------------------------------------
      DO jj = 1, jlat-1
         DO ji = 1, jlon-1
            zlamtm(1,ji,jj) = zlamg(ji  ,jj  )
            zphitm(1,ji,jj) = zphig(ji  ,jj  )
            zlamtm(2,ji,jj) = zlamg(ji+1,jj  )
            zphitm(2,ji,jj) = zphig(ji+1,jj  )
            zlamtm(3,ji,jj) = zlamg(ji+1,jj+1)
            zphitm(3,ji,jj) = zphig(ji+1,jj+1)
            zlamtm(4,ji,jj) = zlamg(ji  ,jj+1)
            zphitm(4,ji,jj) = zphig(ji  ,jj+1)
         END DO
      END DO
      WHERE ( zlamtm(:,:,:) < 0.0_8 )
         zlamtm(:,:,:) = zlamtm(:,:,:) + 360.0_8
      END WHERE
      WHERE ( zlamtm(:,:,:) > 360.0_8 )
         zlamtm(:,:,:) = zlamtm(:,:,:) - 360.0_8
      END WHERE
      !-----------------------------------------------------------------------
      ! Handle case of the wraparound; beware, not working with orca180
      !-----------------------------------------------------------------------
      DO jj = 1, jlat-1
         DO ji = 1, jlon-1
            zlammax = MAXVAL( zlamtm(:,ji,jj) )
            WHERE (zlammax - zlamtm(:, ji, jj) > 180 ) &
               & zlamtm(:,ji,jj) = zlamtm(:,ji,jj) + 360._8
            zphitmax(ji,jj) = MAXVAL(zphitm(:,ji,jj))
            zphitmin(ji,jj) = MINVAL(zphitm(:,ji,jj))
            zlamtmax(ji,jj) = MAXVAL(zlamtm(:,ji,jj))
            zlamtmin(ji,jj) = MINVAL(zlamtm(:,ji,jj))
         END DO
      END DO
      !-----------------------------------------------------------------------
      ! Search for boxes with only land points mark them invalid
      !-----------------------------------------------------------------------
      llinvalidcell(:,:) = .FALSE.
      DO jj = 1, jlat-1
         DO ji = 1, jlon-1
            llinvalidcell(ji,jj) =               &
               & zmskg(ji  ,jj  ) == 0.0_8 .AND. &
               & zmskg(ji+1,jj  ) == 0.0_8 .AND. &
               & zmskg(ji+1,jj+1) == 0.0_8 .AND. &
               & zmskg(ji  ,jj+1) == 0.0_8
         END DO
      END DO

      !------------------------------------------------------------------------
      ! Master loop for grid search
      !------------------------------------------------------------------------

      DO jo = 1+joffset, kobs, jostride

         !---------------------------------------------------------------------
         ! Ensure that all observation longtiudes are between 0 and 360
         !---------------------------------------------------------------------

         IF ( zplam(jo) <   0.0_8 ) zplam(jo) = zplam(jo) + 360.0_8
         IF ( zplam(jo) > 360.0_8 ) zplam(jo) = zplam(jo) - 360.0_8

         !---------------------------------------------------------------------
         ! Find observations which are on within 1e-6 of a grid point
         !---------------------------------------------------------------------

         gridloop: DO jj = 1, jlat-1
            DO ji = 1, jlon-1
               IF ( ABS( zphig(ji,jj) - zpphi(jo) ) < 1e-6 )  THEN
                  zlam = zlamg(ji,jj)
                  IF ( zlam <   0.0_8 ) zlam = zlam + 360.0_8
                  IF ( zlam > 360.0_8 ) zlam = zlam - 360.0_8
                  IF ( ABS( zlam - zplam(jo) ) < 1e-6 ) THEN
                     IF ( llinvalidcell(ji,jj) ) THEN
                        kproc(jo) = kmyproc + 1000000
                        kobsi(jo) = ji + 1
                        kobsj(jo) = jj + 1
                        CYCLE
                     ELSE
                        kproc(jo) = kmyproc
                        kobsi(jo) = ji + 1
                        kobsj(jo) = jj + 1
                        EXIT gridloop
                     ENDIF
                  ENDIF
               ENDIF
            END DO
         END DO gridloop
         
         !---------------------------------------------------------------------
         ! Ensure that all observation longtiudes are between -180 and 180
         !---------------------------------------------------------------------

         IF ( zplam(jo) > 180 ) zplam(jo) = zplam(jo) - 360.0_8

         !---------------------------------------------------------------------
         ! Do coordinate search using brute force.
         ! - For land points kproc is set to number of the processor + 1000000
         !   and we continue the search.
         ! - For ocean points kproc is set to the number of the processor 
         !   and we stop the search.
         !---------------------------------------------------------------------

         IF ( kproc(jo) == -1 ) THEN

            ! Normal case
            gridpoints : DO jj = 1, jlat-1
               DO ji = 1, jlon-1
                  
                  IF ( ( zplam(jo) > zlamtmax(ji,jj) ) .OR. &
                     & ( zplam(jo) < zlamtmin(ji,jj) ) ) CYCLE
                  
                  IF ( ABS( zpphi(jo) ) < 85 ) THEN
                     IF ( ( zpphi(jo) > zphitmax(ji,jj) ) .OR. &
                        & ( zpphi(jo) < zphitmin(ji,jj) ) ) CYCLE
                  ENDIF
                  
                  IF ( linquad( zplam(jo), zpphi(jo), &
                     &          zlamtm(:,ji,jj), zphitm(:,ji,jj) ) ) THEN
                     IF ( llinvalidcell(ji,jj) ) THEN
                        kproc(jo) = kmyproc + 1000000
                        kobsi(jo) = ji + 1
                        kobsj(jo) = jj + 1
                        CYCLE
                     ELSE
                        kproc(jo) = kmyproc
                        kobsi(jo) = ji + 1
                        kobsj(jo) = jj + 1
                        EXIT gridpoints
                     ENDIF
                  ENDIF
                  
               END DO
            END DO gridpoints

         ENDIF

         ! In case of failure retry for obs. longtiude + 360.
         IF ( kproc(jo) == -1 ) THEN
            gridpoints_greenwich : DO jj = 1, jlat-1
               DO ji = 1, jlon-1

                  IF ( ( zplam(jo)+360.0_8 > zlamtmax(ji,jj) ) .OR. &
                     & ( zplam(jo)+360.0_8 < zlamtmin(ji,jj) ) ) CYCLE

                  IF ( ABS( zpphi(jo) ) < 85 ) THEN
                     IF ( ( zpphi(jo) > zphitmax(ji,jj) ) .OR. &
                        & ( zpphi(jo) < zphitmin(ji,jj) ) ) CYCLE
                  ENDIF

                  IF ( linquad( zplam(jo)+360.0_8, zpphi(jo), &
                     &          zlamtm(:,ji,jj), zphitm(:,ji,jj) ) ) THEN
                     IF ( llinvalidcell(ji,jj) ) THEN
                        kproc(jo) = kmyproc + 1000000
                        kobsi(jo) = ji + 1
                        kobsj(jo) = jj + 1
                        CYCLE
                     ELSE
                        kproc(jo) = kmyproc
                        kobsi(jo) = ji + 1
                        kobsj(jo) = jj + 1
                        EXIT gridpoints_greenwich
                     ENDIF
                  ENDIF

               END DO
            END DO gridpoints_greenwich

         ENDIF
      END DO

      !----------------------------------------------------------------------
      ! Synchronize kproc on all processors
      !----------------------------------------------------------------------
      ! IF ( ln_grid_global ) THEN

      ! ELSE
      !    CALL obs_mpp_find_obs_proc( kproc, kobs )
      ! ENDIF

      ! WHERE( kproc(:) >= 1000000 )
      !    kproc(:) = kproc(:) - 1000000
      ! END WHERE

      DEALLOCATE( &
         & zlamg,         &
         & zphig,         &
         & zmskg,         &
         & zphitmax,      &
         & zphitmin,      &
         & zlamtmax,      &
         & zlamtmin,      &
         & llinvalidcell, &
         & zlamtm,        &
         & zphitm,        &
         & zplam,         &
         & zpphi          &
         & )

   END SUBROUTINE obs_grd_bruteforce

   LOGICAL FUNCTION linquad( px, py, pxv, pyv )
      !!----------------------------------------------------------------------
      !!                    ***  FUNCTION linquad ***
      !!
      !! ** Purpose : Determine whether a point P(x,y) lies within or on the
      !!              boundary of a quadrangle (ABCD) of any shape on a plane.
      !!
      !! ** Method  : Check if the vectorial products PA x PC, PB x PA, 
      !!              PC x PD, and PD x PB are all negative.
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  2001-11  (N. Daget, A. Weaver)
      !!        !  2006-08  (A. Weaver) NEMOVAR migration
      !!        !  2006-10  (A. Weaver) Cleanup
      !!----------------------------------------------------------------------

      !! * Arguments
      REAL(KIND=8), INTENT(IN) :: px        ! (lon) of the point P(x,y) 
      REAL(KIND=8), INTENT(IN) :: py        ! (lat) of the point P(x,y)               
      REAL(KIND=8), DIMENSION(4), INTENT(IN) :: &
         & pxv,  &                  ! (lon, lat) of the surrounding cell
         & pyv                     
  
      !! * Local declarations
      REAL(KIND=8) :: zst1
      REAL(KIND=8) :: zst2
      REAL(KIND=8) :: zst3
      REAL(KIND=8) :: zst4

      !-----------------------------------------------------------------------
      ! Test to see if the point is within the cell
      !-----------------------------------------------------------------------
      linquad = .FALSE.
      zst1 =   ( px - pxv(1) ) * ( py - pyv(4) ) &
         &   - ( py - pyv(1) ) * ( px - pxv(4) )
      IF ( zst1 <= 0.0 ) THEN
         zst2 =   ( px - pxv(4) ) * ( py - pyv(3) ) &
         &   - ( py - pyv(4) ) * ( px - pxv(3) )
         IF ( zst2 <= 0.0 ) THEN
            zst3 =   ( px - pxv(3) ) * ( py - pyv(2) ) &
               &   - ( py - pyv(3) ) * ( px - pxv(2) )
            IF ( zst3 <= 0.0) THEN
               zst4 =   ( px - pxv(2) ) * ( py - pyv(1) ) &
                  &   - ( py - pyv(2) ) * ( px - pxv(1) )
               IF ( zst4 <= 0.0 ) linquad = .TRUE.
            ENDIF
         ENDIF
      ENDIF

   END FUNCTION linquad

end module mod_gridsearch