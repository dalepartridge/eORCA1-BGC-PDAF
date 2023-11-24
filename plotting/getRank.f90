!< get the rank of observation in the model foreacst
!< ensembles
subroutine getRank(nlength, Ne, nsize, x, rank)
   integer, intent(in) :: nlength, Ne, nsize
   integer, intent(in) :: x(nsize)
   integer, intent(inout) :: rank(nlength)
   integer :: i, j
   do i = 1, nlength
      do j = 1, Ne+1
         if (x(j + (i-1)*(Ne+1)) == 1) then
               rank(i) = j
               exit
         end if
      end do
   end do
end subroutine getRank
