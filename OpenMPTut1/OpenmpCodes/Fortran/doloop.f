      program doloop
      integer omp_get_thread_num 
!$omp threadprivate (ia) 
      save ia
   
      call omp_set_num_threads(4)
!$omp parallel private(id,i)
      id = omp_get_thread_num()

!$omp do 
      do i=1,8
      ia=id+i
      write(*,*)'ia=',ia, 'for i=',i, 'from id=',id
      end do   
!$omp end do
!$omp do
      do i=1,9
      ia=ia+id
      write(*,*)'updated ia=',ia, 'for i=',i, 'from id=',id
      end do
!$omp end do

!$omp end parallel
      end

