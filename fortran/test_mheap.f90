! Copyright (c) 2017, Omar Richardson
! All rights reserved.

! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:

! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.

! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

program test_mheap
   use mheap
   implicit none
   integer (kind=4), dimension(0:1) :: node
   integer (kind=4), parameter :: nx = 4
   integer (kind=4), parameter :: ny = 5
   real (kind=8), dimension(0:nx-1,0:ny-1) :: field
   integer (kind=4) :: cap, heap_length, tree_length, i,j,k
   integer (kind=4), allocatable, dimension(:,:) :: heap
   integer (kind=4), allocatable, dimension(:) :: indx

   ! Some field according to which the arrays should be indexed.
   do i=0,nx-1
       do j=0,ny-1
           field(i,j) = 1./(i+j+1)
       end do
   end do

   cap = 10
   ! init a heap with the array comparison value, that
   ! compares the nodes' first component to order the heap
   call heap_init(heap,indx,cap)

   ! insert some data
   call heap_insert(heap,indx,cap,heap_length,tree_length, [1,2],field,nx,ny )
   call heap_insert(heap,indx,cap,heap_length,tree_length, [2,0],field,nx,ny )
   call heap_insert(heap,indx,cap,heap_length,tree_length, [3,1],field,nx,ny )
   call heap_insert(heap,indx,cap,heap_length,tree_length, [1,3],field,nx,ny )
   call heap_insert(heap,indx,cap,heap_length,tree_length, [0,4],field,nx,ny )
   call heap_insert(heap,indx,cap,heap_length,tree_length, [3,2],field,nx,ny )

   ! data is kept unordered (except for the root node)
   write(*,*)
   write(*,*) 'unordered traversal, the first element is always the root node'
   do k = 0, heap_length-1
      call heap_peek( heap, indx,cap,heap_length,k,node) 
      write(*,*) node
   enddo

   ! when the root node is popped, a new root node is set.
   ! to traverse it in order just pop all the root elements.
   write(*,*)
   write(*,*) 'ordered traversal'
   do k = 0, heap_length-1
      call heap_pop(heap, indx, cap, heap_length, field, nx, ny, node)
      write(*,*) node
   enddo
   write(*,*) 'now heap is empty:', heap_length

   ! data is not lost from the tree when using pop,
   ! we can reheap the whole tree and start over
   write(*,*) 'tree data is kept after popping, so as long'
   write(*,*) 'as no new insertions are made, the same data'
   write(*,*) 'can be reheaped using the same or another function.'
   write(*,*)
   write(*,*) 'reheap using same function'
   call heap_reheap(heap,indx, cap, heap_length, tree_length, field, nx, ny)
   do k = 0, heap_length-1
      call heap_pop(heap, indx, cap, heap_length, field, nx, ny, node)
      write(*,*) node
   enddo
   write(*,*) 'insert 2 new elements'
   call heap_insert(heap,indx,cap,heap_length,tree_length, [3,4],field,nx,ny )
   call heap_insert(heap,indx,cap,heap_length,tree_length, [1,1],field,nx,ny )
   write(*,*) 'pop all'
   do k = 0, heap_length-1
      call heap_pop(heap, indx, cap, heap_length, field, nx, ny, node)
      write(*,*) node
   enddo

end program test_mheap
