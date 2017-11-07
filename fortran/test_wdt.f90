program test_wdt
    implicit none
    integer (kind=4), parameter :: n_x=4
    integer (kind=4), parameter :: n_y=4
    real (kind=8), dimension(0:n_x-1,0:n_y-1) :: cost_field,wdt_field
    real (kind=8), parameter :: obstacle_value = 200

    ! Example:
    !- - - -
    !* * * -
    !- * - -
    !- - - -

    cost_field = 1
    cost_field(1,1) = 0
    cost_field(0:3,2) = 0
    cost_field(3,1:3)= obstacle_value
    wdt_field = obstacle_value
    wdt_field(2,2) = 0
    call weighted_distance_transform(cost_field,wdt_field,n_x,n_y,obstacle_value)
end program

subroutine weighted_distance_transform(cost_field,wdt_field,n_x,n_y,obstacle_value)
use wdt_module
implicit none
!   Compute the potential in a single cell with a first order upwind method
!f2py intent(in) cell_x,cell_y,pot_field,u_x,u_y,inf
!f2py depend(n_x,n_y) pot_field,u_x,u_y
!f2py intent(out) out_pot
!
integer ::  n_x,n_y
real (kind=8) :: obstacle_value ! value exceeding largest in potential computation. max(pot) + max(u_x,u_y) + 1 will do
real (kind=8) :: unknown_value ! Value exceeding inf, indicated that we do not know yet
integer (kind=4) :: nb_cell_x,nb_cell_y, exists! neighbour indices
real (kind=8):: hor_potential,ver_potential,hor_cost,ver_cost

integer (kind=4) :: direction
integer (kind=4), parameter :: LEFT = 0
integer (kind=4), parameter :: DOWN = 1
integer (kind=4), parameter :: RIGHT = 2
integer (kind=4), parameter :: UP = 3

integer (kind=4), dimension(:,:), allocatable :: nb_values

integer (kind=4) :: exists
integer (kind=4) :: i,j
real (kind=8), dimension(0:n_x-1,0:n_y-1) :: cost_field,wdt_field
real (kind=8), dimension(0:n_x,0:n_y-1) :: costs_x
real (kind=8), dimension(0:n_x-1,0:n_y) :: costs_y
integer (kind=4) :: cell_x,cell_y
real (kind=8) :: pot

!neighbour_pots=  (/ inf, inf, inf, inf /) 
costs_x = obstacle_value
costs_x(1:n_x-1,:) = (cost_field(1:n_x-1,:) + cost_field(0:n_x-2,:))/2
costs_y = obstacle_value
costs_x(:,1:n_y-1) = (cost_field(:,1:n_y-1) + cost_field(:,0:n_y-2))/2

    allocate(nb_values(0:3,0:1))
        nb_values = new_candidate_cells(2,1)
    print*,nb_values
    deallocate(nb_values)

    do i=0,n_x-1
        write(*,*) ( cost_field(j,n_x-i-1), j=0,n_y-1 )
    enddo

end subroutine

