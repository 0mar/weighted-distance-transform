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

! Ideas for optimization:
! - intent in/out/inout
! - [..,..]/dimension(..,..)
! - Increase the size of the array by one so all `exists` can be tossed out
! - Check indices in looping

subroutine weighted_distance_transform(cost_field,wdt_field,n_x,n_y,obstacle_value)
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

integer (kind=4) :: i,j, exists
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
contains
    function new_candidate_cells(cell_x,cell_y) result(cand_cells)
        implicit none
        ! Find new candidate cells
        integer (kind=4) :: direction,nb_x,nb_y
        integer (kind=4), dimension(0:3,0:1) :: cand_cells
        integer (kind=4) :: cell_x,cell_y
        cand_cells = -1
        do direction=0,3
            nb_x = sign(mod(direction+1,2),direction-1) + cell_x
            nb_y = sign(mod(direction,2),direction-2) + cell_y
            if (exists(nb_x,nb_y,n_x,n_y)==1)  then
                if (wdt_field(nb_x,nb_y) == unknown_value) then
                    cand_cells(direction,0) = nb_x
                    cand_cells(direction,1) = nb_y
                end if
            end if
        end do
    end function new_candidate_cells


end subroutine

integer (kind=4) function exists(cell_x,cell_y,max_x,max_y)
! Find out whether the cell exists given the dimensions of the scene
integer (kind=4) :: cell_x,cell_y,max_x,max_y
if (cell_x < 0 .or. cell_y < 0 .or. cell_x>=max_x .or. cell_y>=max_y) then
    exists = 0
    return
end if
exists = 1
return
end function


