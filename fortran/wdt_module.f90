module wdt_module
    implicit none
    public :: exists, new_candidate_cells, add_candidate_cells, propagate_dist
    public :: n_x2,n_y2,obstacle_value,unknown_value
    public :: LEFT, DOWN, RIGHT, UP, NUM_DIRS, KNOWN, UNKNOWN, CANDIDATE, NEW_CANDIDATE
    integer (kind=4) :: n_x2,n_y2,obstacle_value,unknown_value

    integer (kind=4), parameter :: LEFT = 0
    integer (kind=4), parameter :: DOWN = 1
    integer (kind=4), parameter :: RIGHT = 2
    integer (kind=4), parameter :: UP = 3
    integer (kind=4), parameter :: NUM_DIRS = 4

    integer (KIND=4), parameter :: KNOWN = 0
    integer (KIND=4), parameter :: UNKNOWN = 1
    integer (KIND=4), parameter :: CANDIDATE = 2
    integer (KIND=4), parameter :: NEW_CANDIDATE = 3


contains

    function new_candidate_cells(cell_x,cell_y,wdt_field) result(cand_cells)
        implicit none
        ! Find new candidate cells
        integer (kind=4), intent(in) :: cell_x,cell_y
        real (kind=8), dimension(0:n_x2-1,0:n_y2-1), intent(in) :: wdt_field
        integer (kind=4), dimension(0:3,0:1) :: cand_cells !DON'T mark intent out
        integer (kind=4) :: direction,nb_x,nb_y
        cand_cells = -1
        do direction=0,3
            nb_x = sign(mod(direction+1,2),direction-1) + cell_x
            nb_y = sign(mod(direction,2),direction-2) + cell_y
            if (exists(nb_x,nb_y)==1)  then
                if (wdt_field(nb_x,nb_y) == unknown_value) then
                    cand_cells(direction,0) = nb_x
                    cand_cells(direction,1) = nb_y
                end if
            end if
        end do
    end function new_candidate_cells
    subroutine add_candidate_cells(cell_x,cell_y, cell_indicators)
        implicit none
        integer (kind=4) :: cell_x,cell_y
        integer (kind=4), dimension(0:n_x2-1,0:n_y2-1) :: cell_indicators
        integer (kind=4) :: direction
        integer (kind=4) :: nb_x,nb_y

        do direction=0,3
            nb_x = sign(mod(direction+1,2),direction-1) + cell_x
            nb_y = sign(mod(direction,2),direction-2) + cell_y
            if (exists(nb_x,nb_y)==1)  then
				if (cell_indicators(nb_x,nb_y)>=UNKNOWN) then
                    cell_indicators(nb_x,nb_y) = NEW_CANDIDATE
                end if
            end if
        end do
        end subroutine add_candidate_cells

    subroutine propagate_dist(cell_x,cell_y,wdt_field,costs_x,costs_y,out_pot)
    ! Todo: Can we make it faster by passing just the relevant array?
    !   Compute the potential in a single cell with a first order upwind method
    implicit none
    real (kind=8) :: a,b,c,D ! parameters for upwind approximation
    integer (kind=4) :: normal_x,normal_y,face_index_x,face_index_y! Administration
    integer (kind=4) :: nb_cell_x,nb_cell_y! neighbour indices
    real (kind=8):: hor_potential,ver_potential,hor_cost,ver_cost

    integer (kind=4) :: direction

    real (kind=8), dimension(0:n_x2-1,0:n_y2-1) :: wdt_field
    real (kind=8), dimension(0:n_x2,0:n_y2-1) :: costs_x
    real (kind=8), dimension(0:n_x2-1,0:n_y2) :: costs_y
    real (kind=8), dimension(0:3) :: neighbour_pots
    integer (kind=4) :: cell_x,cell_y
    real (kind=8) :: pot, out_pot,cost

    neighbour_pots=  (/ obstacle_value, obstacle_value, obstacle_value, obstacle_value /) 
    hor_cost = obstacle_value
    ver_cost = obstacle_value
        ! Find the minimal directions along a grid cell.
        ! Assume left and below are best, then over!write with right and up if they are better
    do direction=0,3
            normal_x = sign(mod(direction+1,2),direction-1)
            normal_y = sign(mod(direction,2),direction-2)
            ! numerical direction
            nb_cell_x = cell_x + normal_x
            nb_cell_y = cell_y + normal_y
            if (exists(nb_cell_x,nb_cell_y) == 0) then
                cycle
            endif
            pot = wdt_field(nb_cell_x,nb_cell_y)
            ! potential in that neighbour field
            if (direction == RIGHT) then
                face_index_x = nb_cell_x
                face_index_y = nb_cell_y
                ! Unit cost values are defined w.r.t faces, not cells! So the indexing is different with right and up.
                cost = costs_x(face_index_x,face_index_y)
                ! Cost to go from there to here
            elseif (direction == UP) then
                face_index_x = nb_cell_x
                face_index_y = nb_cell_y
                ! Unit cost values are defined w.r.t faces, not cells! So the indexing is different with right and up.
                cost = costs_y(face_index_x,face_index_y)
            elseif (direction == LEFT) then
                face_index_x = nb_cell_x+1
                face_index_y = nb_cell_y
                cost = costs_x(face_index_x,face_index_y)
            elseif (direction == DOWN) then
                face_index_x = nb_cell_x
                face_index_y = nb_cell_y+1
                cost = costs_y(face_index_x,face_index_y)
            else
                !write(*,*) "Exception in compute_potential"
            endif
             
            neighbour_pots(direction) = pot + cost
            ! total potential
            if (neighbour_pots(direction) < neighbour_pots(mod(direction+2,4))) then
                if (mod(direction,2) == 0) then
                    hor_potential = pot
                    hor_cost = cost
                    !write(*,*) "Horizontal changed"
                    ! lowest in horizontal direction
                else
                    ver_potential = pot
                    ver_cost = cost
                    !write(*,*) "Vertical changed"
                    ! lowest in vertical direction
                endif
            endif
    end do

    if (hor_cost >= obstacle_value) then
        !write(*,*) "Horizontal obstacle_valueinite"
        a = 1. / (ver_cost * ver_cost)
        b = -2. * ver_potential / (ver_cost * ver_cost)
        c = (ver_potential / ver_cost) * (ver_potential / ver_cost) -1
    elseif (ver_cost >=obstacle_value) then
        !write(*,*) "Vertical obstacle_valueinite"
        a = 1. / (hor_cost * hor_cost)
        b = -2. * hor_potential / (hor_cost * hor_cost)
        c = (hor_potential / hor_cost) * (hor_potential / hor_cost) -1
    else
        !write(*,*) "All good"
        a = 1. / (hor_cost * hor_cost) + 1. / (ver_cost * ver_cost)
        b = -2. * (hor_potential / (hor_cost * hor_cost) + ver_potential / (ver_cost * ver_cost))
        c = (hor_potential / hor_cost) * (hor_potential / hor_cost) + (ver_potential / ver_cost) * (ver_potential / ver_cost) - 1
    endif

    D = b*b-4.*a*c

    out_pot = (-b + sqrt(D)) / (2.*a)
    end subroutine


    function exists(cell_x,cell_y) result(ex)
        implicit none
        ! Find out whether the cell exists given the dimensions of the scene
        integer (kind=4), intent(in) :: cell_x,cell_y
        integer (kind=4) :: ex
        if (cell_x < 0 .or. cell_y < 0 .or. cell_x>=n_x2 .or. cell_y>=n_y2) then
            ex = 0
        else
            ex = 1
        end if
    end function exists
end module wdt_module

! TODO: Ideas for optimization:
! - intent in/out/inout
! - [..,..]/dimension(..,..)
! - Increase the size of the array by one so all `exists` can be tossed out
! - Check indices in looping

program test_wdt
    implicit none
    integer (kind=4), parameter :: n_x=4
    integer (kind=4), parameter :: n_y=4
    integer (kind=4) :: i,j
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
    do i=0,n_y-1
        write(*,*) ( wdt_field(j,n_x-i-1), j=0,n_x-1 )
    enddo
end program

subroutine weighted_distance_transform(cost_field,wdt_field,n_x,n_y,obs_val)
use wdt_module
implicit none
!   Compute the potential in a single cell with a first order upwind method
!f2py intent(in) cost_field,obs_val,n_x,n_y
!f2py depend(n_x,n_y) cost_field,wdt_field
!f2py intent(out) wdt_field
!
integer (kind=4) :: nb_cell_x,nb_cell_y,n_x,n_y
real (kind=8):: hor_potential,ver_potential,hor_cost,ver_cost,obs_val

integer (kind=4) :: direction
integer (kind=4), dimension(:,:), allocatable :: nb_values

integer (kind=4) :: i,j, min_i,min_j
real (kind=8), dimension(0:n_x-1,0:n_y-1) :: cost_field,wdt_field
real (kind=8), dimension(0:n_x,0:n_y-1) :: costs_x
real (kind=8), dimension(0:n_x-1,0:n_y) :: costs_y
integer (kind=4), dimension(0:n_x-1,0:n_y-1) :: cell_indicators
integer (kind=4) :: cell_x,cell_y
integer (kind=4) :: num_cand_cells
real (kind=8) :: dist, min_dist

n_x2 = n_x
n_y2 = n_y
obstacle_value = obs_val
unknown_value = obs_val + 1
! Cost for moving along horizontal lines
! Todo: I don't think we ever access the outer boundary of the costs_x/y array
costs_x = obstacle_value
costs_x(1:n_x-1,:) = (cost_field(1:n_x-1,:) + cost_field(0:n_x-2,:))/2
! Cost for moving along vertical lines
costs_y = obstacle_value
costs_y(:,1:n_y-1) = (cost_field(:,1:n_y-1) + cost_field(:,0:n_y-2))/2

!Initialize locations (known(exit/obstacles)/unknown)
wdt_field = unknown_value !TODO: Faster to do in loop? Check
cell_indicators = UNKNOWN ! Same!
do i=0,n_x2-1
    do j=0,n_y2-1
        if (cost_field(i,j)==0) then
            ! No cost, so this is an exit
            wdt_field(i,j) = 0
            cell_indicators(i,j) = KNOWN
            ! ! Count the maximum number of candidate cells to allocate later
            ! num_cand_cells = num_cand_cells + NUM_DIRS
            call add_candidate_cells(i,j,cell_indicators)
        elseif (cost_field(i,j)>=obstacle_value) then
            ! 'infinite cost', so this is an obstacle
            wdt_field(i,j) = obstacle_value
            cell_indicators(i,j) = KNOWN
        endif
    enddo
enddo
    do i=0,n_y-1
        write(*,*) ( cell_indicators(j,n_x-i-1), j=0,n_x-1 )
    enddo
! allocate(nb_values(0:num_cand_cells-1,0:1))

! nb_values = new_candidate_cells(2,1,wdt_field)
! print*,nb_values
! deallocate(nb_values)
do while (.true.)
    ! Find the minimal distance to add to the field
    min_dist = obstacle_value
    min_i = -1
    min_j = -1
    do i=0,n_x2-1
        do j=1,n_y2-1
            if (cell_indicators(i,j)>=CANDIDATE) then
                ! Compute distances for all new candidate cells
                if (cell_indicators(i,j) == NEW_CANDIDATE) then
                    call propagate_dist(i,j,wdt_field,costs_x,costs_y,dist)
                    if (dist < wdt_field(i,j)) then
                        ! If the distance is lower, overwrite
                        wdt_field(i,j) = dist
                    end if
                    ! All new candidates now become regular candidate cells
                    cell_indicators(i,j) = CANDIDATE
                end if
                if (min_dist > wdt_field(i,j)) then
                    ! Minimize the distance over all candidate cells
                    min_dist = wdt_field(i,j)
                    min_i = i
                    min_j = j
                end if
            end if
        end do
    end do
    if (min_dist == obstacle_value) then
        exit
    else
    ! Mark the cell as known
        cell_indicators(min_i,min_j) = KNOWN
        call add_candidate_cells(min_i,min_j,cell_indicators)
    end if
end do

end subroutine

