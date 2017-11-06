subroutine compute_potential(cell_x,cell_y,n_x,n_y,pot_field,u_x,u_y,inf,out_pot)

!   Compute the potential in a single cell with a first order upwind method
implicit none
!f2py intent(in) cell_x,cell_y,pot_field,u_x,u_y,inf
!f2py depend(n_x,n_y) pot_field,u_x,u_y
!f2py intent(out) out_pot
!
! this subroutine uses 0 based indexing
!
integer ::  n_x,n_y
real (kind=8) :: inf ! value exceeding largest value in potential computation. max(pot) + max(u_x,u_y) + 1 will do
real (kind=8) :: a,b,c,D ! parameters for upwind approximation
integer (kind=4) :: normal_x,normal_y,face_index_x,face_index_y, exists ! Administration
integer (kind=4) :: nb_cell_x,nb_cell_y! neighbour indices
real (kind=8):: hor_potential,ver_potential,hor_cost,ver_cost

integer (kind=4) :: direction
integer (kind=4), parameter :: LEFT = 0
integer (kind=4), parameter :: DOWN = 1
integer (kind=4), parameter :: RIGHT = 2
integer (kind=4), parameter :: UP = 3

real (kind=8), dimension(0:n_x-1,0:n_y-1) :: pot_field
real (kind=8), dimension(0:n_x,0:n_y-1) :: u_x
real (kind=8), dimension(0:n_x-1,0:n_y) :: u_y
real (kind=8), dimension(0:3) :: neighbour_pots
integer (kind=4) :: cell_x,cell_y
real (kind=8) :: pot, out_pot,cost

neighbour_pots=  (/ inf, inf, inf, inf /) 
hor_cost = inf
ver_cost = inf
    ! Find the minimal directions along a grid cell.
    ! Assume left and below are best, then over!write with right and up if they are better
do direction=0,3

        normal_x = sign(mod(direction+1,2),direction-1)
        normal_y = sign(mod(direction,2),direction-2)
        ! numerical direction
        nb_cell_x = cell_x + normal_x
        nb_cell_y = cell_y + normal_y
        if (exists(nb_cell_x,nb_cell_y,n_x,n_y) == 0) then
            cycle
        endif
        pot = pot_field(nb_cell_x,nb_cell_y)
        ! potential in that neighbour field
        if (direction == RIGHT) then
            face_index_x = nb_cell_x
            face_index_y = nb_cell_y
            ! Unit cost values are defined w.r.t faces, not cells! So the indexing is different with right and up.
            cost = u_x(face_index_x,face_index_y)
            ! Cost to go from there to here
        elseif (direction == UP) then
            face_index_x = nb_cell_x
            face_index_y = nb_cell_y
            ! Unit cost values are defined w.r.t faces, not cells! So the indexing is different with right and up.
            cost = u_y(face_index_x,face_index_y)
        elseif (direction == LEFT) then
            face_index_x = nb_cell_x+1
            face_index_y = nb_cell_y
            cost = u_x(face_index_x,face_index_y)
        elseif (direction == DOWN) then
            face_index_x = nb_cell_x
            face_index_y = nb_cell_y+1
            cost = u_y(face_index_x,face_index_y)
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

if (hor_cost >= inf) then
    !write(*,*) "Horizontal infinite"
    a = 1. / (ver_cost * ver_cost)
    b = -2. * ver_potential / (ver_cost * ver_cost)
    c = (ver_potential / ver_cost) * (ver_potential / ver_cost) -1
elseif (ver_cost >=inf) then
    !write(*,*) "Vertical infinite"
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

