import numpy as np
from pot import compute_potential
import math
import heapq
from scipy.misc import imread
import time

DIR_STRINGS = ["left", "down", "right", "up"]
AXES = ['x', 'y']
DIRS = np.array([[-1, 0], [0, -1], [1, 0], [0, 1]])


def read_image(filename='ex5.png'):
    data = imread(filename, mode='RGB')
    exits = np.where(data[:, :, 0] - data[:, :, 1] > 200)
    obstacles = np.where(data[:, :, 1] < 100)
    u = np.ones(data[:, :, 0].shape)
    u[obstacles] = np.inf
    u[exits] = 0
    return u


# TODO: Change "potentials" to "distances"
def exists(index):
    """
    Checks whether an index exists an array
    :param index: 2D index tuple
    :param max_index: max index tuple
    :return: true if lower than tuple, false otherwise
    """
    return (0 <= index[0] < nx) and (0 <= index[1] < ny)


def get_new_candidate_cells(new_known_cells, unknown_cells):
    new_candidate_cells = set()
    for cell in new_known_cells:
        for direction in DIRS:
            nb_cell = (cell[0] + direction[0], cell[1] + direction[1])
            if nb_cell in unknown_cells:
                new_candidate_cells.add(nb_cell)
    return new_candidate_cells


def compute_potehntial(cell, costs, potential):
    # Find the minimal directions along a grid cell.
    # Assume left and below are best, then overwrite with right and up if they are better
    adjacent_potentials = np.ones(4) * np.inf
    pots_from_axis = [0, 0]  # [x direction, y direction]
    costs_from_axis = [np.inf, np.inf]  #
    for i, dir_s in enumerate(DIR_STRINGS):
        # Direction for which we check the cost
        normal = DIRS[i]
        nb_cell = (cell[0] + normal[0], cell[1] + normal[1])
        if not exists(nb_cell):
            continue
        pot = potential[nb_cell]
        # potential in that neighbour field
        if dir_s == 'left':
            face_index = (nb_cell[0] + 1, nb_cell[1])
        elif dir_s == 'down':
            face_index = (nb_cell[0], nb_cell[1] + 1)
        else:
            face_index = nb_cell
        # Left/right is x, up/down is y
        cost = costs[i % 2][face_index]
        # Proposed cost along this direction
        adjacent_potentials[i] = pot + cost
        # If it is cheaper to go from the opposite direction
        if adjacent_potentials[i] < adjacent_potentials[(i + 2) % 4]:
            pots_from_axis[i % 2] = pot
            costs_from_axis[i % 2] = cost
        hor_pot, ver_pot = pots_from_axis
        hor_cost, ver_cost = costs_from_axis
        # Coefficients of quadratic equation (upwind discretization)
    a = 1 / hor_cost ** 2 + 1 / ver_cost ** 2
    b = -2 * (hor_pot / hor_cost ** 2 + ver_pot / ver_cost ** 2)
    c = (hor_pot / hor_cost) ** 2 + (ver_pot / ver_cost) ** 2 - 1

    D = b ** 2 - 4 * a * c
    x_high = (2 * c) / (-b - math.sqrt(D))
    # Might not be obvious, but why we take the largest root is found in report.
    return x_high


def compute_distance_transform(u):
    """
    Compute the weighted distance transform with cost/image function u using a fast marching algorithm
    We compute the distance transform on a staggered grid: this is more precise
    :param u: nonnegative 2D array with cost in each cell/pixel, infinity is allowed.
    :return: weighted distance transform

    """
    # nx,ny = u.shape

    # Cost for moving along horizontal lines
    u_x = np.ones([nx + 1, ny], order='F') * np.inf
    u_x[1:-1, :] = (u[1:, :] + u[:-1, :]) / 2
    # Cost for moving along vertical lines
    u_y = np.ones([nx, ny + 1], order='F') * np.inf
    u_y[:, 1:-1] = (u[:, 1:] + u[:, :-1]) / 2

    # Initialize locations (known/unknown/exit/obstacle)
    phi = np.ones_like(u, order='F') * np.inf
    exit_locs = np.where(u == 0)
    obstacle_locs = np.where(u == np.inf)
    phi[exit_locs] = 0

    # Initialize Cell structures
    all_cells = {(i, j) for i in range(nx) for j in range(ny)}
    known_cells = {cell for cell in zip(exit_locs[0], exit_locs[1])}
    unknown_cells = all_cells - known_cells - {cell for cell in zip(obstacle_locs[0], obstacle_locs[1])}
    new_candidate_cells = get_new_candidate_cells(known_cells, unknown_cells)
    candidate_cells = {cell: np.inf for cell in new_candidate_cells}
    cand_heap = [(np.inf, cell) for cell in candidate_cells]
    while unknown_cells:
        for cell in new_candidate_cells:
            if True:
                potential = compute_potential(cell[0], cell[1], nx, ny, phi, u_x, u_y, 99999)
            else:
                potential = compute_potehntial(cell, [u_x, u_y], phi)
            candidate_cells[cell] = potential
            # Don't check whether we have the potential already in the heap; check on outcome
            heapq.heappush(cand_heap, (potential, cell))
        popped_new_potential = False
        while not popped_new_potential:
            min_potential, best_cell = heapq.heappop(cand_heap)
            if phi[best_cell] == np.inf:
                popped_new_potential = True
        candidate_cells.pop(best_cell)
        phi[best_cell] = min_potential
        unknown_cells.remove(best_cell)
        known_cells.add(best_cell)
        new_candidate_cells = get_new_candidate_cells({best_cell}, unknown_cells)

    # While there are unknown cells:
    # Compute the value of the candidate cells from the known cells
    # Pick the cheapest candidate cell, make it known
    # Get new candidate cells and continue while
    return phi


u = read_image()
nx, ny = u.shape
import matplotlib.pyplot as plt

time1 = time.time()
phi = compute_distance_transform(u)
time2 = time.time()
print(time2 - time1)
plt.imshow(phi)
plt.show()

# Improvements: -Sweep. -Fortran total. -Downgrade floats
