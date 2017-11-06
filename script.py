import numpy as np

try:
    from fortran_modules.wdt_fortran import propagate_dist as call_propagate_distance

    fortran_lib = True
except ImportError:
    print("No Fortran modules found, falling back on python implementation")
    fortran_lib = False
import math
import heapq
from scipy.misc import imread
import time
import matplotlib.pyplot as plt


class WDT:
    DIR_STRINGS = ["left", "down", "right", "up"]
    DIRS = ((-1, 0), (0, -1), (1, 0), (0, 1))

    def __init__(self, filename='ex2.png'):
        data = imread(filename, mode='RGB')
        self.cost_field = WDT.convert_image_to_cost_field(data)
        self.nx, self.ny = self.cost_field.shape
        self.weighted_distance_transform = None

    @staticmethod
    def convert_image_to_cost_field(data):
        """
        Read image data and convert it to a *landscape*, a 2D array with cost fields.
        This cost field forms the input for the weighted distance transform
        zero costs denote exits, infinite costs denote obstacles.
        For now, we follow previous Mercurial standards: obstacles are in black, exits in red.
        Accessible space is in white, less accessible space has less white.
        :param data: array with RGB image data
        :return: 2D array representing the cost field
        """

        # Exits are present in all red enough places ("R >> BG")
        exits = np.where(data[:, :, 0] - (data[:, :, 1] + data[:, :, 2]) / 2 > 160)
        # Obstacles are in black (so at least G and B must be zero)
        obstacles = np.where(data[:, :, 1] + data[:, :, 2] == 0)
        # Convert image to greyscale
        grey_scales = np.dot(data[..., :3], [0.299, 0.587, 0.114])
        # Boolean index array for places without exits and obstacles
        space = np.ones(grey_scales.shape, dtype=np.bool)
        space[exits] = False
        space[obstacles] = False
        # Cost field: Inversely proportional to greyscale values
        cost_field = np.empty(data[:, :, 0].shape)
        cost_field[obstacles] = np.inf
        cost_field[exits] = 0
        cost_field[space] = 1. / grey_scales[space]
        return cost_field

    def _exists(self, index):
        """
        Checks whether an index exists an array
        :param index: 2D index tuple
        :return: true if lower than tuple, false otherwise
        """
        return (0 <= index[0] < self.nx) and (0 <= index[1] < self.ny)

    def get_new_candidate_cells(self, cell, unknown_cells):
        """
        Compute the new candidate cells (cells for which we have no definite distance value yet
        For more information on the algorithm: check fast marching method
        :param cell: tuple of index; a new cell that has been added to the distance field
        :param unknown_cells: set of tuples; all cells still unknown
        :return: Set of new candidate cells for which to compute the distance
        """
        new_candidate_cells = set()
        for direction in WDT.DIRS:
            nb_cell = (cell[0] + direction[0], cell[1] + direction[1])
            if nb_cell in unknown_cells:
                new_candidate_cells.add(nb_cell)
        return new_candidate_cells

    def propagate_distance(self, cell, costs):
        """
        Compute the weighted distance in a cell using costs and distances in other cells
        :param cell: tuple, index of a candidate cell
        :param costs: list of cost arrays in X and Y direction
        :param wdt_array: the weighted distance transform field up until now
        :return: a approximate distance based on the neighbour cells
        """
        # Find the minimal directions along a grid cell.
        # Assume left and below are best, then overwrite with right and up if they are better
        adjacent_distances = np.ones(4) * np.inf
        pots_from_axis = [0, 0]  # [x direction, y direction]
        costs_from_axis = [np.inf, np.inf]  #
        for i, dir_s in enumerate(WDT.DIR_STRINGS):
            # Direction for which we check the cost
            normal = WDT.DIRS[i]
            nb_cell = (cell[0] + normal[0], cell[1] + normal[1])
            if not self._exists(nb_cell):
                continue
            pot = self.weighted_distance_transform[nb_cell]
            # distance in that neighbour field
            if dir_s == 'left':
                face_index = (nb_cell[0] + 1, nb_cell[1])
            elif dir_s == 'down':
                face_index = (nb_cell[0], nb_cell[1] + 1)
            else:
                face_index = nb_cell
            # Left/right is x, up/down is y
            cost = costs[i % 2][face_index]
            # Proposed cost along this direction
            adjacent_distances[i] = pot + cost
            # If it is cheaper to go from the opposite direction
            if adjacent_distances[i] < adjacent_distances[(i + 2) % 4]:
                pots_from_axis[i % 2] = pot
                costs_from_axis[i % 2] = cost
            hor_pot, ver_pot = pots_from_axis
            hor_cost, ver_cost = costs_from_axis
            # Coefficients of quadratic equation (upwind discretization)
        a = 1. / hor_cost ** 2 + 1. / ver_cost ** 2
        b = -2 * (hor_pot / hor_cost ** 2 + ver_pot / ver_cost ** 2)
        c = (hor_pot / hor_cost) ** 2 + (ver_pot / ver_cost) ** 2 - 1

        D = b ** 2 - 4 * a * c
        x_high = (2 * c) / (-b - math.sqrt(D))
        # Might not be obvious, but why we take the largest root is found in report.
        return x_high

    def get_weighted_distance_transform(self):
        """
        Compute the weighted distance transform with cost/image function u using a fast marching algorithm
        We compute the distance transform with costs defined on a staggered grid for consistency.
        This means that we use costs defined on the faces of cells, found by averaging values of the adjacent cells.

        Starting from the exit, we march over each pixel with the lowest weighted distance,
        until we found values for all pixels.
        :param cost_field: nonnegative 2D array with cost in each cell/pixel, zero and infinity are allowed values.
        :return: weighted distance transform
        """

        # Cost for moving along horizontal lines
        costs_x = np.ones([self.nx + 1, self.ny], order='F') * np.inf
        costs_x[1:-1, :] = (self.cost_field[1:, :] + self.cost_field[:-1, :]) / 2
        # Cost for moving along vertical lines
        costs_y = np.ones([self.nx, self.ny + 1], order='F') * np.inf
        costs_y[:, 1:-1] = (self.cost_field[:, 1:] + self.cost_field[:, :-1]) / 2

        # Initialize locations (known/unknown/exit/obstacle)
        self.weighted_distance_transform = np.ones_like(self.cost_field, order='F') * np.inf
        exit_locs = np.where(self.cost_field == 0)
        obstacle_locs = np.where(self.cost_field == np.inf)
        self.weighted_distance_transform[exit_locs] = 0

        # Initialize Cell structures
        all_cells = {(i, j) for i in range(self.nx) for j in range(self.ny)}
        known_cells = {cell for cell in zip(exit_locs[0], exit_locs[1])}
        unknown_cells = all_cells - known_cells - {cell for cell in zip(obstacle_locs[0], obstacle_locs[1])}
        new_candidate_cells = set()
        for cell in known_cells:
            new_candidate_cells |= self.get_new_candidate_cells(cell, unknown_cells)
        cand_heap = [(np.inf, cell) for cell in new_candidate_cells]
        # Loop until all unknown cells have a distance value
        while unknown_cells:
            # by repeatedly looping over the new candidate cells
            for cell in new_candidate_cells:
                # Compute a distance for each cell based on its neighbour cells
                if fortran_lib:
                    distance = call_propagate_distance(cell[0], cell[1], self.nx, self.ny,
                                                       self.weighted_distance_transform, costs_x, costs_y, 99999)
                else:
                    distance = self.propagate_distance(cell, [costs_x, costs_y])
                # Store this value in the dictionary, and in the heap (for fast lookup)
                # Don't check whether we have the distance already in the heap; check on outcome
                heapq.heappush(cand_heap, (distance, cell))
            # See if the heap contains a good value and if so, add it to the field. If not, finish.
            # Since we can store multiple distance values for one cell, we might need to pop a couple of times
            while True:
                min_distance, best_cell = heapq.heappop(cand_heap)
                if self.weighted_distance_transform[best_cell] == np.inf:
                    # Got a good one: no assigned distance in wdt yet
                    break
                elif min_distance == np.inf:  # No more finite values; done
                    return self.weighted_distance_transform
            # Good value found, add to the wdt and
            self.weighted_distance_transform[best_cell] = min_distance
            unknown_cells.remove(best_cell)
            new_candidate_cells = self.get_new_candidate_cells(best_cell, unknown_cells)
        return self.weighted_distance_transform


wdt = WDT()
phi = wdt.get_weighted_distance_transform()
plt.imshow(phi)
plt.show()
