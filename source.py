    def compute_potential_field(self):
        """
        Compute the potential field as a function of the unit cost.
        Also computes the gradient of the potential field
        Implemented using the fast marching method

        Potential is initialized with zero on exits, and a fixed high value on inaccessible cells.
        :return:
        """
        opposites = {'left': 'right', 'right': 'left', 'up': 'down', 'down': 'up'}
        # This implementation is allowed to be naive: it's costly and should be implemented in FORTRAN
        # But maybe do a heap structure first?
        potential_field = self.initial_interface.copy()
        known_cells = self.exit_cell_set.copy()
        unknown_cells = self.all_cells - known_cells - self.obstacle_cell_set
        # All the inaccessible cells are not required.

        def get_new_candidate_cells(new_known_cells):  # Todo: Finalize list/set interfacing
            new_candidate_cells = set()
            for cell in new_known_cells:
                for direction in ft.DIRECTIONS.values():
                    nb_cell = (cell[0] + direction[0], cell[1] + direction[1])
                    if self._exists(nb_cell) and nb_cell not in known_cells and nb_cell not in self.obstacle_cell_set:
                        new_candidate_cells.add(nb_cell)
            return new_candidate_cells

        def compute_potential_old(cell):
            """
            Computes the potential in one cell, using potential in neighbouring cells.
            """
            # Find the minimal directions along a grid cell.
            # Assume left and below are best, then overwrite with right and up if they are better
            neighbour_pots = {direction: np.Inf for direction in ft.DIRECTIONS}

            hor_potential = ver_potential = 0
            hor_cost = ver_cost = np.Inf

            for direction in ft.DIRECTIONS:
                normal = ft.DIRECTIONS[direction]
                # numerical direction
                nb_cell = (cell[0] + normal[0], cell[1] + normal[1])
                if not self._exists(nb_cell):
                    continue
                pot = potential_field[nb_cell]
                # potential in that neighbour field
                if direction == 'right':
                    face_index = (nb_cell[0] - 1, nb_cell[1])
                elif direction == 'up':
                    face_index = (nb_cell[0], nb_cell[1] - 1)
                    # Unit cost values are defined w.r.t faces, not cells!
                else:
                    face_index = nb_cell
                cost = self.unit_field_dict[opposites[direction]].array[face_index]
                # Cost to go from there to here
                neighbour_pots[direction] = pot + cost
                # total potential
                if neighbour_pots[direction] < neighbour_pots[opposites[direction]]:
                    if direction in ft.HORIZONTAL_DIRECTIONS:
                        hor_potential = pot
                        hor_cost = cost
                        # lowest in horizontal direction
                    elif direction in ft.VERTICAL_DIRECTIONS:
                        ver_potential = pot
                        ver_cost = cost
                        # lowest in vertical direction
                    else:
                        raise ValueError("Direction unknown")
            # Coefficients of quadratic equation
            a = 1 / hor_cost ** 2 + 1 / ver_cost ** 2
            b = -2 * (hor_potential / hor_cost ** 2 + ver_potential / ver_cost ** 2)
            c = (hor_potential / hor_cost) ** 2 + (ver_potential / ver_cost) ** 2 - 1

            D = b ** 2 - 4 * a * c
            x_high = (2 * c) / (-b - math.sqrt(D))
            # Might not be obvious, but why we take the largest root is found in report.
            return x_high

        candidate_cells = {cell: compute_potential(cell[0],cell[1], self.grid_dimension[0],self.grid_dimension[1], potential_field,
                                                   self.unit_field_dict['left'].array,self.unit_field_dict['right'].array,
                                                   self.unit_field_dict['up'].array,self.unit_field_dict['down'].array, 9999)
                           for cell in get_new_candidate_cells(known_cells)}

        new_candidate_cells = get_new_candidate_cells(known_cells)
        # Todo: Proposed improvement: new_candidate_cells = candidate_cells.keys()
        while unknown_cells:
            for candidate_cell in new_candidate_cells:
                if False:
                    potential = compute_potential(candidate_cell)
                else:
                    time1 = time.time()
                    potential = compute_potential(candidate_cell[0],candidate_cell[1], self.grid_dimension[0],self.grid_dimension[1], potential_field,
                                                  self.unit_field_dict['left'].array,self.unit_field_dict['right'].array,
                                                  self.unit_field_dict['up'].array,self.unit_field_dict['down'].array, 9999)
                    time2 = time.time()
                    print("time: %.2e"%(time2-time1))
                candidate_cells[candidate_cell] = potential
            sorted_candidates = sorted(candidate_cells.items(), key=operator.itemgetter(1))  # Todo: Can we reuse this?
            best_cell = sorted_candidates[0][0]
            min_potential = candidate_cells.pop(best_cell)
            potential_field[best_cell] = min_potential
            unknown_cells.remove(best_cell)
            known_cells.add(best_cell)
            new_candidate_cells = get_new_candidate_cells({best_cell})
        self.potential_field.update(potential_field)
