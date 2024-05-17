import random
import matplotlib.pyplot as plt


def hp_format(prot_seq):
    """Return the HP representation of a protein sequence"""
    return "".join(
        ["H" if aa in ["A", "V", "I", "L", "M", "F", "Y", "W"] else "P" for aa in prot_seq]
    )


def random_relative_fold(hp_seq):
    """Return a random fold using relative directions on a 2D lattice"""
    fold = "".join(["-"] + random.choices(["L", "R", "F"], k=len(hp_seq) - 1))
    # Generate only valid folds, i.e., folds that avoid walking on themselves
    while not check_saw(fold, hp_seq):
        fold = "".join(["-"] + random.choices(["L", "R", "F"], k=len(hp_seq) - 1))
    return fold


def generate_lattice(fold, hp_seq, path=None):
    """Return a lattice corresponding to a folding of a sequence"""
    # Maximum axis size needed to represent a fold is twice the fold \
    # length (not counting the origin point)
    axis_size = (len(fold) * 2) - 1

    # A 2D matrix to represent the 2D lattice
    # e.g. lattice[i][j] = [num_AAs_in_cell, cell_status from {Empty, H, P}]
    lattice = [[[0, "E"] for j in range(axis_size)] for i in range(axis_size)]

    # Euclidean space directions decoded into relative directions
    # e.g. Up: [Left: (add_to_x, add_to_y, next_dir)
    #          Right: (add_to_x, add_to_y, next_dir)
    #          Forward: (add_to_x, add_to_y, next_dir)]
    directions = {
        "U": [(0, -1, "L"), (0, 1, "R"), (-1, 0, "U")],
        "D": [(0, 1, "R"), (0, -1, "L"), (1, 0, "D")],
        "R": [(-1, 0, "U"), (1, 0, "D"), (0, 1, "R")],
        "L": [(1, 0, "D"), (-1, 0, "U"), (0, -1, "L")],
    }

    # Start at the middle of the lattice or the origin point
    current_position = (axis_size // 2, axis_size // 2)
    current_direction = "R"

    # Put the first amino acid in its cell according to the folding provided
    # Then update cell status and add its coordinates to path
    lattice[current_position[0]][current_position[1]][0] += 1
    lattice[current_position[0]][current_position[1]][1] = hp_seq[0]
    if path != None:
        path.append((current_position[0], current_position[1]))

    # Do the same thing for the rest of the amino acids
    for hp in range(1, len(fold)):
        if fold[hp] == "L":
            current_position = (
                current_position[0] + directions[current_direction][0][0],
                current_position[1] + directions[current_direction][0][1],
            )
            current_direction = directions[current_direction][0][2]
            lattice[current_position[0]][current_position[1]][0] += 1
            lattice[current_position[0]][current_position[1]][1] = hp_seq[hp]

        elif fold[hp] == "R":
            current_position = (
                current_position[0] + directions[current_direction][1][0],
                current_position[1] + directions[current_direction][1][1],
            )
            current_direction = directions[current_direction][1][2]
            lattice[current_position[0]][current_position[1]][0] += 1
            lattice[current_position[0]][current_position[1]][1] = hp_seq[hp]

        else:
            current_position = (
                current_position[0] + directions[current_direction][2][0],
                current_position[1] + directions[current_direction][2][1],
            )
            current_direction = directions[current_direction][2][2]
            lattice[current_position[0]][current_position[1]][0] += 1
            lattice[current_position[0]][current_position[1]][1] = hp_seq[hp]

        if path:
            path.append((current_position[0], current_position[1]))
    return lattice


def draw_lattice(lattice, path, fold_energy=None):
    """Display a lattice on a mesh grid"""
    # Make the path suitable to be displayed
    dif_view_path = [(cell[1] - len(lattice), len(lattice) - cell[0]) for cell in path]

    # Define the domain of each axis and the colorcodes of their cells
    x_linespace, y_linespace = [], []
    colorcodes = [[0 for j in range(len(lattice))] for i in range(len(lattice))]

    # Fill the above variables
    for i in range(len(lattice)):
        for j in range(len(lattice)):
            y_linespace.append(len(lattice) - i)
            x_linespace.append(j - len(lattice))
            if not lattice[i][j][0]:
                colorcodes[i][j] = "20"
            elif lattice[i][j][0] == 1 and lattice[i][j][1] == "H":
                colorcodes[i][j] = "5"
            elif lattice[i][j][0] == 1 and lattice[i][j][1] == "P":
                colorcodes[i][j] = "1"
            else:
                colorcodes[i][j] = "10"

    # Draw the lattice
    plt.figure(figsize=(10, 10))
    plt.tick_params(axis="x", bottom=False, labelbottom=False)
    plt.tick_params(axis="y", left=False, labelleft=False)
    plt.scatter(x_linespace, y_linespace, c=colorcodes, cmap="Dark2")
    plt.plot([cell[0] for cell in dif_view_path], [cell[1] for cell in dif_view_path])
    if not fold_energy == None:
        plt.title("Energy = " + str(fold_energy))
    plt.show()


def check_saw(fold, hp_seq=None):
    """Function to check if a fold is avoiding self-walk"""
    if not hp_seq:
        hp_seq = "".join(["H" for hp in fold])

    # Decode the fold and then check for SAW
    lattice = generate_lattice(fold, hp_seq)
    for row in lattice:
        for cell in row:
            if cell[0] > 1:
                return False
    return True


def energy_function_1(fold, hp_seq):
    """Function to act as the Berger energy function"""
    path = []
    lattice = generate_lattice(fold, hp_seq, path)
    hh_contacts = 0

    # Go along the path and count the valid contacts
    for ind in range(len(path)):
        i, j = path[ind]
        if lattice[i][j][1] == "H":
            valid_surrounding_cells = [
                cell
                for cell in [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]
                if (cell[0] > -1 and cell[1] > -1)
                and (cell[0] < len(lattice) and cell[1] < len(lattice))
            ]

            # Connected cells do not count as contacts
            for cell in valid_surrounding_cells:
                if ind - 1 >= 0 and ind + 1 < len(path):
                    if path[ind - 1] in valid_surrounding_cells:
                        valid_surrounding_cells.remove(path[ind - 1])
                    if path[ind + 1] in valid_surrounding_cells:
                        valid_surrounding_cells.remove(path[ind + 1])
                elif ind - 1 < 0:
                    if path[ind + 1] in valid_surrounding_cells:
                        valid_surrounding_cells.remove(path[ind + 1])
                else:
                    if path[ind - 1] in valid_surrounding_cells:
                        valid_surrounding_cells.remove(path[ind - 1])

            hh_contact = sum(
                [1 for cell in valid_surrounding_cells if lattice[cell[0]][cell[1]][1] == "H"]
            )
            hh_contacts += hh_contact
    return -1 * (hh_contacts / 2)


def energy_function_2(fold, hp_seq, w1=0, w2=10, w3=40):
    """Function to act as the Custodio energy function"""
    path = []
    lattice = generate_lattice(fold, hp_seq, path)
    hh_contacts = 0
    hp_contacts = 0
    hs_contacts = 0

    # Go along the path and count the valid contacts
    for ind in range(len(path)):
        i, j = path[ind]
        if lattice[i][j][1] == "H":
            valid_surrounding_cells = [
                cell
                for cell in [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]
                if (cell[0] > -1 and cell[1] > -1)
                and (cell[0] < len(lattice) and cell[1] < len(lattice))
            ]

            # Connected cells do not count as contacts
            for cell in valid_surrounding_cells:
                if ind - 1 >= 0 and ind + 1 < len(path):
                    if path[ind - 1] in valid_surrounding_cells:
                        valid_surrounding_cells.remove(path[ind - 1])
                    if path[ind + 1] in valid_surrounding_cells:
                        valid_surrounding_cells.remove(path[ind + 1])
                elif ind - 1 < 0:
                    if path[ind + 1] in valid_surrounding_cells:
                        valid_surrounding_cells.remove(path[ind + 1])
                else:
                    if path[ind - 1] in valid_surrounding_cells:
                        valid_surrounding_cells.remove(path[ind - 1])

            hh_contact = sum(
                [1 for cell in valid_surrounding_cells if lattice[cell[0]][cell[1]][1] == "H"]
            )
            hh_contacts += hh_contact

            hp_contact = sum(
                [1 for cell in valid_surrounding_cells if lattice[cell[0]][cell[1]][1] == "P"]
            )
            hp_contacts += hp_contact

            hs_contact = sum(
                [1 for cell in valid_surrounding_cells if lattice[cell[0]][cell[1]][1] == "E"]
            )
            hs_contacts += hs_contact
    return (w1 * hh_contacts / 2) + (w2 * hp_contacts) + (w3 * hs_contacts)
