import argparse
from hp_model import *
from genetic_algorithm import *


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "seq",
        type=str,
        help="""Input protein sequence. Either provide a sequence or 
                choose from a range samples in sample_seqs.txt""",
    )
    parser.add_argument(
        "--pop_size",
        type=int,
        default=1000,
        help="Population size of each generation of the evolutionary algorithm",
    )
    parser.add_argument(
        "--num_gens",
        type=int,
        default=500,
        help="Number of generations as the stopping criterion of the iterative algorithm",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=1e-5,
        help="Tolerance value as the stopping criterion of the iterative algorithm",
    )
    parser.add_argument(
        "--energy_func",
        type=str,
        choices=["berger", "custodio"],
        default="custodio",
        help="Choice of energy function in the HP model that corresponds to fitness function in the genetic algorithm",
    )
    parser.add_argument(
        "--log",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Log the iterative optimization process or not",
    )
    parser.add_argument(
        "--display",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Display the converged solution on a mesh grid or not",
    )
    args = parser.parse_args()

    try:
        seq_choice = eval(args.seq) - 1
        with open("sample_seqs.txt", "r") as file:
            seqs = [row for row in file.read().splitlines()]
            seq = seqs[seq_choice]
    except:
        seq = args.seq
    hp_seq = hp_format(seq)

    energy_func_map = {"berger": energy_function_1, "custodio": energy_function_2}

    converged_solution = genetic_algorithm(
        hp_seq,
        args.pop_size,
        args.num_gens,
        random_relative_fold,
        energy_func_map[args.energy_func],
        args.epsilon,
        args.log,
    )
    with open("optimal_folding.txt", "w") as file:
        file.write("{}\n{}".format(seq, converged_solution))

    if args.display:
        path = []
        lattice = generate_lattice(converged_solution, hp_seq, path)
        draw_lattice(lattice, path, energy_func_map[args.energy_func](converged_solution, hp_seq))
