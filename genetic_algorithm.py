import numpy as np
from hp_model import check_saw


def genetic_algorithm(_input, pop_size, num_gens, rand_generator, fitness_func, epsilon, log):
    """The Core function of the genetic algorithm"""
    # Create an initial population
    population = [rand_generator(_input) for _ in range(pop_size)]

    # Keep track of convergence
    last_average_score = sum([fitness_func(chrom, _input) for chrom in population]) / len(
        population
    )
    if log:
        print("Gen 1 average score: ", last_average_score)

    # Run the main loop
    for generation in range(num_gens):
        next_population_candidates = []
        while len(set(next_population_candidates)) < pop_size:
            children = crossover(_input, population, fitness_func)
            mutated_children = [mutation(_input, c) for c in children]
            next_population_candidates += mutated_children
        next_population_candidates += population
        population = selection_for_survival(
            next_population_candidates, _input, pop_size, fitness_func
        )

        average_score = sum([fitness_func(chrom, _input) for chrom in population]) / len(population)
        if log and (generation + 1) % 5 == 0:
            print("Gen {} average score: {}".format(generation + 1, last_average_score))

        # Check for convergence
        if abs(last_average_score - average_score) < epsilon:
            if log:
                print("Gen {} average score: {}".format(generation + 1, last_average_score))
            break
        else:
            last_average_score = average_score

    return population[np.argmin([fitness_func(chroms, _input) for chroms in population])]


def selection_for_variation(_input, population, fitness_func, method="tournament", K=2):
    """Function to perform K-tournament selection with a chance for duplicates"""
    if method == "tournament":
        participants = [np.random.choice(population) for _ in range(K)]
        return participants[
            np.argmin([fitness_func(participant, _input) for participant in participants])
        ]


def crossover(_input, population, fitness_func, k=2):
    """Function to perform 1/2-point crossover"""
    # Solutions can also be generated by cloning an existing solution
    # which is analogous to asexual reproduction
    parent1 = selection_for_variation(_input, population, fitness_func)
    parent2 = selection_for_variation(_input, population, fitness_func)

    # Choose the crossover points at random
    crossover_points = sorted(np.random.choice(len(parent1) + 1, k))

    # Generate the children
    if k == 1:
        child1 = parent1[: crossover_points[0]] + parent2[crossover_points[0] :]
        child2 = parent2[: crossover_points[0]] + parent1[crossover_points[0] :]
    elif k == 2:
        child1 = (
            parent1[: crossover_points[0]]
            + parent2[crossover_points[0] : crossover_points[1]]
            + parent1[crossover_points[1] :]
        )
        child2 = (
            parent2[: crossover_points[0]]
            + parent1[crossover_points[0] : crossover_points[1]]
            + parent2[crossover_points[1] :]
        )

    # If invalid children are generated, then parents are copied unchanged
    min_energy_parent = [parent1, parent2][
        np.argmin([fitness_func(p, _input) for p in [parent1, parent2]])
    ]
    if not check_saw(child1, _input) and not check_saw(child2, _input):
        return parent1, parent2
    elif not check_saw(child1, _input) and check_saw(child2, _input):
        return child2, min_energy_parent
    elif check_saw(child1, _input) and not check_saw(child2, _input):
        return child1, min_energy_parent
    else:
        return [child1, child2]


def mutation(_input, chromosome):
    """Function to perform the single-point mutation"""
    # Uniform mutation is applied to all individuals in the children population
    # Each encoding position is randomly mutated with a probability of 1/(Length-1)
    encoding_positions = [position for position in chromosome[1:]]
    mutated_chromosome = []
    mutation_rate = 1 / (len(encoding_positions) - 1)
    for position in encoding_positions:
        if np.random.random() < mutation_rate:
            possible_mutations = ["L", "R", "F"]
            possible_mutations.remove(position)
            mutated_chromosome.append(np.random.choice(possible_mutations))
        else:
            mutated_chromosome.append(position)

    # If mutation of a position leads to an invalid conformation
    # the original value is restored
    mutated_chromosome = "-" + "".join(mutated_chromosome)
    if check_saw(mutated_chromosome, _input):
        return mutated_chromosome
    else:
        return chromosome


def selection_for_survival(next_population_candidates, _input, pop_size, fitness_func):
    """Function to perform selection for survival to the next evolution"""
    # The best individuals from the union of parents and children are selected
    next_population_candidates = list(set(next_population_candidates))
    next_population = sorted(next_population_candidates, key=lambda x: fitness_func(x, _input))
    return next_population[:pop_size]
