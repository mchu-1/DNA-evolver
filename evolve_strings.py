# evolve_strings.py

"""Evolutionary algorithm for finding pairs of non-homologous DNA strings."""

from random import choice, randint, random
from statistics import mode


def get_substrings(full_string: str, k: int) -> list[str]:
    """
    Get all substrings of a string with length k.
    Returns a list in left-to-right order.
    """
    substrings = []
    for i in range(len(full_string)-k+1):
        substrings.append(full_string[i:i+k])

    return substrings


def get_sliding_window(substring: str) -> int:
    """
    Get size of sliding window for a substring search.
    """
    window = len(substring)//2-1 if len(substring) > 1 else 0

    return window


def score_substrings(substrings_1, substrings_2: list[str], sliding_window: int) -> int:
    """
    Score two lists of substrings by number of matches within a sliding window.
    Input two equal length lists of substrings with equal length.
    Output homology score between the lists.
    """
    homology_score = 0
    list_length = len(substrings_1)

    for i in range(list_length):
        left_index = i - sliding_window if i >= sliding_window else 0
        right_index = i + sliding_window + 1 if list_length - i > sliding_window else list_length
        if substrings_1[i] in substrings_2[left_index: right_index]:
            homology_score += 1

    return homology_score


def score_fitness(string_1, string_2: str) -> float:
    """
    Get transposition distance between two strings.
    Based on number of common substrings within size-based sliding windows.
    """
    if len(string_1) != len(string_2):
        raise ValueError("Strings should be equal in length.")

    length = len(string_1)
    score = 0
    max_score = 0

    for k in range(1, length):
        substrings_1 = get_substrings(string_1, k)
        substrings_2 = get_substrings(string_2, k)
        sliding_window = get_sliding_window("." * k)
        homology_score = score_substrings(substrings_1, substrings_2, sliding_window)
        score += homology_score
        # sum homology scores across all substrings
        identity = ["." * k]*(length - k + 1)
        max_score += score_substrings(identity, identity, sliding_window)
        # maximum distance score

    fitness = score/max_score

    return fitness


def generate_string(length: int) -> str:
    """
    Generate random string of some given length.
    Select for strings without RNA polymerase III terminators.
    """
    new_bases = []
    count = 0 # count number of consecutive T's
    for x in range(length):
        new_base = choice(bases)
        if new_base == "T":
            count += 1
        if count == 4:
            new_base = "C"
            count = 0 # reset count
        new_bases.append(new_base)

    new_string = "".join(new_bases)

    return new_string


def generate_population(string_length, size: int) -> list[tuple]:
    """
    Generate population of n pairs of random strings of a given length.
    """
    population = []
    for i in range(size):
        new_strings = []
        for x in range(2):
            new_strings.append(generate_string(string_length)) # generate a new pair of strings
        population.append(tuple(new_strings))

    return population


def mutate_base(base: str) -> str:
    """
    Change one base to a different base.
    """
    new_base = choice([b for b in bases if b != base])

    return new_base


def mutate_string(original_string: str, rate: float) -> str:
    """
    Mutate a string at a given rate (per base).
    """
    if rate < 0 or rate > 1:
        raise ValueError("Mutation rate must be in the range of 0 to 1.")

    new_bases = []
    count = 0 # count T's
    for base in original_string:
        new_base = base
        if random() < rate:
            new_base = mutate_base(base)
        if new_base == "T":
            count += 1
        if count == 4:
            new_base = "C"
            count = 0 # reset count
        new_bases.append(new_base)

    new_string = "".join(new_bases)

    return new_string


def remove_terminators(sequence: str) -> str:
    """
    Remove terminator from a given sequence.
    """
    count = 0
    new_bases = []
    for base in sequence:
        if base == "T":
            count += 1
        if count == 4:
            new_bases.append("A") # flip fourth T
            count = 0 # reset count
        else:
            new_bases.append(base)

    new_sequence = "".join(new_bases)

    return new_sequence


def cross_strings(strings_1, strings_2: tuple[str]) -> list[tuple]:
    """
    Cross two pairs of parent strings to generate two new pairs of child strings.
    """
    child_strings_1 = (strings_1[0], strings_2[1])
    child_strings_2 = (strings_1[1], strings_2[0])

    children = [child_strings_1, child_strings_2]

    return children


def generate_new_strings(population: list[tuple], tournament_size: int, mutation_rate: float) -> list[tuple]:
    """
    Generate new string pairs from a given population of string pairs.
    Create new children based on tournaments of a given size between existing strings.
    """
    fitness = []
    for strings in population:
        new_score = score_fitness(strings[0], strings[1])
        fitness.append(new_score) # build fitness landscape for strings

    suitors = [] # suitors for each string in current population based on fitness
    for i, strings in enumerate(population):
        best_suitor = strings
        for j in range(tournament_size):
            suitor = randint(0, len(population)-1) # index of potential suitor strings
            if fitness[suitor] < fitness[i]:
                best_suitor = population[suitor]
        suitors.append(best_suitor)

    selection = []
    for x in range(len(population)//2):
        new_selection = randint(0, len(population)-1)
        selection.append(new_selection) # select indexes of strings to mate

    new_generation = []
    for i in selection:
        children = cross_strings(population[i], suitors[i])
        mutated_children = []
        for child in children:
            mutant_child = tuple(mutate_string(c, mutation_rate) for c in child)
            mutated_children.append(mutant_child) # mutate children at specified rate
        new_generation.extend(mutated_children)

    return new_generation


def evolve_strings(population_size, string_length, generations, tournament_size:  int, mutation_rate: float) -> list[tuple]:
    """
    Evolve pairs of strings over a given number of generations and tournament size.
    """
    if tournament_size >= population_size:
        raise ValueError("Tournament size must be less than population size.")

    print(f"Initiating population of size {population_size} with {string_length} bp strings ...")
    print(f"Tournament size at {tournament_size} and mutation rate at {mutation_rate} per base ...")
    population = generate_population(string_length, population_size)

    count = 0
    for i in range(generations):
        print(f"Evolving generation {i + 1} ...")
        if len(population) == 0:
            break
        elif tournament_size >= len(population):
            break
        else:
            population = generate_new_strings(population, tournament_size, mutation_rate)
            count += 1

    fitness = [] # final fitness landscape
    for strings in population:
        score = score_fitness(strings[0], strings[1])
        fitness.append(score)

    evolved_strings = []
    for i, strings in enumerate(population):
        evolved_strings.append((fitness[i], strings[0], strings[1]))
        # evolved strings with fitness scores

    evolved_strings.sort(key = lambda s: s[0]) # sort strings based on fitness

    print(f"Final iteration at generation {count}.")

    return evolved_strings


if __name__ == "__main__":
    bases = ["A", "C", "T", "G"]  # DNA bases
    print("###########################################")
    print("Evolve pairs of non-homologous DNA strings.")
    print("###########################################")
    str_length = 60
    evolved_population = evolve_strings(population_size = 1000, string_length = str_length, generations = 1000, tournament_size = 10, mutation_rate = 1/(5*str_length))
    print("Final population:")
    for ix, string in enumerate(evolved_population):
        print(f"{ix+1} - {string}")
    print(f"Mode fitness = {mode(string[0] for string in evolved_population)}") # report modal fitness
    print("Evolution complete!")












