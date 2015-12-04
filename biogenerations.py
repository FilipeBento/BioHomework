#!/usr/bin/python

from __future__ import print_function
import random
import itertools
import math
import matplotlib.pyplot as plt


def should_mutate(mutation_rate):
    return random.random() <= mutation_rate


def should_recombine(recombination_rate):
    return random.random() <= recombination_rate


def mutation_random_pos(sequence):
    return random.randrange(0, len(sequence))


def recombination_random_pos(sequence, rfl):
    return random.randrange(0, len(sequence) - rfl)


def mutate_sequence(sequence):
    """
    Replaces a random char with a value (A, C, T or G) other than its own.
    :param sequence: the string sequence to be altered.
    :return: An altered copy of the original seq.
    """
    pos_to_mutate = mutation_random_pos(sequence)
    nucleotides = ['A', 'C', 'T', 'G']
    nucleotides.remove(sequence[pos_to_mutate])
    mutated_seq = list(sequence)
    mutated_seq[pos_to_mutate] = random.choice(nucleotides)
    return ''.join(mutated_seq)


def recombine_sequence(first_seq, second_seq, rfl):
    """
     Receives a string sequence, a second string sequence and the rfl and alters the string sequence using a
     random subset of size 'rfl' from the second string.
     string X string X int --> string

    :param first_seq: string sequence to be recombined.
    :param second_seq: the second string from which a subsequence will be extracted for recombination in the first seq.
    :param rfl: the recombination fragment length (positive integer)
    :return: The first sequence recombined
    """
    assert rfl > 0
    assert len(first_seq) > 0
    assert len(second_seq) > 0
    if len(second_seq) < rfl:
        raise ValueError("rfl exceeds length of sequence.")
    elif len(second_seq) == rfl:
        return second_seq
    start_pos = recombination_random_pos(first_seq, rfl)
    end = start_pos + rfl
    assert end <= len(first_seq)
    section = second_seq[start_pos:end]
    return replace_interval(first_seq, section, start_pos)


def replace_interval(string, replacement, start):
    """
    Replaces a subset of a 'string' with a certain 'replacement' string, starting at position 'start'.
    string X string X int --> string
    :return: the altered copy of the received string
    """
    assert 0 <= start < len(string)
    replacement_size = len(replacement)
    if start + replacement_size <= len(string):
        return string[:start] + replacement + string[start + replacement_size:]
    else:
        raise ValueError("replacement exceeds string size.")


# Receives a file object representing a FASTA file and parses its contents.
# Returns a dictionary of sequence names containing the corresponding sequences.
# EX: {'seq1': 'AAATTTCCC', 'seq2': 'TTTTGGGG', etc}
def parse_file(file):
    lines = file.readlines()
    sequences = {}
    current_sequence_name = ''
    current_sequence = []
    for line in lines:
        if line.isspace() or line.startswith(';'):
            continue
        if line.startswith('>'):
            # save previous sequence, if there is one
            if current_sequence:
                sequences[current_sequence_name] = current_sequence
            # trim '>' and '\n'
            current_sequence_name = line[1:-1]
            current_sequence = []
            continue
        # store the sequence without the '\n'
        current_sequence.extend(list(line[:-1]))
    if current_sequence_name not in sequences.keys():
        sequences[current_sequence_name] = current_sequence
    return sequences


def create_output_file(sequences, file_name):
    with open(file_name, 'w') as out:
        for seq_name, seq in sequences.iteritems():
            print('>' + seq_name, file=out)
            print(str(seq), file=out)
            print('\n\n', file=out)


def generations(sequences, mutation_rate, recombination_rate, rfl, num_generations):
    """
    Simulates the generations of the sequences, with the values being stored in a dictionary mapping the generation
    number to its sequences.

    Simple example:
    Output (assuming 2 seqs, 2 generations, no recomb and always mutate)
    {
        0: { <--- initial sequences
            'seq1': 'AAAAAAAAA',
            'seq2': 'TTTTTTTTT',
        },
        1: {
            'seq1': 'AAAAAAAGA',
            'seq2': 'TTTTATTTT',
            },
        2: {
            'seq1': 'AACAAAAGA',
            'seq2': 'TTTTATGTT',
        }
    }
    :param sequences: sequences to be recombined
    :param mutation_rate: probability of mutation
    :param recombination_rate: probability of recombination
    :param rfl: fragment length for recombination
    :param num_generations: the number of generations
    :return:
    """
    gens = {}
    prev_generation = dict(sequences)
    for i in xrange(num_generations):
        for seq_name, seq in prev_generation.iteritems():
            if should_mutate(mutation_rate):
                prev_generation[seq_name] = mutate_sequence(seq)
        for seq_name, seq in prev_generation.iteritems():
            if should_recombine(recombination_rate):
                previous_sequences = prev_generation.values()
                previous_sequences.remove(seq)
                random_seq = random.choice(previous_sequences)
                prev_generation[seq_name] = recombine_sequence(seq, random_seq, rfl)
        gens[i] = dict(prev_generation)
    return gens


def create_plot(hamming_values, jukes_values):
    plt.subplot(2, 2, 1)
    plt.plot(list(range(0, len(hamming_values))), hamming_values, marker='o', linestyle='--', color='r')
    plt.ylim([0.0, 1.1])
    plt.subplot(2, 2, 2)
    plt.plot(list(range(0, len(jukes_values))), jukes_values, marker='x', linestyle='--', color='g')
    plt.xlabel('Num Generations')
    plt.ylabel('Distance (Hamming and Jukes)')
    plt.title('Question 1')
    plt.show()


def question1():
    # Generate random seq of size 100 and copy it 100 times
    sequences = {}
    population_size = 100
    sequence_size = 100
    sequence = ''.join(random.choice(('A', 'C', 'T', 'G')) for _ in xrange(sequence_size))
    for i in xrange(population_size):
        sequences["seq%d" % i] = sequence
    # Run the simulation
    gens = generations(sequences, mutation_rate=0.01, recombination_rate=0.01, rfl=5, num_generations=5000)
    hamming_points = []
    for gen_number, gen in gens.iteritems():
        if gen_number == 0:
            continue
        combinations = list(itertools.combinations(gen.itervalues(), 2))
        total_distance = 0.0
        total_combinations = 0
        for s1, s2 in combinations:
            total_distance += float(hamming_distance(s1, s2))
            total_combinations += 1
        hamming_points.append(float(total_distance / total_combinations))
    print('hamming points: ', hamming_points)
    jukes_points = map(jukes_cantor_model, hamming_points)
    print('jukes points: ', jukes_points)
    create_plot(hamming_points, jukes_points)


def jukes_cantor_model(point):
    if point >= 0.75:
        return point
    return -3.0/4 * math.log(1 - 4/float(3) * point)


def hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    mismatches = 0
    for a, b in zip(s1, s2):
        if a != b:
            mismatches += 1
    return mismatches / float(len(s1))


# QUESTION 2 - RUN 5 SIMULATIONS OF RANDOM SEQS AND OUTPUT THE RESULTS TO FASTA FILES
def question2():
    letters = ('A', 'C', 'T', 'G')
    seq_size = 100
    seq1 = ''.join(random.choice(letters) for _ in xrange(seq_size))
    seq2 = ''.join(random.choice(letters) for _ in xrange(seq_size))
    sequences = {
        'seq1_1': seq1,
        'seq1_2': seq1,
        'seq1_3': seq1,

        'seq2_1': seq2,
        'seq2_2': seq2,
        'seq2_3': seq2,
    }
    # Run the five simulations and
    results = [
        generations(sequences, mutation_rate=0.1,   recombination_rate=0,     rfl=5, num_generations=1000),
        generations(sequences, mutation_rate=0.1,   recombination_rate=0.1,   rfl=5, num_generations=1000),
        generations(sequences, mutation_rate=0.1,   recombination_rate=0.001, rfl=5, num_generations=1000),
        generations(sequences, mutation_rate=0.001, recombination_rate=0.1,   rfl=5, num_generations=1000),
        generations(sequences, mutation_rate=0.1,   recombination_rate=0.01,  rfl=5, num_generations=2500)
    ]
    # Output the last gen of each simulation to individual files, adding the initial sequences
    for i, result in enumerate(results, start=1):
        # find last gen
        keys = result.keys()
        last_gen = max(keys)
        # adds the initial seqs so the neighbour-joining tree can be created
        result[last_gen]['seq1_initial'] = seq1
        result[last_gen]['seq2_initial'] = seq2
        # fasta file with the last generation
        create_output_file(result[last_gen], "sim%d.fasta" % i)


def main():
    print('COMPUTATIONAL BIOLOGY - Mutations over generations homework')
    print('Running...')
    question2()
    print('Done. Exiting')

if __name__ == '__main__':
    main()
