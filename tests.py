import random
import biogenerations as bio


def test_generations():
    letters = ('A', 'C', 'T', 'G')
    seq_size = 100
    seq1 = ''.join(random.choice(letters) for _ in xrange(seq_size))
    seq2 = ''.join(random.choice(letters) for _ in xrange(seq_size))
    print 'seq1 is ', seq1
    print 'seq2 is ', seq2
    sequences = {
        'seq1_1': seq1,
        'seq1_2': seq1,
        'seq1_3': seq1,

        'seq2_1': seq2,
        'seq2_2': seq2,
        'seq2_3': seq2,
    }
    result = bio.generations(sequences, 1, 1, 4, 4)
    print result


def mutate_sequences():
    for i in xrange(1000):
        seq = ''.join(random.choice(('A', 'C', 'T', 'G')) for _ in xrange(random.randint(1, 100)))
        assert bio.mutate_sequence(seq) != seq
        assert len(bio.mutate_sequence(seq)) == len(seq)


def recombine_sequences():
    print 'RECOMB --------------------------'
    for i in xrange(1000):
        rfl = 3
        size = random.randint(3, 20)
        seq = ''.join(random.choice(('A', 'C', 'T', 'G')) for _ in xrange(size))
        seq2 = ''.join(str(i) for i in xrange(size))
        recombined_seq = bio.recombine_sequence(seq, seq2, rfl)
        print 'original = ', seq
        print 'recomb   = ', recombined_seq
        assert len(recombined_seq) == len(seq) and recombined_seq != seq


def replace_interval(string, replacement, start):
    assert 0 <= start < len(string)
    replacement_size = len(replacement)
    if start + replacement_size <= len(string):
        return string[:start] + replacement + string[start + replacement_size:]
    else:
        raise ValueError("replacement exceeds string size.")


def test_replaces():
    test = "Hello, all, I'm Jamie Cullum. Today I'm going to ruin two pianos."
    test2 = "XXXXX"
    print test
    print replace_interval(test, test2, 0)

test_generations()