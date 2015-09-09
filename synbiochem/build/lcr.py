# import collections
import itertools


def get_dominos(melting_temp, sequences, plasmid_seq=None, shuffle=False):
    num_sequences = len(sequences)
    orderings = sorted([list(a) for a in set(itertools.permutations(range(num_sequences)))]) \
        if shuffle \
        else [range(num_sequences)]

    if plasmid_seq is not None:
        sequences.append(plasmid_seq)
        for ordering in orderings:
            ordering.insert(0, num_sequences) 
            ordering.append(num_sequences)

    print orderings

    pairs = []

    for ordering in orderings:
        pairs.extend(pairwise(ordering))

    print pairs
    # print collections.Counter(pairs)
    print [a for a in sorted(set(pairs))]
    
    for pair in sorted(set(pairs)):
        print pair[0], pair[1], sequences[pair[0]][-3:], sequences[pair[1]][:3], sequences[pair[0]][-3:] + sequences[pair[1]][:3]


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)


get_dominos(60, ['saaaaae', 'sccccce', 'sggggge'], 'sttttte', False)