import re


def kmer_composition(text, k):
    """
    generate all kmers of a DNA string (including repeats)
    in alphabetical order
    :param text: DNA string
    :param k: length of desired kmers
    :return:
    """
    l = len(text)
    kmerlist = [text[i:i+k] for i in range(l-k+1)]
    sortedlist = sorted(kmerlist)
    return sortedlist


def string_recomposition(kmerlist):
    """
    Given a list of kmers that cascade exactly to form a
    DNA string, returns that string
    :param kmerlist: list of equal length kmers
    :return: string of DNA
    """
    l = len(kmerlist)
    searchlength = len(kmerlist[0])-1
    beginnings = set([x[:-1] for x in kmerlist])
    endings = set([x[1:] for x in kmerlist])
    startsegment = list(beginnings - endings)
    startstring = [x for x in kmerlist if re.match(startsegment[0], x)][0]
    print(startstring)
    for i in range(l-1):
        print(i)
        next = [x for x in kmerlist if re.match(startstring[-searchlength:], x)][0]
        print(next)
        startstring = startstring + next[-1]
        print(startstring)
    return startstring


def pathtogenome(path):
    """
    given a sequential path of tiled kmers,
    reconstruct genome sequence
    :param path: list of kmers in order
    :return: string of DNA
    """
    startingstring = path[0]
    for i in range(1, len(path)):
        startingstring = startingstring + path[i][-1]
    return startingstring


def overlap_graph()