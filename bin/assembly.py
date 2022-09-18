import re
from collections import defaultdict


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


def overlap_graph(overlaplist):
    """
    prints out sets of kmers where the first k-1
    basepairs overlap with the last k-1 basepairs of another
    :param overlaplist: list of kmers of the same length
    :return: prints overlap of each kmer separated by ->
    """
    overlapdict = defaultdict(list)
    for i in range(len(overlaplist)):
        for j in range(len(overlaplist)):
            if i == j: continue
            if overlaplist[i][1:] != overlaplist[j][:-1]:
                continue
            overlapdict[overlaplist[i]].append(overlaplist[j])
    for k, v in overlapdict.items():
        print(f"{k} -> {','.join(v)}")


def binary_kmers(repeat=4):
    pools = [tuple('01')] * repeat
    result = [[]]
    for pool in pools:
        result = [x + [y] for x in result for y in pool]
    resultstrings = ["".join(x) for x in result]
    return resultstrings

def debruijn_graph(dnastring, k):
    """
    similar to overlap graph but the nodes are length k-1
    and the kmers are the edges
    :param dnastring: string of DNA
    :return: prints debruijn graph
    """
    dbdict = defaultdict(list)
    for i in range(len(dnastring)-k+1):
        kmer = dnastring[i:i+k]
        dbdict[kmer[:-1]].append(kmer[1:])
    for k, v in dbdict.items():
        print(f"{k} -> {','.join(v)}")
