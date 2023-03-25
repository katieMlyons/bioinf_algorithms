import time

from src.ori import better_frequent_kmer, frequency_table
from src.dna import immediate_neighbors, hamming
import numpy as np
import math
import random

# there is a 9-mer in each of these strings
stringlist = ["atgaccgggatactgataaaaaaaagggggggggcgtacacattagataaacgtatgaagtacgttagactcggcgccgccg",
 "acccctattttttgagcagatttagtgacctggaaaaaaaatttgagtacaaaacttttccgaataaaaaaaaaggggggga",
 "tgagtatccctgggatgacttaaaaaaaagggggggtgctctcccgatttttgaatatgtaggatcattcgccagggtccga",
 "gctgagaattggatgaaaaaaaagggggggtccacgcaatcgcgaaccaacgcggacccaaaggcaagaccgataaaggaga",
 "tcccttttgcggtaatgtgccgggaggctggttacgtagggaagccctaacggacttaataaaaaaaagggggggcttatag",
 "gtcaatcatgttcttgtgaatggatttaaaaaaaaggggggggaccgcttggcgcacccaaattcagtgtgggcgagcgcaa",
 "cggttttggcccttgttagaggcccccgtaaaaaaaagggggggcaattatgagagagctaatctatcgcgtgcgtgttcat",
 "aacttgagttaaaaaaaagggggggctggggcacatacaagaggagtcttccttatcagttaatgctgtatgacactatgta",
 "ttggcccattggctaaaagcccaacttgacaaatggaagatagaatccttgcataaaaaaaagggggggaccgaaagggaag",
 "ctggtgagcaacgacagattcttacgtgcattagctcgcttccggggatctaatagcacgaagcttaaaaaaaaggggggga"]

masterstring = "".join(stringlist)

findmotif = better_frequent_kmer(masterstring, 15) # 'aaaaaaaaggggggg'

def motif_enumeration(dna, k, d):
    """
    Based on first string in the list, find k-d motifs of kmers across the other strings
    :param dna: list of strings
    :param k: integer length of k-mer
    :param d: integer maximum number of mismatches
    :return:
    """
    patterns = set()
    for kmer in frequency_table(dna[0], k):
        for kdmotif in immediate_neighbors(kmer, d):
            neighborstotest = immediate_neighbors(kdmotif, d)
            for otherstring in dna[1:]:
                for neighb in neighborstotest:
                    if otherstring.find(neighb) > -1:
                        break
                else:
                    break
            else:
                patterns.add(kdmotif)
    return patterns

# try motif matrix as numpy array
motiflist = [
    ['T', 'C', 'G', 'G', 'G', 'G', 'G', 'T', 'T', 'T', 'T', 'T'],
    ['C', 'C', 'G', 'G', 'T', 'G', 'A', 'C', 'T', 'T', 'A', 'C'],
    ['A', 'C', 'G', 'G', 'G', 'G', 'A', 'T', 'T', 'T', 'T', 'C'],
    ['T', 'T', 'G', 'G', 'G', 'G', 'A', 'C', 'T', 'T', 'T', 'T'],
    ['A', 'A', 'G', 'G', 'G', 'G', 'A', 'C', 'T', 'T', 'C', 'C'],
    ['T', 'T', 'G', 'G', 'G', 'G', 'A', 'C', 'T', 'T', 'C', 'C'],
    ['T', 'C', 'G', 'G', 'G', 'G', 'A', 'T', 'T', 'C', 'A', 'T'],
    ['T', 'C', 'G', 'G', 'G', 'G', 'A', 'T', 'T', 'C', 'C', 'T'],
    ['T', 'A', 'G', 'G', 'G', 'G', 'A', 'A', 'C', 'T', 'A', 'C'],
    ['T', 'C', 'G', 'G', 'G', 'T', 'A', 'T', 'A', 'A', 'C', 'C']
]

motifarray = np.array(motiflist)
byposition = np.transpose(motifarray)
# iterate over transposed array to do columns
profile = np.empty((4,motifarray.shape[1]))
for idx in range(len(byposition)):
    profile[0][idx] = sum(np.char.count(byposition[idx], 'A'))/byposition.shape[1]
    profile[1][idx] = sum(np.char.count(byposition[idx], 'C'))/byposition.shape[1]
    profile[2][idx] = sum(np.char.count(byposition[idx], 'G'))/byposition.shape[1]
    profile[3][idx] = sum(np.char.count(byposition[idx], 'T'))/byposition.shape[1]
# four rows for each nucleotide, 12 cols for each sequence, proportions of each nucleotide
# calculate entropy: -(sum of p*log2 p) where p is freq of each nucleotide
entropies = np.empty((1,profile.shape[1]))
profilebyposition = np.transpose(profile)
for i in range(len(profilebyposition)):
    subsum = 0
    for bp in profilebyposition[i]:
        if bp==0:
            continue
        else:
            subsum += bp*math.log2(bp)
    entropies[0][i] = -1 * subsum


def distancebetweenpatternandstrings(pattern, dna):
    """
    total distance (sites added up) between a pattern and its closest
    kmer on each dna string
    :param pattern: dna string to calculate against
    :param dna: list of lists of sample dna
    :return: minimum distance score between pattern and kmers of dna strings
    """
    k = len(pattern)
    distance = 0
    for dnastring in dna:
        hammingdistance = math.inf
        for substring in frequency_table(dnastring, k):
            if hamming(pattern, substring) < hammingdistance:
                hammingdistance = hamming(pattern, substring)
        distance += hammingdistance
    return distance


def medianstring(dna,k):
    """

    :param dna: list of lists of sample dna
    :param k: length of kmer to search for
    :return: median string with the shortest distance from a kmer in each dna string
    """
    distance = math.inf
    patterns = immediate_neighbors("A"*k, k)  # all possible combos of length k
    for pattern in patterns:
        distbetween = distancebetweenpatternandstrings(pattern, dna)
        if distance > distbetween:
            distance = distbetween
            median = pattern
    return median


def profilemost(text, profile):
    """
    find the most probable kmer in a string
    according to a profile matrix of probabilities
    :param text: dna string
    :param profile: numpy array of 4 rows (A,C,G,T) and k length
    :return: kmer string
    """
    k = profile.shape[1]
    prob = 0
    mostprob = text[:k]
    for kmer in frequency_table(text, k):
        score = 1
        for order, bp in enumerate(kmer):
            for idx, nuc in enumerate(["A", "C", "G", "T"]):
                if bp == nuc:
                    score *= profile[idx][order]
                    break
        if score > prob:
            prob = score
            mostprob = kmer
    return mostprob


def profilemost_random(text, profile):
    """
    sample a kmer according to probabilities
    from a profile
    :param text: dna string
    :param profile: motif matrix with same number of columns
                    as length of desired kmer
    :return: random kmer
    """
    k = profile.shape[1]
    slice = lambda x, y: x[y:y+k]
    problist = []
    for kmer in [slice(text, idx) for idx in range(len(text)-k+1)]:
        score = 1
        for order, bp in enumerate(kmer):
            for idx, nuc in enumerate(["A", "C", "G", "T"]):
                if bp == nuc:
                    score *= profile[idx][order]
                    break
        problist.append(score)
    newidx = gibbs_random(problist)
    return text[newidx:newidx+k]


def motifmatrix(dna, pseudo=True):
    """
    creates motif matrix for set of kmers
    :param dna: list of equal length strings
    :return: numpy array 4 rows by length of strings
    """
    array = np.array([[bp for bp in seq] for seq in dna])
    profile = np.empty((4,array.shape[1]))
    byposition = np.transpose(array)
    for i in range(len(byposition)):
        for idx, nuc in enumerate(("A","C","G","T")):
            profile[idx][i] = sum(np.char.count(byposition[i],nuc))
    if pseudo:
        if 0 in profile:
            profile = (profile + 1) / (byposition.shape[1] + 1)
        else:
            profile = profile / byposition.shape[1]
    else:
        profile = profile / byposition.shape[1]
    return profile


def greedymotifsearch(dna, k):
    """
    greedy algorithm for most probable motif matrix
    :param dna: list of dna strings
    :param k: length of desired motif
    :return: motif matrix
    """
    bestmotifs = [i[0:k] for i in dna] # first kmer in each dna string
    for kmer in frequency_table(dna[0],k):
        motifs = [kmer]
        for i in range(1,len(dna)):
            prof = motifmatrix(motifs)
            motifs.append(profilemost(dna[i], prof))
        if scoreprofile(motifmatrix(motifs), len(dna)) < scoreprofile(motifmatrix(bestmotifs), len(dna)):
            bestmotifs = motifs
    return bestmotifs


def scoreprofile(prof, num):
    score = 0
    for site in np.transpose(prof):
        highest = (1 - max(site)) * num
        score += highest
    return score


def randomizedmotifsearch(dna, k):
    """
    searches for best motif set of length k
    by starting with a random motif set, getting the profile,
    finding the profile-most motifs and repeating
    :param dna: list of strings of DNA
    :param k: length of motifs
    :return: best motif set
    """
    t = len(dna[0])
    slice = lambda x, y: x[y:y+k]
    bestmotifs = [slice(line, random.randint(0, t-k)) for line in dna]
    while True:
        prof = motifmatrix(bestmotifs, pseudo=False)
        bestscore = scoreprofile(prof, len(dna))
        motifs = [profilemost(line, prof) for line in dna]
        if scoreprofile(motifmatrix(motifs, pseudo=False), len(dna)) < bestscore:
            bestmotifs = motifs
        else:
            return bestmotifs, bestscore


def gibbs_random(freqs):
    """
    selects an index in [0,t-1] where t is the length of the frequency list
    with probabilities corresponding to the frequencies
    :param freqs: list of numbers, not necessarily adding to 1
    :return: integer
    """
    t = len(freqs)
    array = np.array(freqs)
    if sum(array) != 1:
        array = array/sum(array)
        return np.random.choice(t, 1, p=array)[0]
    else:
        return np.random.choice(t, 1, p=array)[0]


def gibbs_sampler(dna, k, N):
    """
    N iterations of Gibbs sampling. Starting with random motif set,
    pick a string randomly and get a random new kmer based off the
    probabilities in the profile using gibbs_random function
    :param dna: list of equal length DNA strings
    :param k: length of desired kmers
    :param N: number of iterations of Gibbs sampler
    :return: best motif set
    """
    t = len(dna[0])
    slice = lambda x, y: x[y:y+k]
    bestmotifs = [slice(line, random.randint(0, t - k)) for line in dna]
    bestprof = motifmatrix(bestmotifs)
    bestscore = scoreprofile(bestprof, len(dna))
    motifs = bestmotifs.copy()
    for j in range(N):
        i = random.randint(0, len(dna)-1)
        motifs.pop(i)
        prof = motifmatrix(motifs)
        randkmer = profilemost_random(dna[i], prof)
        motifs.insert(i, randkmer)
        profafterinsert = motifmatrix(motifs)
        if scoreprofile(profafterinsert, len(dna)) < bestscore:
            bestmotifs = motifs
            bestprof = profafterinsert
            bestscore = scoreprofile(bestprof, len(dna))
    return bestmotifs, bestscore