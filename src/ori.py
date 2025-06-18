# finding origin of replication
import math
from typing import Dict, Any


def pattern_count(text, pattern):
    """
    Counts the instances of a kmer pattern in a string.
    Includes overlapping patterns.
    :param text: string to search
    :param pattern: pattern to find in string
    :return: count of pattern repeats in string
    """
    count = 0
    plen = len(pattern)
    tlen = len(text)
    for i in range((tlen - plen) + 1):
        if text[i:i + plen] == pattern:
            count += 1
    return count


def most_frequent(text, k):
    """
    identifies the most frequent pattern of length k
    in a string. can return multiple patterns
    :param text: string to search
    :param k: integer length of kmer to find
    :return: list of most frequent patterns
    """
    freqs = {}
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        if kmer not in freqs:  # check that kmer hasn't been counted
            subcount = pattern_count(text, kmer)
            freqs[kmer] = subcount
    maxfreq = max(freqs.values())
    maxlist = [p for p in freqs if freqs[p] == maxfreq]
    return maxlist


def frequency_table(text, k):
    """
    For faster calculation of frequent kmers.
    Makes frequency table of patterns of length k
    :param text: string to search
    :param k: integer length of kmers
    :return: dictionary of kmers and their frequencies
    """
    freqmap = {}
    n = len(text)
    for i in range((n - k) + 1):
        pattern = text[i:i + k]
        freqmap[pattern] = freqmap.get(pattern, 0) + 1
    return freqmap


def better_frequent_kmer(text, k):
    """
    possibly faster way of finding most frequent kmers
    :param text: string to search
    :param k: integer length of kmers
    :return: list of most frequent kmers
    """
    freqmap = frequency_table(text, k)
    maxval = max(freqmap.values())
    mostcommon = [i for i in freqmap if freqmap[i] == maxval]
    return mostcommon


# Probability functions
# Calculating the probability of certain patterns occurring in a string

def pr(N, a, pattern, t=1):
    """
    Finds the probability of a pattern occurring in a random string
    :param N: random string length
    :param a: number of letters in the alphabet
    :param pattern: pattern string to match
    :param t: minimum number of matches
    :return: probability
    """
    letters = [str(x) for x in range(a)]
    pools = [tuple(letters)] * N  # generate all possible combos
    combos = [[]]
    counter = 0
    for pool in pools:
        combos = [x + [y] for x in combos for y in pool]
    combolist = ["".join(p for p in prod) for prod in combos]
    numbercombos = len(combolist)
    for combo in combolist:
        if pattern_count(combo, str(pattern)) >= t:
            counter += 1
    print(f"{numbercombos} combos, {counter} counts, {counter / numbercombos} proportion")


def count_patterns(N, p, t=1):
    """
    number of ways that t copies of a non-overlapping pattern
    can appear in a string of length n
    :param N: length of text
    :param p: length of pattern
    :param t: number of times to match
    :return:
    """
    nonplength = N - (p * t)  # set up binomial
    binom = math.factorial(t + nonplength) / (math.factorial(t) * math.factorial(nonplength))
    return binom


# DNA functions
def reverse_complement(string):
    dnadict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    string = string.upper()
    complement = ''
    for i in range(len(string)-1, -1, -1):
        complement += dnadict[string[i]]
    return complement

def pattern_index(text, pattern):
    indexlist = []
    plen = len(pattern)
    tlen = len(text)
    for i in range(tlen - plen + 1):
        if text[i:i + plen] == pattern:
            indexlist.append(i)
    return indexlist


def dna_combinations(repeat=3):
    """
    based on itertools product
    :param repeat: length of substring combos
    :return: returns generator of substrings
    """
    pools = [tuple('ACGT')] * repeat
    result = [[]]
    for pool in pools:
        result = [x + [y] for x in result for y in pool]
    for prod in result:
        yield "".join(p for p in prod)


def frequent_kmer(text, k=5):
    subgen = dna_combinations(repeat=k)
    freqlist = []
    for pattern in subgen:
        count = pattern_count(text, pattern)
        freqlist.append((count, pattern))
    mostfreq = max([c[0] for c in freqlist])
    freqlist = sorted(freqlist, reverse=True)
    for t in freqlist:
        if t[0] == mostfreq:
            print(t)


def clump_finder(genome, k, L, t):
    """
    finds any patterns that form an (L, t) clump
    :param genome: string to search on
    :param k: size of k-mers
    :param L: window to look for clumps
    :param t: number of repeats within window
    :return: distinct k-mers over threshold
    """
    numwindows = len(genome) - L + 1
    print(numwindows)
    totalset = set()
    for i in range(numwindows):
        window = genome[i:i + L]
        freqmap = frequency_table(window, k)
        for key in freqmap:
            if freqmap[key] >= t:
                totalset.add(key)
    return totalset


### got this from github
def ClumpFinder(genome, L, k, t):
    """
    We defined a k-mer as a "clump" if it appears many times within a short interval of the genome.
    More formally, given integers L and t, a k-mer Pattern forms an (L, t)-clump inside a (longer)
    string Genome if there is an interval of Genome of length L in which this k-mer appears at least t times.
    """
    pattern_map = {}
    result = []
    for i in range(len(genome) - k + 1):
        pattern = genome[i:i + k]
        if pattern in pattern_map:
            pattern_map[pattern].append(i)
        else:
            pattern_map[pattern] = [i]

    for pattern, arr in pattern_map.items():
        if len(arr) < t:
            continue
        for i in range(len(arr) - t + 1):
            if arr[i + t - 1] - arr[i] + k <= L:
                result.append(pattern)
                break

    return set(result)


## also from github
def BetterClumpFinder(genome, k, L, t):
    """Clump finder that is efficient enough to search the entire genome of E Coli."""

    n = len(genome) - k + 1
    m = L - k + 1
    pattern_count = {}
    pattern_list = []
    kmers = []

    for i in range(n):
        pattern_list.append(genome[i:i + k])
        pattern_count[genome[i:i + k]] = 0

    for i in range(m):
        pattern_count[genome[i:i + k]] += 1

    for i in pattern_count:
        if pattern_count[i] == 3:
            kmers.append(i)

    for i in range(m, n):
        pattern_count[genome[i:i + k]] += 1
        pattern = pattern_list[i - m]
        pattern_count[pattern] -= 1
        if pattern_count[genome[i:i + k]] == 3:
            kmers.append(genome[i:i + k])

    return len(list(dict.fromkeys(kmers)))
