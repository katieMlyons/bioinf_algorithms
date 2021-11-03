# Chapter 1: looking at a bacterial ori with string functions
import math


def pattern_count(text, pattern):
    count = 0
    plen = len(pattern)
    tlen = len(text)
    for i in range((tlen - plen) + 1):
        if text[i:i + plen] == pattern:
            count += 1
    return count


def pattern_index(text, pattern):
    indexlist = []
    plen = len(pattern)
    tlen = len(text)
    for i in range(tlen - plen + 1):
        if text[i:i + plen] == pattern:
            indexlist.append(i)
    for i in indexlist:
        print(i, end=' ')
    print('\n')


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


def frequency_table(text, k=5):
    freqmap = {}
    n = len(text)
    for i in range((n - k) + 1):
        pattern = text[i:i + k]
        freqmap[pattern] = freqmap.get(pattern, 0) + 1
    return freqmap


def better_frequent_kmer(text, k=5):
    mostcommon = []
    freqmap = frequency_table(text, k)
    maxval = max(freqmap.values())
    for substring in freqmap:
        if freqmap[substring] == maxval:
            mostcommon.append(substring)
    return mostcommon


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


def pr(N, a, pattern, t=1):
    """
    Finds the probability of a pattern occurring in a random string
    :param N: random string length
    :param a: number of letters in the alphabet
    :param pattern: pattern to match
    :param t: minimum number of matches
    :return: probability
    """
    letters = [str(x) for x in list(range(a))]
    pools = [tuple(letters)] * N  # generate all possible combos
    combos = [[]]
    counter = 0
    for pool in pools:
        combos = [x + [y] for x in combos for y in pool]
    combolist = []
    for prod in combos:
        combolist.append("".join(p for p in prod))
    numbercombos = len(combolist)
    for combo in combolist:
        if pattern_count(combo, str(pattern)) >= t:
            counter += 1
    print(f"{numbercombos} combos, {counter} counts, {counter / numbercombos} proportion")


def pr_approx(N, a, pattern, t=1):
    """
    approximates the probability of a pattern appearing in a random string
    :param N: length of text
    :param a: number of letters in the alphabet
    :param pattern: pattern to match
    :param t: number of times to match
    :return:
    """
    plength = len(pattern)
    nonplength = N - (plength * t)  # set up binomial
    binom = math.factorial(t + nonplength) / (math.factorial(t) * math.factorial(nonplength))
    return binom


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
