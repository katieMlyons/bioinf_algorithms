from src.ori import frequency_table, reverse_complement
from typing import Union


def skew(string:str) -> list[int]:
    '''running score of G - C skew across genome'''
    skewcount = 0
    skewlist = []
    skewlist.append(skewcount)
    for char in string:
        if char == "G":
            skewcount += 1
        elif char == "C":
            skewcount -= 1
        skewlist.append(skewcount)
    return skewlist


def min_skew(genome:str) -> list[int]:
    '''indices of minimum G-C skew in a genomic sequence'''
    skewlist = skew(genome)
    minskew = min(skewlist)
    minlist = []
    for i, val in enumerate(skewlist):
        if val == minskew:
            minlist.append(i)
    return minlist


def hamming_distance(p:str, q:str) -> int:
    if len(p) != len(q):
        return False
    else:
        return sum(c1 != c2 for c1, c2 in zip(p, q))


def approximate_pattern_matching(pattern:str, text:str, d:int) -> list[int]:
    '''list of indices of kmers within hamming distance d of a pattern'''
    kmers = set() #don't need to recalculate good kmers
    idxlist = []
    k = len(pattern)
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        if kmer in kmers:
            idxlist.append(i)
        elif hamming_distance(kmer, pattern) <= d:
            kmers.add(kmer)
            idxlist.append(i)
    return idxlist


def approximate_pattern_count(text:str, pattern:str, d:int) -> int:
    '''count of kmers in a sequence within hamming distance d of a pattern'''
    kmers = set() # don't need to recalculate good kmers
    count = 0
    k = len(pattern)
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        if kmer in kmers:
            count += 1
        elif kmer == pattern or hamming_distance(kmer, pattern) <= d:
            kmers.add(kmer)
            count += 1
    return count


def neighbors(pattern:str, d:int) -> list[str]:
    '''set of all strings within d Hamming distance of a pattern'''
    if d == 0:
        return [pattern]
    elif len(pattern) == 1:
        return ['A', 'C', 'G', 'T']
    else:
        neighborhood = set()
        suffixneighbors = neighbors(pattern[1:], d)
        for text in suffixneighbors:
            if hamming_distance(pattern[1:], text) < d:
                for nuc in ['A', 'C', 'G', 'T']:
                    neighborhood.add(nuc + text)
            else:
                neighborhood.add(pattern[0] + text)
        return list(neighborhood)

def frequent_words_with_mismatches(text:str, k:int, d:int) -> list[str]:
    '''returns a list of kmers with the maximum number of d-neighbors
    in a sequence'''
    freqmap = {}
    for i in range(len(text) + k - 1):
        kmer = text[i: i+k]
        neighborhood = neighbors(kmer, d)
        for neighbor in neighborhood:
            freqmap[neighbor] = freqmap.get(neighbor, 0) + 1
    m = max(freqmap.values())
    patterns = [p for p in freqmap if freqmap[p] == m]
    return patterns
    

def frequent_words_with_mismatches_rc(text, k, d):
    freqmap = frequency_table(text, k)
    pool = set()
    resultdict = {}
    for pattern in freqmap:
        if d > 0:
            pool.update(neighbors(pattern, d))
            pool.update(neighbors(reverse_complement(pattern), d))
        else:
            pool.add(neighbors(pattern, d))
            pool.add(neighbors(reverse_complement(pattern), d))
    for tester in pool:
        for map in freqmap:
            if (hamming_distance(tester, map) <= d): # or (hamming(tester, reverse_complement(map)) <= d):
                resultdict[tester] = resultdict.get(tester, 0) + (1*freqmap[map])
    finallist = []
    for result in resultdict:
        if reverse_complement(result) in resultdict:
            finallist.append((resultdict[result]+resultdict[reverse_complement(result)], result, reverse_complement(result)))
        else:
            finallist.append((resultdict[result], result, " "))
    maxval = max(x[0] for x in finallist)
    condensed = [x for x in finallist if x[0] == maxval]
    return condensed


