from bin.ori import frequency_table
dnadict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def reverse_complement(string):
    string = string.upper()
    complement = ''
    for i in range(len(string)-1, -1, -1):
        complement += dnadict[string[i]]
    return complement


def skew(string):
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


def min_skew(string):
    skewlist = skew(string)
    minskew = min(skewlist)
    for i, val in enumerate(skewlist):
        if val == minskew:
            print(i, end=" ")
    print("\n")


def hamming(str1, str2):
    if len(str1) != len(str2):
        return False
    counter = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            counter += 1
    return counter


def count_hamming(text, pattern, d):
    patlength = len(pattern)
    counter = 0
    for i in range(len(text)-patlength+1):
        teststring = text[i:i+patlength]
        if hamming(teststring, pattern) <= d:
            counter += 1
    return counter


def frequent_words_with_mismatches(text, k, d):
    freqmap = frequency_table(text, k)
    patterndict = {}
    for pattern in freqmap:
        patterndict[pattern] = freqmap[pattern]
        for pattern2 in freqmap:
            if pattern != pattern2:
                if hamming(pattern, pattern2) <= d:
                    patterndict[pattern] += freqmap[pattern2]
    maxval = max(patterndict.values())
    print(maxval)
    for entry in patterndict:
        if patterndict[entry] == maxval:
            print(entry, end = " ")
    print("\n")


def immediate_neighbors(pattern, d):
    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return ['A', 'C', 'G', 'T']
    neighborhood = set()
    suffixneighbors = immediate_neighbors(pattern[1:], d)
    for stri in suffixneighbors:
        if hamming(pattern[1:], stri) < d:
            for nuc in ['A', 'C', 'G', 'T']:
                neighborhood.add(nuc+stri)
        else:
            neighborhood.add(pattern[0]+stri)
    return neighborhood


def frequent_words_with_mismatches_rc(text, k, d):
    freqmap = frequency_table(text, k)
    pool = set()
    resultdict = {}
    for pattern in freqmap:
        if d > 0:
            pool.update(immediate_neighbors(pattern, d))
            pool.update(immediate_neighbors(reverse_complement(pattern), d))
        else:
            pool.add(immediate_neighbors(pattern, d))
            pool.add(immediate_neighbors(reverse_complement(pattern), d))
    for tester in pool:
        for map in freqmap:
            if (hamming(tester, map) <= d): # or (hamming(tester, reverse_complement(map)) <= d):
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


