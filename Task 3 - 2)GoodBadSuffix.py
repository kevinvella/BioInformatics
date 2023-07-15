import argparse
import sys

def badCharacterTable(pattern):
    table = {}
    for i in range(len(pattern)):
        table[pattern[i]] = i
    return table


def goodSuffixTable(pattern):
    m = len(pattern)
    table = [0] * (m + 1)
    suffix = [0] * (m + 1)

    for i in range(m):
        suffix[i] = m
    j = m
    for i in range(m - 1, -1, -1):
        if j == m - i:
            while j < m and pattern[(i + j) - 1] == pattern[j]:
                j += 1
            suffix[j] = m - i - j
        j -= 1

    for i in range(m):
        table[i] = m - suffix[i]

    for i in range(m - 1):
        table[m - suffix[i]] = m - i - 1

    return table


def boyerMoore(text, pattern):
    m = len(pattern)
    n = len(text)
    bcTable = badCharacterTable(pattern)
    gsTable = goodSuffixTable(pattern)
    i = 0

    while i <= n - m:
        j = m - 1
        while j >= 0 and pattern[j] == text[i + j]:
            j -= 1

        if j == -1:
            # Match found
            return i

        if text[i + j] in bcTable:
            bc_offset = bcTable[text[i + j]]
        else:
            bc_offset = -1

        gs_offset = gsTable[j + 1]

        i += max(gs_offset, j - bc_offset)

    return -1


text = ""
pattern = ""
if sys.gettrace() is None:
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Good Bad Suffix Search')
    parser.add_argument('pattern', help='The pattern that you want to search')
    parser.add_argument('text', help='The text you want to search in')
    args = parser.parse_args()

    text = args.text #"Lorem dolor sit amet. ipsu"
    pattern = args.pattern # "ipsum"
else:
    text = "ABACADABACABA"
    pattern = "DABA"


try:
    result = boyerMoore(text, pattern)
    if result != -1:
        print(f"Pattern found at index {result}")
        print(f"Matched substring: {text[result:result+len(pattern)]}")
    else:
        print("Pattern not found")
except:
    print("An error has occured while doing the native text search")
