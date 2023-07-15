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


def boyerMoore(text, pattern, n):
    m = len(pattern)
    bc_table = badCharacterTable(pattern)
    gs_table = goodSuffixTable(pattern)
    i = 0

    matches = []
    while i <= n - m:
        j = m - 1
        mismatch_count = 0
        while j >= 0 and pattern[j] == text[i + j]:
            j -= 1

        if j == -1:
            # Match found
            matches.append(i)
            i += gs_table[0]
        else:
            if text[i + j] in bc_table:
                bc_offset = bc_table[text[i + j]]
            else:
                bc_offset = -1

            gs_offset = gs_table[j + 1]

            i += max(gs_offset, j - bc_offset)

        mismatch_count += 1
        if mismatch_count > n:
            break

    return matches


text = ""
pattern = ""
mismatches = 1
if sys.gettrace() is None:
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Good Bad Suffix Search')
    parser.add_argument('pattern', help='The pattern that you want to search')
    parser.add_argument('text', help='The text you want to search in')
    parser.add_argument('--mismatches', type=int, help='Threshold for mismatches. Default is 1')
    args = parser.parse_args()

    text = args.text #"Lorem dolor sit amet. ipsu"
    pattern = args.pattern # "ipsum"
    
    if args.mismatches is not None:
        mismatches = (args.mismatches)
else:
    text = "ABACADABACABA"
    pattern = "DABA"


try:
    result = boyerMoore(text, pattern, mismatches)
    if len(result) > 0:
        print(f"Pattern found at index {result}")
        print(f"Matched substring: {text[result:result+len(pattern)]}")
    else:
        print("Pattern not found")
except Exception as e:
    print("An error has occured while doing the pigeon hole search")
    print(e)