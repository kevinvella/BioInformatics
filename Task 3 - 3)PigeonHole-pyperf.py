import pyperf
import matplotlib.pyplot as plt

def bad_character_table(pattern):
    table = {}
    for i in range(len(pattern)):
        table[pattern[i]] = i
    return table


def good_suffix_table(pattern):
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


def boyer_moore(text, pattern, n):
    m = len(pattern)
    bc_table = bad_character_table(pattern)
    gs_table = good_suffix_table(pattern)
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

def benchmark_naiveTextSearch(loops):

    for _ in range(loops):
        # Main execution
        # Example usage
        text = "ABACADABACABA"
        pattern = "DAB"
        mismatches = 1
        boyer_moore(text, pattern, mismatches)

runner = pyperf.Runner()
results = runner.bench_func('benchmark_pigeonHoleTextSearch', benchmark_naiveTextSearch, 10)

results.dump('benchmark_pigeonHoleTextSearch.json',replace=True)