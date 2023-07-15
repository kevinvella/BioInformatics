import pyperf
import matplotlib.pyplot as plt

def naiveTextSearch(text, pattern):
    """
    Perform a naive text search for exact matching.

    Args:
        text (str): The text to search within.
        pattern (str): The pattern to search for.

    Returns:
        list: A list of starting positions of matches.
    """
    matches = []
    textLength = len(text)
    patternLength = len(pattern)

    for i in range(textLength):
        j = 0
        while j < patternLength and text[i + j] == pattern[j]:
            j += 1

        if j == patternLength:
            matches.append(i)

    return matches

def benchmarkNaiveTextSearch(loops):
    text = "Lorem dolor sit amet."
    pattern = "ipsum"

    for _ in range(loops):
        naiveTextSearch(text, pattern)

runner = pyperf.Runner()
results = runner.bench_func('benchmark_naiveTextSearch', benchmarkNaiveTextSearch, 10)

results.dump('benchmark_naiveTextSearch.json',replace=True)