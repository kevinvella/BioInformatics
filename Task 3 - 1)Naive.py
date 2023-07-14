import argparse
import sys

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

text = ""
pattern = ""
if sys.gettrace() is None:
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Naive Search')
    parser.add_argument('pattern', help='The pattern that you want to search')
    parser.add_argument('text', help='The text you want to search in')
    args = parser.parse_args()

    text = args.text #"Lorem dolor sit amet. ipsu"
    pattern = args.pattern # "ipsum"
else:
    text = "Lorem ipsum dolor sit amet"
    pattern = "ipsum"

try:
    matches = naiveTextSearch(text, pattern)
    print(matches)
except:
    print("An error has occured while doing the native text search")