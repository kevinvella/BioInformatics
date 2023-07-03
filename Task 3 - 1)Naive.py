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

text = "ipsumipsumLorem ipsum dolor sit amet. ipsum"
pattern = "ipsum"

matches = naiveTextSearch(text, pattern)
print(matches)