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
    text_len = len(text)
    pattern_len = len(pattern)

    for i in range(text_len - pattern_len + 1):
        j = 0
        while j < pattern_len and text[i + j] == pattern[j]:
            j += 1

        if j == pattern_len:
            matches.append(i)

    return matches

text = "Hello, how are you? Hello, nice to meet you."
pattern = "Hello"

matches = naiveTextSearch(text, pattern)
print(matches)  # Output: [0, 19]