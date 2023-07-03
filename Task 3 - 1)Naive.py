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

text = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer sodales, nibh ut laoreet congue, nulla eros tempus ante, quis mattis ipsum odio nec quam. Fusce et magna imperdiet, convallis libero nec, efficitur neque. Praesent porttitor, velit sit amet laoreet posuere, felis augue rhoncus quam, quis venenatis est tellus nec nulla. Nullam volutpat est ac lacus consequat pellentesque. Sed aliquam orci id bibendum ultrices. Nulla facilisis elit nec neque cursus rutrum. Nunc tincidunt ullamcorper cursus. Orci varius natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus."
pattern = "ipsum"

matches = naiveTextSearch(text, pattern)
print(matches)  # Output: [6, 134]