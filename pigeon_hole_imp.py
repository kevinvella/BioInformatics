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
    n = len(text)
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


# Example usage
text = "ABACADABACABA"
pattern = "DAB"
mismatches = 1
matches = boyer_moore(text, pattern, mismatches)
if matches:
    print(f"Pattern found at indices: {matches}")
    for match in matches:
        print(f"Matched substring: {text[match:match + len(pattern)]}")
else:
    print("Pattern not found")
