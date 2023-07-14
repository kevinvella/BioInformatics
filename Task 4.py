import argparse
import sys

def load(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence)
    return sequences


def preSequences(sequences):
    preprocessed_sequences = []
    for sequence in sequences:
        # Remove any unwanted characters, convert to uppercase, etc.
        preprocessed_sequence = sequence.upper()
        preprocessed_sequences.append(preprocessed_sequence)
    return preprocessed_sequences


def createWordIndex(sequences, word_length):
    word_index = {}
    for i, sequence in enumerate(sequences):
        for j in range(len(sequence) - word_length + 1):
            word = sequence[j:j+word_length]
            if word in word_index:
                word_index[word].append((i, j))
            else:
                word_index[word] = [(i, j)]
    return word_index


def extractWords(query_sequence, word_length):
    query_words = []
    for i in range(len(query_sequence) - word_length + 1):
        query_words.append(query_sequence[i:i+word_length])
    return query_words


def matchWords(query_words, word_index):
    matching_sequences = set()
    for word in query_words:
        if word in word_index:
            matching_sequences.update([match[0] for match in word_index[word] if match[0] not in matching_sequences])
    return matching_sequences


def scoreSequences(matching_sequences, word_index, query_words):
    scores = {}
    for sequence_id in matching_sequences:
        sequence_matches = []
        for word in query_words:
            if word in word_index:
                sequence_matches.extend([match[1] for match in word_index[word] if match[0] == sequence_id])
        score = 0
        for query_match in query_words:
            for sequence_match in sequence_matches:
                if sequence_match <= sequence_match + len(query_match):
                    score += 1
        scores[sequence_id] = score
    return scores




def rank(scores):
    ranked_sequences = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    return ranked_sequences


# Example usage
query_file = 'query.fa'
target_file = 'target.fa'
word_length = 3
if sys.gettrace() is None:
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='MINI BLAST Arguments')
    parser.add_argument('query', help='The query file you want to search in')
    parser.add_argument('target', help='The target file you want to search in')
    parser.add_argument('--wordlength', type=int, help='Word length. Default is 3')
    args = parser.parse_args()

    query_file = args.query
    target_file = args.target
    
    if args.wordlength is not None:
        word_length = (args.wordlength)

try:
    # Load sequences
    query_sequence = load(query_file)[0]
    target_sequences = load(target_file)

    # Preprocess sequences
    query_sequence = preSequences([query_sequence])[0]
    target_sequences = preSequences(target_sequences)

    # Create word index
    word_index = createWordIndex(target_sequences, word_length)

    # Extract query words
    query_words = extractWords(query_sequence, word_length)

    # Match query words with target sequences
    matching_sequences = matchWords(query_words, word_index)

    # Score sequences based on word matches
    scores = scoreSequences(matching_sequences, word_index, query_words)

    # Rank sequences based on scores
    ranked_sequences = rank(scores)

    # Output the ranking
    print("Ranking of target sequences:")
    for sequence_id, score in ranked_sequences:
        print(f"Sequence ID: {sequence_id}, Score: {score}")
except Exception as e:
    print("An error has occured while doing the BLAST search")
    print(e)