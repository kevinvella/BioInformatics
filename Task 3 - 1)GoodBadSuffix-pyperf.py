import pyperf
import matplotlib.pyplot as plt

def load_sequences(file_path):
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


def preprocess_sequences(sequences):
    preprocessed_sequences = []
    for sequence in sequences:
        # Remove any unwanted characters, convert to uppercase, etc.
        preprocessed_sequence = sequence.upper()
        preprocessed_sequences.append(preprocessed_sequence)
    return preprocessed_sequences


def create_word_index(sequences, word_length):
    word_index = {}
    for i, sequence in enumerate(sequences):
        for j in range(len(sequence) - word_length + 1):
            word = sequence[j:j+word_length]
            if word in word_index:
                word_index[word].append((i, j))
            else:
                word_index[word] = [(i, j)]
    return word_index


def extract_query_words(query_sequence, word_length):
    query_words = []
    for i in range(len(query_sequence) - word_length + 1):
        query_words.append(query_sequence[i:i+word_length])
    return query_words


def match_query_words(query_words, word_index):
    matching_sequences = set()
    for word in query_words:
        if word in word_index:
            matching_sequences.update([match[0] for match in word_index[word] if match[0] not in matching_sequences])
    return matching_sequences


def score_sequences(matching_sequences, word_index, query_words):
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




def rank_sequences(scores):
    ranked_sequences = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    return ranked_sequences

def benchmark_naiveTextSearch(loops):
    text = "Lorem dolor sit amet."
    pattern = "ipsum"

    for _ in range(loops):
        # Main execution
        query_file = 'query.fa'
        target_file = 'target.fa'
        word_length = 3

        # Load sequences
        query_sequence = load_sequences(query_file)[0]
        target_sequences = load_sequences(target_file)

        # Preprocess sequences
        query_sequence = preprocess_sequences([query_sequence])[0]
        target_sequences = preprocess_sequences(target_sequences)

        # Create word index
        word_index = create_word_index(target_sequences, word_length)

        # Extract query words
        query_words = extract_query_words(query_sequence, word_length)

        # Match query words with target sequences
        matching_sequences = match_query_words(query_words, word_index)

        # Score sequences based on word matches
        scores = score_sequences(matching_sequences, word_index, query_words)

        # Rank sequences based on scores
        ranked_sequences = rank_sequences(scores)

runner = pyperf.Runner()
results = runner.bench_func('benchmark_goodBadSuffixTextSearch', benchmark_naiveTextSearch, 10)

results.dump('benchmark_goodBadSuffixTextSearch.json',replace=True)