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
    preprocessedSequences = []
    for sequence in sequences:
        # Remove any unwanted characters, convert to uppercase, etc.
        preprocessedSequence = sequence.upper()
        preprocessedSequences.append(preprocessedSequence)
    return preprocessedSequences


def createWordIndex(sequences, wordLength):
    wordIndex = {}
    for i, sequence in enumerate(sequences):
        for j in range(len(sequence) - wordLength + 1):
            word = sequence[j:j+wordLength]
            if word in wordIndex:
                wordIndex[word].append((i, j))
            else:
                wordIndex[word] = [(i, j)]
    return wordIndex


def extractWords(querySequence, wordLength):
    query_words = []
    for i in range(len(querySequence) - wordLength + 1):
        query_words.append(querySequence[i:i+wordLength])
    return query_words


def matchWords(queryWords, wordIndex):
    matchingSequences = set()
    for word in queryWords:
        if word in wordIndex:
            matchingSequences.update([match[0] for match in wordIndex[word] if match[0] not in matchingSequences])
    return matchingSequences


def scoreSequences(matchingSequences, wordIndex, queryWords):
    scores = {}
    for sequenceId in matchingSequences:
        sequenceMatches = []
        for word in queryWords:
            if word in wordIndex:
                sequenceMatches.extend([match[1] for match in wordIndex[word] if match[0] == sequenceId])
        score = 0
        for queryMatch in queryWords:
            for sequenceMatch in sequenceMatches:
                if sequenceMatch <= sequenceMatch + len(queryMatch):
                    score += 1
        scores[sequenceId] = score
    return scores




def rank(scores):
    ranked_sequences = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    return ranked_sequences


query_file = 'query.fa'
target_file = 'target.fa'
wordLength = 3
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
        wordLength = (args.wordlength)

try:
    # Load sequences
    querySequence = load(query_file)[0]
    target_sequences = load(target_file)

    # Preprocess sequences
    querySequence = preSequences([querySequence])[0]
    target_sequences = preSequences(target_sequences)

    # Create word index
    wordIndex = createWordIndex(target_sequences, wordLength)

    # Extract query words
    queryWords = extractWords(querySequence, wordLength)

    # Match query words with target sequences
    matchingSequences = matchWords(queryWords, wordIndex)

    # Score sequences based on word matches
    scores = scoreSequences(matchingSequences, wordIndex, queryWords)

    # Rank sequences based on scores
    rankedSequences = rank(scores)

    # Output the ranking
    print("Ranking of target sequences:")
    for sequenceId, score in rankedSequences:
        print(f"Sequence ID: {sequenceId}, Score: {score}")

except Exception as e:
    print("An error has occured while doing the BLAST search")
    print(e)