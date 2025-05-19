import sys
from itertools import product
import numpy as np
import random



# O(k^2)
# output: k^d
def neighbors(pattern, distance):
    if distance == 0:
        return {pattern}
    if len(pattern) == 1:
        return {"A", "C", "G", "T"}
    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], distance)
    for text in suffix_neighbors:
        if (hamming_distance(pattern[1:], text) < distance):
            for x in {"A","C","G","T"}:
                neighborhood.add(x + text)
        else:
            neighborhood.add(pattern[0] + text)
    return neighborhood



# O(n*k)
def approx_pattern_occurs(text, pattern, distance):
    for i in range(0,len(text) - len(pattern) + 1):
        window = text[i : i + len(pattern)]
        if hamming_distance(pattern, window) <= distance:
            return True
    return False



# O(k)
def hamming_distance(str1, str2):
    if (len(str1) != len(str2)):
        return -1
    mismatches = 0
    for i in range(0, len(str1)):
        if (str1[i] != str2[i]):
            mismatches += 1
    return mismatches



# O(n*k^d*t*n*k) = O(n^2*t*k^(d+1))
# checks all kmers in strand 1 to see if they are in the other strands
def motif_enumeration(dna, k, d):
    all_patterns = set()
    for i in range(0, len(dna[0]) - k + 1):
        pattern = dna[0][i : i + k]
        for approx_pattern in neighbors(pattern, d): # O(k^2) k^d
            for strand in dna: # t
                if approx_pattern_occurs(strand, approx_pattern, d):# O(n*k)
                    all_patterns.add(approx_pattern)
                    print(approx_pattern)
    return all_patterns



# O(4^k)
def generate_all_kmers(k):
    alphabet = ["A", "G", "C", "T"]
    return [ ''.join(x) for x in product(alphabet, repeat=k) ]



# O(t*n*k)
# finds the motif in each strand that is closest to pattern and returns the score
def total_min_score(pattern, dna):
    score = 0
    k = len(pattern)
    for text in dna:
        min_distance = hamming_distance(pattern, text[0:k])# O(k)
        for i in range(1, len(text) - k + 1):
            distance = hamming_distance(pattern, text[i:i+k])# O(k)
            if distance < min_distance:
                min_distance = distance
        score += min_distance
    return score



# O(n^2*t^2*k)
# brute force: checks all patterns in dna
def median_string(dna, k):
    min_score = sys.maxsize
    median = ""
    dna_string = "".join(dna)
    for i in range(0, len(dna_string) - k + 1):# t*n
        pattern = dna_string[i : i + k]
        score = total_min_score(pattern, dna)# O(t*n*k)
        if score < min_score:
            min_score = score
            median = pattern
    return median



# O(4^k*t*n*k)
def median_string_all_kmers(dna, k):
    min_score = sys.maxsize
    median = ""
    all_kmers = generate_all_kmers(k)# O(4^k)
    for kmer in all_kmers:
        score = total_min_score(kmer, dna)# O(t*n*k)
        if score < min_score:
            min_score = score
            median = kmer
    return median



# O(t*k)
# creates a count matrix for given motifs
def count_matrix(motifs):
    profile_matrix = np.zeros((4, len(motifs[0])))
    for i in range(0, len(motifs[0])):
        for j in range(0, len(motifs)):
            match (motifs[j][i]):
                case "A":
                    profile_matrix[0][i] += 1
                case "C":
                    profile_matrix[1][i] += 1
                case "G":
                    profile_matrix[2][i] += 1
                case "T":
                    profile_matrix[3][i] += 1
                case _:
                    print("not a valid character")
                    return
    return profile_matrix



# O(t*k)
# creates a profile matrix for given motifs
def profile(motifs):
    profile_matrix = count_matrix(motifs)# O(t*k)
    for i in range(0, len(profile_matrix)):
        for j in range(0, len(profile_matrix[0])):
            profile_matrix[i][j] = (profile_matrix[i][j] + 1) / (len(motifs) + 4)
    return profile_matrix



# O(t*k)
# creates a consensus string for given motifs
def consensus(motifs):
    consensus_string = ""
    scores = count_matrix(motifs)# O(t*k)
    for i in range(0, len(scores[0])):
        max_score = 0
        max_index = 0
        for j in range(0, len(scores)):
            if scores[j][i] > max_score:
                max_score = scores[j][i]
                max_index = j
        match (max_index):
            case 0:
                consensus_string += "A"
            case 1:
                consensus_string += "C"
            case 2:
                consensus_string += "G"
            case 3:
                consensus_string += "T"
    return consensus_string
        


# O(k)
# generates the probability of a pattern based on the profile matrix
def string_probability(pattern, profile):
    probability = 1
    for i in range(0, len(pattern)):
        match (pattern[i]):
            case "A":
                probability = probability * profile[0][i]
            case "C":
                probability = probability * profile[1][i]
            case "G":
                probability = probability * profile[2][i]
            case "T":
                probability = probability * profile[3][i]
            case _:
                print("not a valid character")
    return probability



# O(n*k)
# generates the profile most probable kmer in one strand
def profile_most_probable_kmer(text, k, profile):
    max_probability = 0
    best_string = text[:k]
    for i in range(0, len(text) - k + 1):
        probability = string_probability(text[i:i + k], profile)# O(k)
        if probability > max_probability:
            max_probability = probability
            best_string = text[i:i+k]
    return best_string



# O(t*k)
# generates the score of a list of motifs
def score(motifs):
    score = 0
    consensus_string = consensus(motifs)# O(t*k)
    for motif in motifs:
        score += hamming_distance(consensus_string, motif)# O(k)
    return score



# O(n^2*t*k)
# for every kmer in the first strand, this adds the profile most probable kmer then adds that to the profile matrix
def greedy_motif_search(dna, k):
    best_motifs = []
    for line in dna:
        best_motifs.append(line[:k])

    for i in range(0, len(dna[0]) - k + 1):# n
        motif1 = dna[0][i:i+k]
        motifs = [motif1]
        for j in range(1,len(dna)):# t
            profile_matrix = profile(motifs)# O(t*k)
            new_motif = profile_most_probable_kmer(dna[j], k, profile_matrix)# O(n*k)
            motifs.append(new_motif)
        if score(motifs) < score(best_motifs):# O(t*k)
            best_motifs = motifs
    return best_motifs



# O(n*t*k)
# creates a list of profile most probable motifs in a group of strands
def motif_method(profile, dna):
    k = len(profile[0])
    motif_list = []
    for text in dna:# t
        motif_list.append(profile_most_probable_kmer(text, k, profile))# O(n*k)
    return motif_list



# O(r*n*t*k) where r is repeats, which is random and different every time, honestly though it can just be treated as a constant
# repeatedly generates all profile most probable motifs, then makes a new profile from that, and this repeats
def randomized_motif_search(dna, k):
    motif_list = []
    for line in dna:# t
        random_index = random.randrange(0, len(dna[0]) - k + 1)
        motif_list.append(line[random_index : random_index + k])
    best_motifs = motif_list[:]
    while True:
        profile_matrix = profile(motif_list)# O(t*k)
        motif_list = motif_method(profile_matrix, dna)# O(n*t*k)
        if score(motif_list) < score(best_motifs):# O(t*k)
            best_motifs = motif_list[:]
        else:
            return best_motifs



# O(i*r*n*t*k)
def iterate_randomized_motif_search(dna, k, iterations):
    best_score = score(randomized_motif_search(dna, k))# O(r*n*t*k)
    best_motifs = []
    for i in range(iterations):
        current_motifs = randomized_motif_search(dna, k)# O(r*n*t*k)
        current_score = score(current_motifs)# O(t*k)
        if current_score < best_score:
            best_score = current_score
            best_motifs = current_motifs
    return best_motifs



# O(n)
def random_with_probability_distribution(probabilities):
    total = sum(probabilities)
    new_probabilities = [probability / total for probability in probabilities]
    random_number = random.choices(range(0,len(new_probabilities)), weights=new_probabilities)
    return random_number[0]



# O(r*n*k)
def gibbs_sampler(dna, k, repeats):
    #generate an initial random set of motifs from each string in dna
    motif_list = []
    for line in dna:
        random_index = random.randrange(0, len(dna[0]) - k + 1)
        motif_list.append(line[random_index: random_index + k])

    #repeat many times to converge on best motifs
    best_motifs = motif_list[:]
    for i in range(0, repeats):# r
        #remove a random string from dna and create a new profile matrix without it
        index = random.randrange(0, len(dna))
        motif_list.pop(index)
        removed_string = dna[index]
        profile_matrix = profile(motif_list)# O(t*k)

        #generate probabilities of all substrings of the string removed from dna using profile matrix
        probabilities = []
        for j in range(0, len(removed_string) - k + 1):# n
            pattern = removed_string[j:j+k]
            probabilities.append(string_probability(pattern, profile_matrix))# O(k)

        #generate a random substring of the removed string using the probability distribution, and add it to motif_list
        random_motif_index = random_with_probability_distribution(probabilities)# O(n)
        random_motif = removed_string[random_motif_index : random_motif_index + k]
        motif_list.insert(index, random_motif)
        if score(motif_list) < score(best_motifs):# O(t*k)
            best_motifs = motif_list[:]
    return best_motifs



# O(i*r*n*k)
def iterate_gibbs_sampler(dna, k, repeats, iterations):
    best_motifs = gibbs_sampler(dna, k, repeats)# O(r*n*k)
    best_score = score(best_motifs)
    for i in range(iterations):
        current_motifs = gibbs_sampler(dna, k, repeats)# O(r*n*k)
        current_score = score(current_motifs)# O(t*k)
        if current_score < best_score:
            best_score = current_score
            best_motifs = current_motifs
    return best_motifs