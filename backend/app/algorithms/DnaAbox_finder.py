# O(n)
def pattern_occurrences(text, pattern):
    occurrence_indices = []
    pattern_length = len(pattern)
    for i in range(0, len(text) - pattern_length + 1):
        if text[i : i + pattern_length] == pattern:
            occurrence_indices.append(i)
    return occurrence_indices



# O(n*k)
def approximate_pattern_occurrences(distance, pattern, text):
    distance = int(distance)
    occurrences = []
    pattern_length = len(pattern)
    for i in range(0,len(text) - pattern_length + 1):
        window = text[i : i + pattern_length]
        if hamming_distance(pattern, window) <= distance: # O(k)
            occurrences.append([i, text[i : i + pattern_length]])
    return occurrences



# O(n)
def generate_frequency_table(text, k):
    frequency_map = {}
    for i in range(0,len(text) - k + 1):
        kmer = text[i : i + k]
        if (kmer in frequency_map):
            frequency_map[kmer] += 1
        else:
            frequency_map[kmer] = 1
    return frequency_map



# O(n)
def frequent_words(text, k):
    frequency_map = generate_frequency_table(text, k) # O(n)
    frequent_patterns = []
    max_frequency = 0
    for kmer in frequency_map:
        if (frequency_map[kmer] > max_frequency):
            max_frequency = frequency_map[kmer]
    for kmer in frequency_map: #this might be redundant, add a variable max_kmer and just append that
        if (frequency_map[kmer] == max_frequency):
            frequent_patterns.append(kmer)
    return frequent_patterns



# O(n*k^2)
def frequent_words_with_mismatches(text, k, distance):
    frequency_map = {}
    frequent_patterns = []
    for i in range(0, len(text) - k + 1):
        pattern = text[i : i + k]
        neighborhood = neighbors(pattern, distance) # O(k^2)
        for neighbor in neighborhood:
            if neighbor in frequency_map:
                frequency_map[neighbor] += 1
            else:
                frequency_map[neighbor] = 1
    max_frequency = 0
    for kmer in frequency_map:
        if (frequency_map[kmer] > max_frequency):
            max_frequency = frequency_map[kmer]
    for kmer in frequency_map:
        if (frequency_map[kmer] == max_frequency):
            frequent_patterns.append(kmer)
    return frequent_patterns



# O(n*k^2)
def frequent_words_mismatch_reverse_complement(distance, k, text):
    distance = int(distance)
    k = int(k)
    frequency_map = {}
    frequent_patterns = []
    for i in range(0, len(text) - k + 1):
        pattern = text[i : i + k]
        neighborhood = neighbors(distance, pattern) # O(k^2)
        for neighbor in neighborhood:
            if neighbor in frequency_map:
                frequency_map[neighbor] += 1
                frequency_map[reverse_complement(neighbor)] += 1 # O(k)
            else:
                frequency_map[neighbor] = 1
                frequency_map[reverse_complement(neighbor)] = 1 # O(k)
    max_frequency = 0
    for kmer in frequency_map:
        if (frequency_map[kmer] > max_frequency):
            max_frequency = frequency_map[kmer]
    for kmer in frequency_map:
        if (frequency_map[kmer] == max_frequency):
            frequent_patterns.append(kmer)
    return frequent_patterns



# O(k)
def reverse_complement(pattern):
    complement_pattern = ""
    for char in pattern:
        match char:
            case "A":
                complement_char = "T"
            case "T":
                complement_char = "A"
            case "C":
                complement_char = "G"
            case "G":
                complement_char = "C"
            case _:
                print("invalid character")
                return
        complement_pattern = complement_char + complement_pattern
    return complement_pattern



# O(n^2)
def find_clumps(k, window_size, minimum_occurrences, genome):
    k = int(k)
    window_size = int(window_size)
    minimum_occurrences = int(minimum_occurrences)

    clump_patterns = []
    for i in range(0, len(genome) - window_size + 1):
        frequency_map = generate_frequency_table(genome[i:i + window_size], k) # O(n)
        for kmer in frequency_map:
            if (frequency_map[kmer] >= minimum_occurrences):
                clump_patterns.append((i, kmer))
    return list(set(clump_patterns))



# O(n)
def find_skew(genome):
    skew_array = []
    skew_amount = 0
    for char in genome:
        if char == "C":
            skew_amount -= 1
        elif char == "G":
            skew_amount += 1
        skew_array.append(skew_amount)
    return skew_array



# O(n)
def min_skew_indices(genome):
    skew_array = find_skew(genome) # O(n)
    min = skew_array[0]
    min_index_array = []
    for i in range(0, len(genome)):
        if (skew_array[i] < min):
            min = skew_array[i]
    for i in range(0, len(genome)):
        if (skew_array[i] == min):
            min_index_array.append(i)
    return min_index_array



# O(k)
def hamming_distance(str1, str2):
    if (len(str1) != len(str2)):
        raise Exception('length of strings are not equal')
    mismatches = 0
    for i in range(0, len(str1)):
        if (str1[i] != str2[i]):
            mismatches += 1
    return mismatches



# O(k^2)
def neighbors(distance, pattern):
    distance = int(distance)
    if distance == 0:
        return {pattern}
    if len(pattern) == 1:
        return {"A", "C", "G", "T"}
    neighborhood = set()
    suffix_neighbors = neighbors(distance, pattern[1:]) # T(k-1)
    for text in suffix_neighbors:
        if (hamming_distance(pattern[1:], text) < distance): # O(k)
            for x in {"A","C","G","T"}:
                neighborhood.add(x + text)
        else:
            neighborhood.add(pattern[0] + text)
    return neighborhood



# O(n)
def find_potential_dnaa_boxes(k, distance, window_size, window_location, genome):
    k = int(k)
    distance = int(distance)
    window_size = int(window_size)

    min_skew = min_skew_indices(genome) # O(n)
    first_index = min_skew[0]
    last_index = min_skew[len(min_skew) - 1]
    match (window_location):
        case "after":
            ori_region = genome[first_index : last_index + window_size]
        case "before":
            ori_region = genome[first_index - window_size: last_index]
        case "centered":
            ori_region = genome[first_index - int(window_size / 2): last_index + int(window_size / 2)]
        case _:
            raise Exception("invalid window location")
    return frequent_words_mismatch_reverse_complement(distance, k, ori_region)# O(s*k^2)