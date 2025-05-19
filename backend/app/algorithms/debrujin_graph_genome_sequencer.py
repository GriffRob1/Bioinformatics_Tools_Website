import random
import copy



# O(k)
def hamming_distance(str1, str2):
    if (len(str1) != len(str2)):
        return -1
    mismatches = 0
    for i in range(0, len(str1)):
        if (str1[i] != str2[i]):
            mismatches += 1
    return mismatches



# O(n)
def kmer_composition(k, text):
    kmers = []
    for i in range(len(text) - k + 1):
        kmers.append(text[i:i+k])
    return kmers



# O(n)
def path_to_genome(kmers):
    genome = ''
    for i in range(len(kmers) - 1):
        genome = genome + kmers[i][0]
    genome = genome + kmers[len(kmers) - 1][:]
    return genome



# O(n)
# translates a sequence of read pairs to a single genome string, checking for correct alignment.
def read_pairs_path_to_genome(all_pairs_string, d):
    genome = ''
    all_pairs_array = read_pairs_string_to_array(all_pairs_string)# O(n)
    num_pairs = len(all_pairs_array)
    k = len(all_pairs_array[0][0]) # here k represents the length of one side of the read pair, which is technically 1 less than the actual k for the data

    for i in range(k+d+1):
        genome += all_pairs_array[i][0][0]

    for i in range(k+d+1, num_pairs-1):
        if all_pairs_array[i][0][0] != all_pairs_array[i - (k+d+1)][1][0]:
            raise Exception('not a valid path')
        else:
            genome += all_pairs_array[i][0][0]

    last_pair = all_pairs_array[num_pairs - 1]
    #the corresponding pair, where the second pattern is the same as the first pattern in last_pair
    corresponding_pair = all_pairs_array[num_pairs - (2 + k + d)]
    if last_pair[0] != corresponding_pair[1]:
        raise Exception('not a valid path')
    else:
        genome += last_pair[0]

    for i in range(num_pairs - (d+1), num_pairs - 1):
        genome += all_pairs_array[i][1][0]
    genome += last_pair[1]
    return genome



# O(n)
# converts from string format to array format
# ex. 'AGGCT|GCTAG' -> [['AGGCT'], ['GCTAG']]
def read_pairs_string_to_array(all_pairs_string):
    all_pairs_array = []
    halfway_index = all_pairs_string[0].index('|')
    for pair in all_pairs_string:
        all_pairs_array.append([pair[:halfway_index], pair[halfway_index + 1:]])
    return all_pairs_array



# O(1)
def prefix(pattern):
    return pattern[:len(pattern) - 1]



# O(1)
def suffix(pattern):
    return pattern[1:]



# O(n^2)
def overlap_graph_from_text(k, text):
    kmers = kmer_composition(k, text)# O(n)
    graph = {}
    for i in range(len(kmers)):# n
        for j in range(i,len(kmers)):# n
            if suffix(kmers[i]) == prefix(kmers[j]):# O(1)
                add_adjacent_node(graph, kmers[i], kmers[j])# O(1)
    return graph



# O(n^2)
def overlap_graph_from_kmers(kmers):
    graph = {}
    for kmer1 in kmers:# n
        for kmer2 in kmers:# n
            if suffix(kmer1) == prefix(kmer2):# O(1)
                add_adjacent_node(graph, kmer1, kmer2)# O(1)
    return graph



# O(n)
def debrujin_graph_from_text(k, text):
    graph = {}
    for i in range(len(text) - k + 1):# n
        graph[text[i:i+k]] = []
    for i in range(len(text) - k):# n
        kmer = text[i:i+k+1]
        graph[prefix(kmer)].append(suffix(kmer))# O(1)
    return graph



# O(n)
def debrujin_graph_from_kmers(kmers):
    graph = {}
    for kmer in kmers:# n
        add_adjacent_node(graph, prefix(kmer), suffix(kmer))# O(1)
    return graph



# O(n)
def paired_debrujin_graph(all_pairs_string):
    graph = {}
    all_pairs_array = read_pairs_string_to_array(all_pairs_string)# O(n)
    for pair in all_pairs_array:  # n
        start_node = prefix(pair[0]) + '|' + prefix(pair[1])
        end_node = suffix(pair[0]) + '|' + suffix(pair[1])
        add_adjacent_node(graph, start_node, end_node)  # O(1)
    return graph



# O(1)
def add_adjacent_node(graph, node, adjacent_node):
    if node not in graph:
        graph[node] = {"indegree": 0,
                       "outdegree": 0,
                       "adjacent_nodes": []}
    if adjacent_node not in graph:
        graph[adjacent_node] = {"indegree": 0,
                                "outdegree": 0,
                                "adjacent_nodes": []}
    graph[node]["adjacent_nodes"].append(adjacent_node)
    graph[node]["outdegree"] += 1
    graph[adjacent_node]["indegree"] += 1



# O(n)
# removes the degrees from the graph object for visualization
# ex: {"AGGCT": {"indegree": 0, "outdegree": 1, "adjacent_nodes": ["GCTAG"]}}
#  -> {"AGGCT": ["GCTAG"]}
def remove_degrees(graph):
    graph_copy = copy.deepcopy(graph)
    for node in graph_copy:
        graph_copy[node] = graph_copy[node]["adjacent_nodes"]
    return graph_copy



# O(n)
# finds the start and end of the eulerian path
def find_all_possible_endpoints(graph):
    endpoints = [[],[]]
    for node in graph:
        if graph[node]["outdegree"] > graph[node]["indegree"]:
            endpoints[0].append(node)
        if graph[node]["indegree"] > graph[node]["outdegree"]:
            endpoints[1].append(node)
    return endpoints



# O(n)
# chooses a random start and end from the set of all possible start nodes and end nodes
def find_random_endpoints(graph):
    all_possible_endpoints = find_all_possible_endpoints(graph)# O(n)
    return [random.choice(all_possible_endpoints[0]),
            random.choice(all_possible_endpoints[1])]



# O(n)
# counts the number of edges in the graph
# PRECONDITION: graph is in no degree form
def find_edge_count(graph):
    count = 0
    for node in graph:
        count += len(graph[node])
    return count



# O(V+E)
# finds a eulerian path through the graph
def find_eulerian_path(graph):
    start_node, end_node = find_random_endpoints(graph)# O(n)
    graph_copy = copy.deepcopy(graph)
    add_adjacent_node(graph_copy, end_node, start_node)# O(1)
    return find_eulerian_cycle(remove_degrees(graph_copy), start_node)# O(V+E)



# O(V+E)
# finds a eulerian cycle through the graph
# PRECONDITION: graph is in no degree form
def find_eulerian_cycle(graph, start_node):
    path = [start_node]
    edge_count = find_edge_count(graph)
    temp_node = graph[start_node].pop()
    edge_count -= 1
    while temp_node != start_node:
        path.append(temp_node)
        temp_node = graph[temp_node].pop()
        edge_count -= 1

    while edge_count > 0:
        index = 0
        while len(graph[path[index]]) == 0:
            index += 1
        temp_start_node = path[index]
        temp_node = graph[temp_start_node].pop()
        edge_count -= 1
        while temp_node != temp_start_node:
            path.insert(index+1, temp_node)
            temp_node = graph[temp_node].pop()
            edge_count -= 1
            index += 1
        path.insert(index+1, temp_start_node)
    return path



# O(V+E)
def debrujin_assembler(kmers):
    debrujin_graph = debrujin_graph_from_kmers(kmers)# O(n)
    eulerian_path = find_eulerian_path(debrujin_graph)# O(V+E)
    return path_to_genome(eulerian_path)# O(n)



# O(V+E)
def debrujin_assembler_read_pairs(all_pairs, d):
    read_pair_graph = paired_debrujin_graph(all_pairs)# O(n)
    all_pairs_string = find_eulerian_path(read_pair_graph)# O(V+E)
    return read_pairs_path_to_genome(all_pairs_string, d)# O(n)