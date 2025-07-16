from graph import Graph
from motif_finder import hamming_distance
import numpy as np
import random

DNA_KEY = ['A', 'C', 'G', 'T']

DNA_MATRIX = [[1, -1, -1, -1],
              [-1, 1, -1, -1],
              [-1, -1, 1, -1],
              [-1, -1, -1, 1]]

BLOSUM62_KEY = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G','H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X']

BLOSUM62 = [[4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, -1, -1],
            [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, -2, 0, -1],
            [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 4, -3, 0, -1],
            [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, -3, 1, -1],
            [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1],
            [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, -2, 4, -1],
            [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, -3, 4, -1],
            [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -4, -2, -1],
            [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, -3, 0, -1],
            [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, 3, -3, -1],
            [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, 3, -3, -1],
            [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, -3, 1, -1],
            [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, 2, -1, -1],
            [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, 0, -3, -1],
            [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -3, -1, -1],
            [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, -2, 0, -1],
            [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, -1, -1],
            [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -2, -2, -1],
            [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -1, -2, -1],
            [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, 2, -2, -1],
            [-2, -1, 4, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, -3, 0, -1],
            [-1, -2, -3, -3, -1, -2, -3, -4, -3, 3, 3, -3, 2, 0, -3, -2, -1, -2, -1, 2, -3, 3, -3, -1],
            [-1, 0, 0, 1, -3, 4, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -2, -2, -2, 0, -3, 4, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]]

PAM250_KEY = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

PAM250 = [[2, -2, 0, 0, -2, 0, 0, 1, -1, -1, -2, -1, -1, -3, 1, 1, 1, -6, -3, 0, 0, 0, 0],
          [-2, 6, 0, -1, -4, 1, -1, -3, 2, -2, -3, 3, 0, -4, 0, 0, -1, 2, -4, -2, -1, 0, -1],
          [0, 0, 2, 2, -4, 1, 1, 0, 2, -2, -3, 1, -2, -3, 0, 1, 0, -4, -2, -2, 2, 1, 0],
          [0, -1, 2, 4, -5, 2, 3, 1, 1, -2, -4, 0, -3, -6, -1, 0, 0, -7, -4, -2, 3, 3, -1],
          [-2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3, 0, -2, -8, 0, -2, -4, -5, -3],
          [0, 1, 1, 2, -5, 4, 2, -1, 3, -2, -2, 1, -1, -5, 0, -1, -1, -5, -4, -2, 1, 3, -1],
          [0, -1, 1, 3, -5, 2, 4, 0, 1, -2, -3, 0, -2, -5, -1, 0, 0, -7, -4, -2, 3, 3, -1],
          [1, -3, 0, 1, -3, -1, 0, 5, -2, -3, -4, -2, -3, -5, 0, 1, 0, -7, -5, -1, 0, 0, -1],
          [-1, 2, 2, 1, -3, 3, 1, -2, 6, -2, -2, 0, -2, -2, 0, -1, -1, -3, 0, -2, 1, 2, -1],
          [-1, -2, -2, -2, -2, -2, -2, -3, -2, 5, 2, -2, 2, 1, -2, -1, 0, -5, -1, 4, -2, -2, -1],
          [-2, -3, -3, -4, -6, -2, -3, -4, -2, 2, 6, -3, 4, 2, -3, -3, -2, -2, -1, 2, -3, -3, -1],
          [-1, 3, 1, 0, -5, 1, 0, -2, 0, -2, -3, 5, 0, -5, -1, 0, 0, -3, -4, -2, 1, 0, -1],
          [-1, 0, -2, -3, -5, -1, -2, -3, -2, 2, 4, 0, 6, 0, -2, -2, -1, -4, -2, 2, -2, -2, -1],
          [-3, -4, -3, -6, -4, -5, -5, -5, -2, 1, 2, -5, 0, 9, -5, -3, -3, 0, 7, -1, -4, -5, -2],
          [1, 0, 0, -1, -3, 0, -1, 0, 0, -2, -3, -1, -2, -5, 6, 1, 0, -6, -5, -1, -1, 0, -1],
          [1, 0, 1, 0, 0, -1, 0, 1, -1, -1, -3, 0, -2, -3, 1, 2, 1, -2, -3, -1, 0, 0, 0],
          [1, -1, 0, 0, -2, -1, 0, 0, -1, 0, -2, 0, -1, -3, 0, 1, 3, -5, -3, 0, 0, -1, 0],
          [-6, 2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4, 0, -6, -2, -5, 17, 0, -6, -5, -6, -4],
          [-3, -4, -2, -4, 0, -4, -4, -5, 0, -1, -1, -4, -2, 7, -5, -3, -3, 0, 10, -2, -3, -4, -2],
          [0, -2, -2, -2, -2, -2, -2, -1, -2, 4, 2, -2, 2, -1, -1, -1, 0, -6, -2, 4, -2, -2, -1],
          [0, -1, 2, 3, -4, 1, 3, 0, 1, -2, -3, 1, -2, -4, -1, 0, 0, -5, -3, -2, 3, 2, -1],
          [0, 0, 1, 3, -5, 3, 3, 0, 2, -2, -3, 0, -2, -5, 0, 0, -1, -6, -4, -2, 2, 3, -1],
          [0, -1, 0, -1, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, 0, 0, -4, -2, -1, -1, -1, -1]]



#import numpy as np
#data = np.loadtxt('PAM250.txt', dtype=int).tolist()
#print(data)



def longest_path_to_alignment(longest_path, sequence1, sequence2):
    longest_path = id_to_tuple(longest_path)
    aligned_sequence1 = ''
    sequence1_index = 0
    aligned_sequence2 = ''
    sequence2_index = 0

    node = longest_path[0] # tuple of indices ex. ('0', '1')
    next_node = longest_path[1]
    # accounts for local alignment and sets the correct starting indices
    if next_node[0] > (node[0] + 1) or next_node[1] > (node[1] + 1):
        sequence1_index = next_node[0]
        sequence2_index = next_node[1]


    for i in range(len(longest_path) - 1):
        node = longest_path[i] # tuple of indices
        next_node = longest_path[i+1]
        if (node[0] + 1) == next_node[0] and (node[1] + 1) == next_node[1]: # match or mismatch
            aligned_sequence1 += sequence1[sequence1_index]
            aligned_sequence2 += sequence2[sequence2_index]
            sequence1_index += 1
            sequence2_index += 1
        elif (node[0] + 1) == next_node[0]: # deletion
            aligned_sequence1 += sequence1[sequence1_index]
            aligned_sequence2 += '-'
            sequence1_index += 1
        elif (node[1] + 1) == next_node[1]: # insertion
            aligned_sequence2 += sequence2[sequence2_index]
            aligned_sequence1 += '-'
            sequence2_index += 1
    return aligned_sequence1, aligned_sequence2



def id_to_tuple(ids):
    ids_as_tuples = []
    for id in ids:
        indexes = id.split('-')
        ids_as_tuples.append((int(indexes[-2]), int(indexes[-1])))
    return ids_as_tuples






def create_basic_global_alignment_graph(sequence1, sequence2):
    graph = Graph()
    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            graph.add_node(f'{i}-{j}')

    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            if i < len(sequence1):
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j}', weight=0)
            if j < len(sequence2):
                graph.add_edge(f'{i}-{j}', f'{i}-{j+1}', weight=0)
            if (i < len(sequence1)) and (j < len(sequence2)):
                weight = 1 if sequence1[i] == sequence2[j] else 0
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j+1}', weight=weight)

    return graph






def create_scored_global_alignment_graph(scoring_matrix, scoring_key, indel_score, sequence1, sequence2):
    graph = Graph()
    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            graph.add_node(f'{i}-{j}')

    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            if i < len(sequence1):
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j}', weight=indel_score) # deletion
            if j < len(sequence2):
                graph.add_edge(f'{i}-{j}', f'{i}-{j+1}', weight=indel_score) # insertion
            if (i < len(sequence1)) and (j < len(sequence2)):
                char1_index = scoring_key.index(sequence1[i])
                char2_index = scoring_key.index(sequence2[j])
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j+1}', weight=scoring_matrix[char1_index][char2_index]) # match or mismatch
    return graph






def create_affine_gap_global_alignment_graph(scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost, sequence1, sequence2):
    graph = Graph()
    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            graph.add_node(f'{i}-{j}')   # match/mismatch layer
            graph.add_node(f'D-{i}-{j}') # deletion layer
            graph.add_node(f'I-{i}-{j}') # insertion layer

    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            if i < len(sequence1):
                graph.add_edge(f'{i}-{j}', f'D-{i+1}-{j}', weight=initial_gap_cost) # initial deletion cost
                graph.add_edge(f'D-{i}-{j}', f'D-{i+1}-{j}', weight=additional_gap_cost) # additional deletion cost
            if j < len(sequence2):
                graph.add_edge(f'{i}-{j}', f'I-{i}-{j+1}', weight=initial_gap_cost) # initial insertion cost
                graph.add_edge(f'I-{i}-{j}', f'I-{i}-{j+1}', weight=additional_gap_cost) # additional insertion cost
            if (i < len(sequence1)) and (j < len(sequence2)):
                char1_index = scoring_key.index(sequence1[i])
                char2_index = scoring_key.index(sequence2[j])
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j+1}', weight=scoring_matrix[char1_index][char2_index]) # match or mismatch
            graph.add_edge(f'D-{i}-{j}', f'{i}-{j}', weight=0) # returning to match/mismatch layer costs nothing
            graph.add_edge(f'I-{i}-{j}', f'{i}-{j}', weight=0) # returning to match/mismatch layer costs nothing

    # make 0-0 a source node with no in edges
    graph.remove_edge('D-0-0', '0-0', 0)
    graph.remove_edge('I-0-0', '0-0', 0)
    return graph






def create_scored_local_alignment_graph(scoring_matrix, scoring_key, indel_score, sequence1, sequence2):
    graph = Graph()
    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            graph.add_node(f'{i}-{j}')

    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            # add 0 weight edge from 0-0 to every node
            graph.add_edge('0-0', f'{i}-{j}', weight=0)
            # add a 0 weight edge from every node to the sink node
            graph.add_edge(f'{i}-{j}', f'{len(sequence1)}-{len(sequence2)}', weight=0)
            if i < len(sequence1):
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j}', weight=indel_score) # deletion
            if j < len(sequence2):
                graph.add_edge(f'{i}-{j}', f'{i}-{j+1}', weight=indel_score) # insertion
            if (i < len(sequence1)) and (j < len(sequence2)):
                char1_index = scoring_key.index(sequence1[i])
                char2_index = scoring_key.index(sequence2[j])
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j+1}', weight=scoring_matrix[char1_index][char2_index]) # match or mismatch

    # remove self loops in source and sink nodes
    graph.remove_edge('0-0', '0-0', 0)
    graph.remove_edge(f'{len(sequence1)}-{len(sequence2)}', f'{len(sequence1)}-{len(sequence2)}', 0)
    return graph






def create_affine_gap_local_alignment_graph(scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost, sequence1, sequence2):
    graph = Graph()
    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            graph.add_node(f'{i}-{j}')   # match/mismatch layer
            graph.add_node(f'D-{i}-{j}') # deletion layer
            graph.add_node(f'I-{i}-{j}') # insertion layer

    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            # add 0 weight edge from 0-0 to every match/mismatch layer node
            graph.add_edge('0-0', f'{i}-{j}', weight=0)
            # add a 0 weight edge from every match/mismatch layer node to the sink node
            graph.add_edge(f'{i}-{j}', f'{len(sequence1)}-{len(sequence2)}', weight=0)
            if i < len(sequence1):
                graph.add_edge(f'{i}-{j}', f'D-{i+1}-{j}', weight=initial_gap_cost) # initial deletion cost
                graph.add_edge(f'D-{i}-{j}', f'D-{i+1}-{j}', weight=additional_gap_cost) # additional deletion cost
            if j < len(sequence2):
                    graph.add_edge(f'{i}-{j}', f'I-{i}-{j+1}', weight=initial_gap_cost) # initial insertion cost
                    graph.add_edge(f'I-{i}-{j}', f'I-{i}-{j+1}', weight=additional_gap_cost) # additional insertion cost
            if (i < len(sequence1)) and (j < len(sequence2)):
                char1_index = scoring_key.index(sequence1[i])
                char2_index = scoring_key.index(sequence2[j])
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j+1}', weight=scoring_matrix[char1_index][char2_index]) # match or mismatch
            graph.add_edge(f'D-{i}-{j}', f'{i}-{j}', weight=0) # returning to match/mismatch layer costs nothing
            graph.add_edge(f'I-{i}-{j}', f'{i}-{j}', weight=0) # returning to match/mismatch layer costs nothing

    # make 0-0 a source node with no in edges
    graph.remove_edge('D-0-0', '0-0', 0)
    graph.remove_edge('I-0-0', '0-0', 0)

    # remove self loops in source and sink nodes
    graph.remove_edge('0-0', '0-0', 0)
    graph.remove_edge(f'{len(sequence1)}-{len(sequence2)}', f'{len(sequence1)}-{len(sequence2)}', 0)
    return graph






def create_scored_fitting_alignment_graph(scoring_matrix, scoring_key, indel_score, container_sequence, fitting_sequence):
    graph = Graph()
    fitting_length = len(fitting_sequence)
    container_length = len(container_sequence)
    for i in range(container_length + 1):
        for j in range(fitting_length + 1):
            graph.add_node(f'{i}-{j}')

    for i in range(container_length + 1):
        for j in range(fitting_length + 1):
            if i < container_length:
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j}', weight=indel_score) # deletion
            if j < fitting_length:
                graph.add_edge(f'{i}-{j}', f'{i}-{j+1}', weight=indel_score) # insertion
            if (i < container_length) and (j < fitting_length):
                char1_index = scoring_key.index(container_sequence[i])
                char2_index = scoring_key.index(fitting_sequence[j])
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j+1}', weight=scoring_matrix[char1_index][char2_index]) # match or mismatch

    for i in range(container_length + 1):
        # add 0 weight edges from 0-0 to the first column nodes
        graph.add_edge('0-0', f'{i}-0', 0)
        # add 0 weight edges from the last column nodes to the sink node
        graph.add_edge(f'{i}-{fitting_length}', f'{container_length}-{fitting_length}', 0)

    # remove self loops in source and sink nodes
    graph.remove_edge('0-0', '0-0', 0)
    graph.remove_edge(f'{container_length}-{fitting_length}', f'{container_length}-{fitting_length}', 0)
    return graph






def create_affine_gap_fitting_alignment_graph(scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost, container_sequence, fitting_sequence):
    container_length = len(container_sequence)
    fitting_length = len(fitting_sequence)
    graph = Graph()
    for i in range(container_length + 1):
        for j in range(fitting_length + 1):
            graph.add_node(f'{i}-{j}')   # match/mismatch layer
            graph.add_node(f'D-{i}-{j}') # deletion layer
            graph.add_node(f'I-{i}-{j}') # insertion layer

    for i in range(container_length + 1):
        for j in range(fitting_length + 1):
            if i < container_length:
                graph.add_edge(f'{i}-{j}', f'D-{i+1}-{j}', weight=initial_gap_cost) # initial deletion cost
                graph.add_edge(f'D-{i}-{j}', f'D-{i+1}-{j}', weight=additional_gap_cost) # additional deletion cost
            if j < fitting_length:
                graph.add_edge(f'{i}-{j}', f'I-{i}-{j+1}', weight=initial_gap_cost) # initial insertion cost
                graph.add_edge(f'I-{i}-{j}', f'I-{i}-{j+1}', weight=additional_gap_cost) # additional insertion cost
            if (i < container_length) and (j < fitting_length):
                char1_index = scoring_key.index(container_sequence[i])
                char2_index = scoring_key.index(fitting_sequence[j])
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j+1}', weight=scoring_matrix[char1_index][char2_index]) # match or mismatch
            graph.add_edge(f'D-{i}-{j}', f'{i}-{j}', weight=0) # returning to match/mismatch layer costs nothing
            graph.add_edge(f'I-{i}-{j}', f'{i}-{j}', weight=0) # returning to match/mismatch layer costs nothing

    for i in range(container_length + 1):
        # add 0 weight edges from 0-0 to the first column nodes
        graph.add_edge('0-0', f'{i}-0', 0)
        # add 0 weight edges from the last column nodes to the sink node
        graph.add_edge(f'{i}-{fitting_length}', f'{container_length}-{fitting_length}', 0)

    # make 0-0 a source node with no in edges
    graph.remove_edge('D-0-0', '0-0', 0)
    graph.remove_edge('I-0-0', '0-0', 0)

    # remove self loops in source and sink nodes
    graph.remove_edge('0-0', '0-0', 0)
    graph.remove_edge(f'{container_length}-{fitting_length}', f'{container_length}-{fitting_length}', 0)
    return graph






def create_scored_overlap_alignment_graph(scoring_matrix, scoring_key, indel_score, sequence1, sequence2):
    graph = Graph()
    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            graph.add_node(f'{i}-{j}')

    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            if i < len(sequence1):
                weight = indel_score
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j}', weight=weight) # deletion
            if j < len(sequence2):
                weight = indel_score
                graph.add_edge(f'{i}-{j}', f'{i}-{j+1}', weight=weight) # insertion
            if (i < len(sequence1)) and (j < len(sequence2)):
                char1_index = scoring_key.index(sequence1[i])
                char2_index = scoring_key.index(sequence2[j])
                weight = scoring_matrix[char1_index][char2_index]
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j+1}', weight=weight) # match or mismatch

    # add 0 weight edges from 0-0 to the first column nodes
    for i in range(len(sequence1) + 1):
        graph.add_edge('0-0', f'{i}-0', 0)

    # add 0 weight edges from the last row to the sink node
    for i in range(len(sequence2) + 1):
        graph.add_edge(f'{len(sequence1)}-{i}',f'{len(sequence1)}-{len(sequence2)}', 0)

    # remove self loops in source and sink nodes
    graph.remove_edge('0-0', '0-0', 0)
    graph.remove_edge(f'{len(sequence1)}-{len(sequence2)}', f'{len(sequence1)}-{len(sequence2)}', 0)
    return graph






def create_affine_gap_overlap_alignment_graph(scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost, sequence1, sequence2):
    graph = Graph()
    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            graph.add_node(f'{i}-{j}')   # match/mismatch layer
            graph.add_node(f'D-{i}-{j}') # deletion layer
            graph.add_node(f'I-{i}-{j}') # insertion layer

    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            if i < len(sequence1):
                graph.add_edge(f'{i}-{j}', f'D-{i+1}-{j}', weight=initial_gap_cost) # initial deletion cost
                graph.add_edge(f'D-{i}-{j}', f'D-{i+1}-{j}', weight=additional_gap_cost) # additional deletion cost
            if j < len(sequence2):
                graph.add_edge(f'{i}-{j}', f'I-{i}-{j+1}', weight=initial_gap_cost) # initial insertion cost
                graph.add_edge(f'I-{i}-{j}', f'I-{i}-{j+1}', weight=additional_gap_cost) # additional insertion cost
            if (i < len(sequence1)) and (j < len(sequence2)):
                char1_index = scoring_key.index(sequence1[i])
                char2_index = scoring_key.index(sequence2[j])
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j+1}', weight=scoring_matrix[char1_index][char2_index]) # match or mismatch
            graph.add_edge(f'D-{i}-{j}', f'{i}-{j}', weight=0) # returning to match/mismatch layer costs nothing
            graph.add_edge(f'I-{i}-{j}', f'{i}-{j}', weight=0) # returning to match/mismatch layer costs nothing

    # add 0 weight edges from 0-0 to the first column nodes
    for i in range(len(sequence1) + 1):
        graph.add_edge('0-0', f'{i}-0', 0)

    # add 0 weight edges from the last row to the sink node
    for i in range(len(sequence2) + 1):
        graph.add_edge(f'{len(sequence1)}-{i}',f'{len(sequence1)}-{len(sequence2)}', 0)

    # make 0-0 a source node with no in edges
    graph.remove_edge('D-0-0', '0-0', 0)
    graph.remove_edge('I-0-0', '0-0', 0)
    # remove self loops in source and sink nodes
    graph.remove_edge('0-0', '0-0', 0)
    graph.remove_edge(f'{len(sequence1)}-{len(sequence2)}', f'{len(sequence1)}-{len(sequence2)}', 0)
    return graph






def basic_pairwise_global_alignment(sequence1, sequence2):
    graph = create_basic_global_alignment_graph(sequence1, sequence2)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(longest_path, sequence1, sequence2), score



def scored_pairwise_global_alignment(scoring_matrix, scoring_key, indel_score, sequence1, sequence2):
    graph = create_scored_global_alignment_graph(scoring_matrix, scoring_key, indel_score, sequence1, sequence2)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(longest_path, sequence1, sequence2), score



def scored_pairwise_global_alignment_affine_gap_penalty(scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost, sequence1, sequence2):
    graph = create_affine_gap_global_alignment_graph(scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost, sequence1, sequence2)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(longest_path, sequence1, sequence2), score



def scored_pairwise_local_alignment(scoring_matrix, scoring_key, indel_score, sequence1, sequence2):
    graph = create_scored_local_alignment_graph(scoring_matrix, scoring_key, indel_score, sequence1, sequence2)
    longest_path,  max_score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(longest_path, sequence1, sequence2), max_score



def scored_pairwise_local_alignment_affine_gap_penalty(scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost, sequence1, sequence2):
    graph = create_affine_gap_local_alignment_graph(scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost, sequence1, sequence2)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(longest_path, sequence1, sequence2), score



def scored_pairwise_fitting_alignment(scoring_matrix, scoring_key, indel_score, container_sequence, fitting_sequence):
    graph = create_scored_fitting_alignment_graph(scoring_matrix, scoring_key, indel_score, container_sequence, fitting_sequence)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(container_sequence)}-{len(fitting_sequence)}')
    return longest_path_to_alignment(longest_path, container_sequence, fitting_sequence), score



def scored_pairwise_fitting_alignment_affine_gap_penalty(scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost, container_sequence, fitting_sequence):
    graph = create_affine_gap_fitting_alignment_graph(scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost, container_sequence, fitting_sequence)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(container_sequence)}-{len(fitting_sequence)}')
    return longest_path_to_alignment(longest_path, container_sequence, fitting_sequence), score



def scored_pairwise_overlap_alignment(scoring_matrix, scoring_key, indel_score, sequence1, sequence2):
    graph = create_scored_overlap_alignment_graph(scoring_matrix, scoring_key, indel_score, sequence1, sequence2)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(longest_path, sequence1, sequence2), score



def scored_pairwise_overlap_alignment_affine_gap_penalty(scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost, sequence1, sequence2):
    graph = create_affine_gap_overlap_alignment_graph(scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost, sequence1, sequence2)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(longest_path, sequence1, sequence2), score






def consensus_sequence(sequences, character_list):
    # creates a matrix that counts character frequencies by position
    count_matrix = np.zeros((len(character_list), len(sequences[0])))
    for i in range(len(sequences)):
        for j in range(len(sequences[0])):
            char = sequences[i][j]
            if char != '-':
                index = character_list.index(char)
                count_matrix[index][j] += 1

    # creates a consensus string based on the most frequent character from each column
    consensus_pattern = ''
    for j in range(0, len(count_matrix[0])):
        column = count_matrix[:, j]
        max_value = float('-inf')
        for i in range(len(column)):
            if column[i] > max_value:
                max_value = column[i]
        max_indices = []
        for i in range(len(column)):
            if column[i] == max_value:
                max_indices.append(i)
        index = random.choices(max_indices)
        consensus_pattern += character_list[index[0]]

    return consensus_pattern






def edit_distance(sequence1, sequence2):
    sequence1_length = len(sequence1)
    sequence2_length = len(sequence2)
    distance_array = np.zeros((sequence1_length + 1, sequence2_length + 1))

    for i in range(sequence1_length + 1):
        distance_array[i][0] = i
    for j in range(sequence2_length + 1):
        distance_array[0][j] = j

    for i in range(1, sequence1_length + 1):
        for j in range(1, sequence2_length + 1):
            if sequence1[i - 1] == sequence2[j - 1]:
                distance_array[i][j] = distance_array[i - 1][j - 1]
            else:
                distance_array[i][j] = 1 + min(int(distance_array[i][j - 1]), int(distance_array[i - 1][j]), int(distance_array[i - 1][j - 1]))

    return int(distance_array[sequence1_length][sequence2_length])






def greedy_multiple_alignment(scoring_matrix, scoring_key, indel_score, sequences):
    sequences = sequences.copy()
    aligned_sequences = []

    least_distance = edit_distance(sequences[0], sequences[1])
    closest_pair = (sequences[0], sequences[1])
    for i in range(len(sequences) - 1):
        for j in range(i+1, len(sequences)):
            distance = edit_distance(sequences[i], sequences[j])
            if distance < least_distance:
                least_distance = distance
                closest_pair = (sequences[i], sequences[j])

    pairwise_alignment, score = scored_pairwise_global_alignment(scoring_matrix, scoring_key, indel_score, closest_pair[0], closest_pair[1])
    aligned_sequences.append(pairwise_alignment[0])
    aligned_sequences.append(pairwise_alignment[1])
    sequences.remove(closest_pair[0])
    sequences.remove(closest_pair[1])
    consensus = consensus_sequence(aligned_sequences, scoring_key)

    while len(sequences) > 0:
        least_distance = edit_distance(consensus, sequences[0])
        closest_sequence = sequences[0]
        for j in range(len(sequences)):
            distance = edit_distance(consensus, sequences[j])
            if distance < least_distance:
                least_distance = distance
                closest_sequence = sequences[j]

        pairwise_alignment, score = scored_pairwise_global_alignment(scoring_matrix, scoring_key, indel_score, consensus, closest_sequence)
        aligned_sequences.append(pairwise_alignment[1])
        sequences.remove(closest_sequence)
        consensus = consensus_sequence(aligned_sequences, scoring_key)

    return aligned_sequences










seq1 = 'CAAGGAACTCAATCGTCGAGTCCACGGGGGGCAGAACACGCTATATTTAATCTTGATGAGGAACGCAAATAACCATGGTTGCACGTGAGGATTTTCTTTAGTGAGTTGGGTTGCTTGGTAACTTATCCACTGCTATCTTAAGGGGGTTACTTCGGGATGAACGGCTTATGACAATCACAGTGAGGTCCGTCCCGGCCGATATGAGTTCTATGTTTTAACAGCGTCACCAGTGTCACGTACGGGGCCACCTCAGGCCCTGACCAGGGAATAGAGCGATTTGGGGACTTTCCCGGGTGATGTCTACCAGGAAGTTCGGTACCACTGACTTTGAATAATACTGTCAAAGGGGCTGCACCTTCCCGAGTTCGTCGTCATTACACAGCGCATATATTACACGTTAAGCCGTTTATCCGCATGTTATGCCAATTCGCGTCTTGCCAGGTGCCAACGAGCCTGATAAAGCAGTGGGTAGCGCCGGCACAGTATGTAGCAAGTTCCCCGCCGCGCGTTGAAAGCGTTACGTACAGGCGGCTAAGCGACGTTAAAATTGTCGCTTGCCTAACCCATCTCCCTGACACGGAACATAGCGAATAGTAGTCAACGGAGTTATGGTACAAAGCCTGAAAGCGACCTCAGACGAAGGGTCTGCCCGCAGGACGTGGGCTCTAATCCTCGGGGGCCTCGCCTACGTAGCACATCCCCAATAGCACTAAGAAGATGTGAACGAAACGCCGCTGTCGGATTCCAATTCTGAAATAGATAGTACCGGGTCCGAGGCGATGGAGGGTGGCGAAACCCCCATTTACGCATAGCGGTAACTTGGTCCCGGACTATTTATCAGTTGGTACCCTCGGCCCTGGTGGATGTGTTTTACGATGCTATAGCGCGTATCGATTAGCTATGCTATCTATATTGCGCGCATATGCTAGGCTATGCTAGCTCTAGAGCAGCACACATATTCGATCGTATACGTACGTACGTACGTACGTCGTCGATGCTAGCTATCGATCGACTAGCTGGAATGATGTATAGACATCGCCTAC'
seq2 = 'CGCCCTACCTTTCGTGCCTATCAAAGATCTAGCTAGCTAGCGACGTAGCTATTATCTATCTCGAGCTACTATCGTACGAGCCGCGAGCTCTTCATCGTATATCGGCTATCGACTAGCAGCTAGCTAGCAGCCAGATCACTATTTCAGTTACGCATCGGTAGGCGCATTCCGAGCTCCGACTAGCGAGAAAACCAGCACAATGAAGTGCGCCTCCATATGTTAGTATACGGGTTGCCTCAGCTTTGGGCGGGCCTAAGGGCGGGATGAACGGCTTATTCCTGTGCAGTAGGTCGGTCCCGCCGATATGAGTTCTGGAGGTATTTAACAGCAGGGGCACGTGCACTGGGCCGCCTCAGGGGCTCACCATGACCAGGACATCGAGCGGGATCTCAACTTTGGGGACTTTCCCGGGTGATGTCTGCCACTGATGTTCGGTACCAGCGACCATACTGTGAAAAGGGGCTGCAAACATTATCTTCCCGAGTTGTCATTACACAGCATATATTACGCCGAAACTACGTCTCCGCATGTTTAACTCGCGTTTTGCCAGAGCCTAGCTAGGTTATAAAGGCGGCACAGTTTGTTCTAAAACCCGCCCTCGCCGAAGCGTTACGTACAGCGGCAAAGCGACGTTGAAATTGTCGATTGCCTATGGTTGACGCGGAACATAGCGTGACATAGTAGTGATACTCCAAGCGTGAAAGTGACCTAGCGGGTCTGAGCAGGACGTTGGCTCGAGCACGTGCAAGTTAGGGTTGTTCATGCCTGGGGAAAATAGTTTACGATTATTCGGCCATGCGTTAGCGTGCTGCTGTCCAGGCGCAGGGGCTTTCGAAAGTTTACTCCGTGATTGCGATTACCCTGAAGCCAGAGCTACACTTCCACGGACCGAACGCCGAATATCCGGATCCTGCCTGGCTTTCCGGCTTCGGAACTGAGAGAGGATCCGAGAGATAGCTGGTAA'

'''alignment3, max_score1 = scored_pairwise_global_alignment_affine_gap_penalty(DNA_MATRIX, DNA_KEY, -10, -5, seq1, seq2)
print(max_score1)
print(alignment3[0])
print(alignment3[1])
print(hamming_distance(alignment3[0], alignment3[1]))'''



'''seq3 = 'ENLHDMCNDHMIHGAVFKKQYKSGEPKEQIFWYPMIQFLISVCRSIVRNQKDADWPPHKHWERCYHNSTPTCCLFYGQHHHFDVAGMHKIAPMHLLVDGKTQMTILLRLQDHGAHWPYCIVTTRRAMGTEAEKSDYFNKTCKLRAVSWLYDQYLEPFNGYYSVCKNCGADRGECLMRYIAHPKMFCYYWNPNLRPPGHLMGNPRIWTMPAELFSYTFCTYYLDKLRWIQCNAMDWIHSYWETVDTNCMYITMFGLMQDNTVNTACCIPFIMGFIFSDMSNQGNSVEGMWPMMTNLQRARCNFHRSPKGEKYFTWGYDKAHGIGQAYPAKMMNQSMDQQTCRGAWYPLFGIWAVASFHSCLWESKKPHLAVDCAENMCAYHHMWLYSPGDGHVDARLDYTLTGKSKMPCWPKKQWADACDVRVCYAIFIKWLICNAEFGALMYQMCDKPVCFLMACKIQYTWMARCWADMWWNFRIEAQLDYDWVWRMLTNMMTPLCVHANYNQKYDEDSWCDSANSQFAWAWYVQWFNEWAQLTYMSEMQVWSFLSTGPIEIWMPTYEMGDASTNTDNLRWNYMFVLKVTHYVHFNDVHPYSRHNDPSLHIMWDMGSHCKFNVDMDYMIWMNHCRFCKLGTINYNNAPNEDTFNTRDWMVFSPYHSAVCEMYNWNGRPDGLVIDEAAPQSSMGINSDNRSSWNHRCNRYQCFDWYYMKDEYHERHNSSTFFYWGPIDIFQFKMPYCPRDKLGMGKWKQDNPDVVGCWGIGGVFFKCYKADLMVRSPAYHYGGCPHGVIQEPEVHRCRKCFLITHAFRGVINDWIFSMTQLDSSERHINIIDKKSHATCVYARAPLLPHVQFAVKTCWKWRDTVDVVAWQQYLITAAKGLAKIIMICWIYHLHEIHDQCPVYCGPIMPGNWQFNIEGQCCIMACKDTNQATAWHQLGSWGKHTEVRLAISNPHIHFSGSIVRLSVHCHHVSYHSIRALFETNMSCMSELQVTNLFCQFGRAANYKQSQEGKKCHNPVQQWHVLMWAMFTNEHPDEMKCAVMDFFPYPILNQNFADERHMQHFINVVTDYHQGEQYGYIVIEETEARIIVGCKSHDRLDDDLSSHVNHLYYLDLIRMDVMGDGMCMPYCPYCSQRSPPQISPMVHMGCFVWITFICKVQYIVCWIKTGPSLRYETYPWDSAKWLSWFAQINAWLVVIPCDFLMDTNNSCDHTIQLFALINMRFLSYLEWRMPAAIWAARFNCPETPFIPMYACFEMSCWLGNFPFVPHQDTTIYNTFNFSHHTFMKRKGHGVDCTKWTTREQAKCATASEENQTEEMNFDNYSFRIMRANFSLRFIFSYIERKLTSMNKTWWIGWSCRNPDQHSWLQMYFTHILCAEFEAFTAPPHYLHPMGEVMHIIHVHPRPPGMMDYLHHHMPQPFHQGDGDSIQSVIKSAAPYCFCQYGKIHYWACHIYDHDFAAMNTYITTEKAFHNYKYSQSPIMDIRGGPASTCSHRSNMAWNVRFGSNRKIYGYACMTCFSVIMMKMHRQKWKHGNCMQGAPRVCRLVAAGVFECRYHNDYQDVCTCNNRWTWEWHMMRHKDALLMSMPTTAIQPERSVISGEGENEVFAAKNQMVPQCPPTWETYRRNVNVQRSSCPRWTKVQTPHLSILVYNNRKQFMTNSSWMFAGTPTIDQMDDEPNYYKEWGNCTYAAMSQQKKVQIEFACCINHDNRFQLVAFHRAQAFNNIVYYVWFAKKRVQCPALDEHGAKMMMFYIIASYLFPNQIDFEWMKSPMEPCNMMDRKGRIRAQSKHFWMGRFYQPKAETARCSKWVNWMTIVAYWGSDHQVAKTPVRHWHSIRTVPMYVKFYCNIQPLSPEVAAPWNNIGMLFGNEHEMVARNALDKQDAPAFVAKFGHDYHAKEVFEKWAGTQPKTTNAGFISYWCDCYKRMAGCMSIWWGPVHRVPKFDWTCNLLCEWCRSMKQVDQTAIAYMTIMWCIRHFHTSKGWTHRMPAWEAIMCYVCELCDNLRERMRPFSLRLPCYEATCSDPCNEIWYFIAHPAAYEDQKTSVLKRGEDHVCSNVTHNKMKVQLSGVRHFEHYRTMQYKTRSCFKMCMRKMHFEMEQMDMIMWGVPDTYNMATWQVPHQMMSEHGHNDWRNSQHTPYHWIGRVFSGEWDQHNHAIDRRMQPNGRPDIPLSQRTEPDSLGSPAIMTLFKLISNKTNWTSQSGVNESYQFDHKNEDFNVEKSYPDIYNMQFFMQVLIRKIHKISFLEDYAVSPDWQMRKWSSHQQNNLRLPVRDRAADPSSSDFGFLDMPYDENMPNYNHYLFRQWADWMQPMEDEDSYRMQMYQNIVITEFPSFYIMCFHQFWSCVMWERTYGQLHLKNMFSRFPYGLQDHGFHPMHNSGNGPGLWARKRTYLGLQACPISVFYVRPIAPRWEGRNMTTTSHMWQAGQCWKTHPEPQPPLKGNDCSKEEWRTVFISGEEHQWWCKWRFRRPSMHSPCKYKLCWRLLDTRSKQSAPNLKGLAYNGKFLWDYEYHGLHNLEDDDFSLVIMQDRFMWQWWRRAGVMYDNTLPHFYRTWSSARLWHIVWACCMCQKWLPHNGRLGMWQMCFYNMPMMDNYEDAALLTGTPHLPIHMMMEIHFGWCGRLWHFGMQLVRMLPDLNWALYNMLCVTEWTMGHSHDHGWMNDRQVPRQCGWCTDQVMTPFLYFLNARKVHVLSCCD'
seq4 = 'LNQQIVLSGHKFRMKQELLSNNNSQMETMKPHEGACQDCLALERSNSTQCPQGWAIVWYDKLKMMMFAQDLVYEDMMSRRTNLMAWFQHSVCIWNIALNVGRNDEDSMDRSEHNVPSETQIQTMQLYMNHADLFWWPATFREWHGEKMYKTYDTWLWTIEINNCHMFAWVCSECPWWCDDWNMWDYADSYNIYVASFAFRNFGKGSAVDVMSWVPFMTIIKLRVTIVLYIQYTHHMTNVISHYTWMWMWRIMVCLMETQLIINIWFQQWVCSEACVCAAHQMMLTMNVHTVLQIVFMYMPCWPRVRPPMESFQDCRSGLCPRQTSFICRLADWRCRCVCPFHRTKQANCICCYDQSAITGPYRIIVATCGCMIKLAWVYNQESCANCTRRQWPGMAIETIIGKDFDANQVQIKRLLSNAPFYVIQEYNKAIHNYYVMCMFHLSAWWCEIIRTIQPPEDEHGMHNFWFGYPFSGNFEHLNEGAQKMVVQAPHYVYNRKWHHMYCWHYDYSLCRVAFRQKIKECQKVVWWLRPAITWAMETEQPWEYNCFKNFQNWLVWCTYPTHYHATCCYDHSCSTDVNSCNLITKFDCCINWMEIKHFTINRHPWIDNFIETHFPVHPERGMMLPCRESKDHRILERQSQQERDRHLFPMSLGCPNRKGLFYWCEKIPHEPCHGFQYNRNLAHTQWGYQVTGSNMQAHCTDWNWNAWQKHVIRIGMMNVHWMLATSVGRQRVDCHARLATTWYKTQWKWKEKEMWCMACHKTHNFYYDNQEFASPLWHCPWWPIESWWWNPKHGCNTCRGVGFRKEKQMKGRPHTGCHFQCMSMVSSERHNNIIDKKSHATCVYARAPLLSMHHFGYKHVQFAVKTCAPYTDRHQKWQCYDVDTTDVVAWQQYGITAAKGLAKIIISCPICWIYHLHYCRFVEHHDQCPVRDFLGQWGVSEGFQHIMPGNVRTSQFNIEGQMCIMACRDRNQATAWHQLGSWGKHIEVRLAVRLSVHCHNKIMFILVSYHSIRALFETNMSCMSELSVTNLFCQFGRAANPPSQEGKKCHNPVQDWHVLMWAMFTEESPDEMWMAMREDFVPYYMILERHMQHFINVVDKDYHQGQYGCIVIENIAARYCKSHDRLDDDEIYPWWYFSSHVNHLYEPGCIMLDLIRIPKLCPLDVMGDGMCMPYCPYCSQRNPPQISPMVHCFVWITFILKVQYFVCWIKTGPSLRKEKYAKWLSVFAQIIPCDFFRMDTNQALINMRQLSELTYRMPAAIWAQRFWCPETPFPMYASFEMKCWLGNFPFVPHQDTTIHHTFMVHPMRWPNYGHNNYWRSQAKCTPHRTASEENQTEAMNFDNYSKRIMRANFSLRFIFLTSMNKTWWDQHTWLQDMTFDVFQHILCAEFEAFTANTMSPPHYLHPMGEVMHIIHYPRPPGMMDVLHHHMPQPSWKKFWDRCNYCFCQYGKIPHWACHIGDHWFFAGNTYINTEKAFHNSKYSQSPIMYNIPFASTCHKMTGSHRSNMAWVVRFGKNRKIYGYTCMTCFSVITMKMHKRGCTLGFWWKHGNCMQSATAHIETWHCRLVAAFVFECRYHNYIGVMGYQDVCTCNNRWTWEQPTCHMMRHKQIHMWPDEASDFQKERMTNKVPRCPPWWPPTGCAVITYRRNVNHMQCQRSSCPKVQTPQDSILVYNNRKQLANAMFAGTPTPEEHYDRFDQMDDEPNYRHKERGNCTYPSKMGSANAMRQQGDCIRHDNEFQLVAFHRAQAFNNNVYNPWDNFQCPALDEHGAEMMMFYIIFSFQEMQYVEHAMCCNMRAQSKHFWMGRFYQPKAETRRCSKWVNWMTIVSNYWGSDHQVAKTPARTVYCDAQPLSMEECAPYVTLFYNEHEMVNALQKQDFGHDYHAKEVFDYIFLMKIGVTQLLPKDVTIISRPNNKAASDSIHKFTAKGLYVMMHTEAEGKFYEDKNQGWEIGFHIWGYQQFADSVVSVPTSISNWTVIAPQQKYYKVMRIIVSHKRDTSLILPNIMWIVWVNIFEDSMNCTLNHHNNPKAYAVPMGRIDYTIPIDIKPNVFFIMNCECEFDLEHHWCKALASCSDFLYADSLEGTGLGHRKSIVCCMYDYGECDCAEGPHTLGSYKRTAPHYTQVFERFMYPINMYMLTEAMSTGTVGHAIPDHDGECVQGSFDYNEIVMEHHIHRWRILDTKWKIRSALALKIGMYRWCCSKMCGFKPLHLGPPGSPPQYANPMEQYWLAWYYPCPYKQLLESYWMDANSQETPQVYNMGAEPTRFSIWFQNPYMENTWVHSHSERPCGYQAPYFKANCLRTCGSLWGFSTYPIYWGQSFYIQREMNQTWSNLWMTVMLLRYKGDTAECMVFCVWNHDEIFQTNKLQGNMALFMNGIWDTIYTELGHWWIHWPKGPPMCLPEPNMAWIMGSVEPSERQAGNKLTIQAYTFKSASYACHLVLYCRLCDHPASYRAFSVILMMSEDYYCQKDHKKFNEMWNIEMQNESRINVGPTPANNEFCADEMETELFCVGRDYAVNGDQQIGCVKYYHQNLTIAMVREANRECFQSIDWSRYWCRVGGQKLSPVKEWWIKVYEGQPCDMFCRTPEKEVHLFAQNPTTKWHVCSYLKDLARIVINSAMNASIAKDSSCIARQAKKPMPTYHPGEIPMGVPGDQGNYAFHKKERTHHEPWMENACFASNWCM'
alignment4, max_score2 = scored_pairwise_global_alignment_affine_gap_penalty(BLOSUM62, BLOSUM62_KEY, -10, -2, seq3,seq4)
print(max_score2)
print(alignment4[0])
print(alignment4[1])
print(hamming_distance(alignment4[0], alignment4[1]))
'''


'''container_sequence1 = 'KQHFYGGYSNNRVEKQNTPGWSTYTCENCQMCGVLMTNEQLLSMCPMQQYMGMRIKQVWFTWVRMWPYHESCMRHQGYQFVNHVNWVPYSKTIHEVNKFWDFMSFDCIKAYHFLYHSWRMHLWCFIGSPYIRCVICYTCPGWANWQHHWYNHEYSEWTMCQSGYNYVSCMMSDDMAHDKFIFHPHFGSGFIYHGAQEDSMRLQHMIRWSSFWEEIPVNVIIIEHYTIGVEHMGFCYCVNTCFDMICHWRNAGRLELCEQWLRGFKWPQHMKDEVSQSKEMIDLMSEMGWMHHETTQYHSKGFREFYIIQRWCSCQEVEPQWFPKMYPGQRPFQTPVCPLNYDPKRVNEPTPLHFCCQEFPVANWATMERKQTYVWDKNQSAVWMFGNLETIFTPIITPISETSAYNGYGSDMCSVIRAVAWCEVNVTNRMQDLRASWDMMLCGMEHSRIGHNTVDYLKWMQSYNTWPAFFNRPFLKYSPWKHMNHHSQWEYHHDFEYYYLGRANHYHRIESTNQEVWYGDEMMKGGVREWKDLPCGFKHGNVQDNLAKKHKWSSMVYSGEEEMLFTWQFGGSEKHKKWFYSAVWPKPWIDWRQLHGKEAFRHWRFYYYWVNPDFIYFQCVAGTPVYTMFNPGHICCNEFVCVYPICLQAGLYVETKGCIMQQQDRVKMSQNRDTYCPMGHHKWVPDDLNMNLAQWREANEQAVIDWINEHDEMDFKISAMYWWSKIFASQKMRIERKMQCLLFWPEPRVRVLFTQMYLMETHNHVGLKWHHYPLFSVKAKKLFVISYQMQHQLRLQNFAHAVGFVTSPLINAVDNNRNTWMTHMWLHGVKTKIRSGFVLIRASDFCVNYHLHDPMTMPLPDEWIKVTQLKDTMLQNHTAPYNYRKMACLNDDACPGDNIRAKPHPGETAHAHSYHSELREMQMIQGAWDFAANGMCSVYHVVTKSPRHTLPFVHCWQMKNGEWQGFWRWWSGITISMTHIMVPCIGLWLTTLIQDGNGDY'
fitting_sequence1 = 'FGFKHGNVQDSSMVYSEMLFTWQFGGNEKHKKWFGSWVWPKPWIDGPWRQLHGKEAFRHWYYWVNPDFIYHAMRVQVFRFQCVAGTPVYTMFNPGHIC'
alignment5, max_score_3 = scored_pairwise_fitting_alignment_affine_gap_penalty(container_sequence1, fitting_sequence1, BLOSUM62, BLOSUM62_KEY, -10, -5)
print(max_score_3)
print(alignment5[0])
print(alignment5[1])'''

'''seq5 = 'CAAGGAACTCAATCGTCGAGTCCACGGGGGGCAGAACACGCTATATTTAATCTTGATGAGGAACGCAAATAACCATGGTTGCACGTGAGGATTTTCTTTAGTGAGTTGGGTTGCTTGGTAACTTATCCACTGCTATCTTAAGGGGGTTACTTCGGGATGAACGGCTTATGACAATCACAGTGAGGTCCGTCCCGGCCGATATGAGTTCTATGTTTTAACAGCGTCACCAGTGTCACGTACGGGGCCACCTCAGGCCCTGACCAGGGAATAGAGCGATTTGGGGACTTTCCCGGGTGATGTCTACCAGGAAGTTCGGTACCACTGACTTTGAATAATACTGTCAAAGGGGCTGCACCTTCCCGAGTTCGTCGTCATTACACAGCGCATATATTACACGTTAAGCCGTTTATCCGCATGTTATGCCAATTCGCGTCTTGCCAGGTGCCAACGAGCCTGATAAAGCAGTGGGTAGCGCCGGCACAGTATGTAGCAAGTTCCCCGCCGCGCGTTGAAAGCGTTACGTACAGGCGGCTAAGCGACGTTAAAATTGTCGCTTGCCTAACCCATCTCCCTGACACGGAACATAGCGAATAGTAGTCAACGGAGTTATGGTACAAAGCCTGAAAGCGACCTCAGACGAAGGGTCTGCCCGCAGGACGTGGGCTCTAATCCTCGGGGGCCTCGCCTACGTAGCACATCCCCAATAGCACTAAGAAGATGTGAACGAAACGCCGCTGTCGGATTCCAATTCTGAAATAGATAGTACCGGGTCCGAGGCGATGGAGGGTGGCGAAACCCCCATTTACGCATAGCGGTAACTTGGTCCCGGACTATTTATCAGTTGGTACCCTCGGCCCTGGTGGATGTGTTTTACGATGCTATAGCGCGTATCGATTAGCTATGCTATCTATATTGCGCGCATATGCTAGGCTATGCTAGCTCTAGAGCAGCACACATATTCGATCGTATACGTACGTACGTACGTACGTCGTCGATGCTAGCTATCGATCGACTAGCTGGAATGATGTATAGACATCGCCTAC'
seq6 = 'CGCCCTACCTT----TCGTGCCTATCAAAGATCTAGCTAGCTAGCGACGTAGCTATTATCTATCTCGAGCTACTATCGTACGAGCCGCGAGCTCTTCATCGTATATCGGCTATCG-----ACTAGCAGCTAGCTAGC--AGCCAGATCACTATTTCAGTTACGCATCGGTAGGCGCATTCCGAGCTCCGACTAGCGAGAAAACCAGCACAA---TGAAGTGCGCCTCCATATGTTAGTATACGGGTTGCCTCAGCTTTGGGCGGGCCTAAGGGCGGGATGAACGGCTTAT----TCCTGTGCAGTAGGTCGGTCCCGCCGATATGAGTTCTGGAGGTATTTAACAGCAGGGGCACGTGCACTGGGCCGCCTCAGGGGCTCACCATGACCAGGACA--TCGAGCGGGATCTCAACTTTGGGGACTTTCCCGGGTGATGTC---TGCCACTGATGTTCGGTAC-CAGCGACCATACTGTGAAAAGGGGCTGCAAACATTATCTTCCCGAGTTGTCA---TTACACAGCATATATTACGCCGAAACTACGTCTCCGCATGTTTAAC---TCGCGTTTTGCCAGAGCCTAGC---TAGGTTATAAAGGCGGCACAGTTTGTT--CTAAAACCCGCCCTCGCCGAAGCGT-TACGTACAGCGGCAAAGCGACGTTGAAATTGTCGATTGCCTATGGTTGACGCGGAACATAGCGTGACATAGTAGTGAT-----------------ACTCCAAGCGTGAAAGTGACC-TAGCGGGTCTGAGCAG---GACGTTGGCTCGAGCACGTGCAAGTTAGGGTTGTTCATGCCTGGGGAAAATAGTTTACGATTAT----TCGGCCATGCG---------TTAGCGTGCTGCTGTCCAGGCGCAGGGGCTTTCGAAAGTTTACTCCGTGATTGCGAT----TACCCTGAAGCCAGAGCTACACTTCCACG-GACCGAACGC---CGAATATCCGGATCCTGCCTG-GCTTTCCGGCTTCGGAACTGAGAGAGGATCCGAGAGATAGCTGGTAA'
sequences1 = [seq5, seq6]
matrix = consensus_sequence(sequences1, DNA_KEY)
print(seq5)
print(seq6)
print(matrix)'''




sequences2 = []
with open('sequences.txt', 'r') as sequences_txt:
    for line in sequences_txt:
        sequences2.append(line[:-1].upper())

multiple_alignment1 = greedy_multiple_alignment(BLOSUM62, BLOSUM62_KEY, -5, sequences2[0:10])
for alignment in multiple_alignment1:
    print(alignment)
