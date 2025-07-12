from graph import Graph
from motif_finder import hamming_distance

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



def longest_path_to_alignment(sequence1, sequence2, longest_path):
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



def longest_path_to_alignment_affine_gap(sequence1, sequence2, longest_path):
    longest_path = id_to_tuple_affine_gap(longest_path)
    aligned_sequence1 = ''
    sequence1_index = 0
    aligned_sequence2 = ''
    sequence2_index = 0

    node = longest_path[0] # tuple of indices ex. ('M', '0', '1')
    next_node = longest_path[1]
    # accounts for local alignment and sets the correct starting indices
    if next_node[1] > (node[1] + 1) or next_node[2] > (node[2] + 1):
        sequence1_index = next_node[1]
        sequence2_index = next_node[2]


    for i in range(len(longest_path) - 1):
        node = longest_path[i] # tuple of indices
        next_node = longest_path[i+1]
        if node[0] == next_node[0] == 'M': # match or mismatch
            aligned_sequence1 += sequence1[sequence1_index]
            aligned_sequence2 += sequence2[sequence2_index]
            sequence1_index += 1
            sequence2_index += 1
        elif next_node[0] == 'D': # deletion
            aligned_sequence1 += sequence1[sequence1_index]
            aligned_sequence2 += '-'
            sequence1_index += 1
        elif next_node[0] == 'I': # insertion
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



def id_to_tuple_affine_gap(ids):
    ids_as_tuples = []
    for id in ids:
        indexes = id.split('-')
        if indexes[0] == 'D':
            ids_as_tuples.append(('D', int(indexes[1]), int(indexes[2])))
        elif indexes[0] == 'I':
            ids_as_tuples.append(('I', int(indexes[1]), int(indexes[2])))
        else:
            ids_as_tuples.append(('M',int(indexes[-2]), int(indexes[-1])))
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






def create_scored_global_alignment_graph(sequence1, sequence2, scoring_matrix, scoring_key, indel_score):
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

    return graph






def create_scored_local_alignment_graph(sequence1, sequence2, scoring_matrix, scoring_key, indel_score):
    graph = Graph()
    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            graph.add_node(f'{i}-{j}')

    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            graph.add_edge('0-0', f'{i}-{j}', weight=0) # add 0 weight edge from 0-0 to every node
            if i != len(sequence1)  and j != len(sequence2):
                # add a 0 weight edge from every node to the sink node
                graph.add_edge(f'{i}-{j}', f'{len(sequence1)}-{len(sequence2)}', weight=0)
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
    graph.remove_edge('0-0', '0-0', 0) # remove self cycle from 0-0 to itself
    return graph






def create_affine_gap_global_alignment_graph(sequence1, sequence2, scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost):
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






def create_affine_gap_local_alignment_graph(sequence1, sequence2, scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost):
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

    # remove self cycle from 0-0 to itself
    graph.remove_edge('0-0', '0-0', 0)
    # remove self cycle from sink node to itself
    graph.remove_edge(f'{len(sequence1)}-{len(sequence2)}', f'{len(sequence1)}-{len(sequence2)}', 0)
    return graph






def create_scored_fitting_alignment_graph(container_sequence, fitting_sequence, scoring_matrix, scoring_key, indel_score):
    graph = Graph()
    fitting_length = len(fitting_sequence)
    container_length = len(container_sequence)
    for i in range(container_length + 1):
        for j in range(fitting_length + 1):
            graph.add_node(f'{i}-{j}')

    for i in range(container_length + 1):
        for j in range(fitting_length + 1):
            if i < container_length:
                weight = indel_score
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j}', weight=weight) # deletion
            if j < fitting_length:
                weight = indel_score
                graph.add_edge(f'{i}-{j}', f'{i}-{j+1}', weight=weight) # insertion
            if (i < container_length) and (j < fitting_length):
                char1_index = scoring_key.index(container_sequence[i])
                char2_index = scoring_key.index(fitting_sequence[j])
                weight = scoring_matrix[char1_index][char2_index]
                graph.add_edge(f'{i}-{j}', f'{i+1}-{j+1}', weight=weight) # match or mismatch

    for i in range(container_length + 1):
        # add 0 weight edges from 0-0 to the first column nodes
        graph.add_edge('0-0', f'{i}-0', 0)
        # add 0 weight edges from the last column nodes to the sink node
        graph.add_edge(f'{i}-{fitting_length}', f'{container_length}-{fitting_length}', 0)

    # remove self loops in source and sink nodes
    graph.remove_edge('0-0', '0-0', 0)
    graph.remove_edge(f'{container_length}-{fitting_length}', f'{container_length}-{fitting_length}', 0)
    return graph






def create_scored_overlap_alignment_graph(sequence1, sequence2, scoring_matrix, scoring_key, indel_score):
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






def basic_pairwise_global_alignment(sequence1, sequence2):
    graph = create_basic_global_alignment_graph(sequence1, sequence2)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(sequence1, sequence2, longest_path), score



def scored_pairwise_global_alignment(sequence1, sequence2, scoring_matrix, scoring_key, indel_score):
    graph = create_scored_global_alignment_graph(sequence1, sequence2, scoring_matrix, scoring_key, indel_score)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(sequence1, sequence2, longest_path), score



def scored_pairwise_local_alignment(sequence1, sequence2, scoring_matrix, scoring_key, indel_score):
    graph = create_scored_local_alignment_graph(sequence1, sequence2, scoring_matrix, scoring_key, indel_score)
    longest_path,  max_score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(sequence1, sequence2, longest_path), max_score



def scored_pairwise_global_alignment_affine_gap_penalty(sequence1, sequence2, scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost):
    graph = create_affine_gap_global_alignment_graph(sequence1, sequence2, scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(sequence1, sequence2, longest_path), score



def scored_pairwise_local_alignment_affine_gap_penalty(sequence1, sequence2, scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost):
    graph = create_affine_gap_local_alignment_graph(sequence1, sequence2, scoring_matrix, scoring_key, initial_gap_cost, additional_gap_cost)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment_affine_gap(sequence1, sequence2, longest_path), score



def scored_pairwise_fitting_alignment(container_sequence, fitting_sequence, scoring_matrix, scoring_key, indel_score):
    graph = create_scored_fitting_alignment_graph(container_sequence, fitting_sequence, scoring_matrix, scoring_key, indel_score)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(container_sequence)}-{len(fitting_sequence)}')
    return longest_path_to_alignment(container_sequence, fitting_sequence, longest_path), score



def scored_pairwise_overlap_alignment(sequence1, sequence2, scoring_matrix, scoring_key, indel_score):
    graph = create_scored_overlap_alignment_graph(sequence1, sequence2, scoring_matrix, scoring_key, indel_score)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(sequence1, sequence2, longest_path), score


'''seq1 = 'CAAGGAACTCAATCGTCGAGTCCACGGGGGGCAGAACACGCTATATTTAATCTTGATGAGGAACGCAAATAACCATGGTTGCACGTGAGGATTTTCTTTAGTGAGTTGGGTTGCTTGGTAACTTATCCACTGCTATCTTAAGGGGGTTACTTCGGGATGAACGGCTTATGACAATCACAGTGAGGTCCGTCCCGGCCGATATGAGTTCTATGTTTTAACAGCGTCACCAGTGTCACGTACGGGGCCACCTCAGGCCCTGACCAGGGAATAGAGCGATTTGGGGACTTTCCCGGGTGATGTCTACCAGGAAGTTCGGTACCACTGACTTTGAATAATACTGTCAAAGGGGCTGCACCTTCCCGAGTTCGTCGTCATTACACAGCGCATATATTACACGTTAAGCCGTTTATCCGCATGTTATGCCAATTCGCGTCTTGCCAGGTGCCAACGAGCCTGATAAAGCAGTGGGTAGCGCCGGCACAGTATGTAGCAAGTTCCCCGCCGCGCGTTGAAAGCGTTACGTACAGGCGGCTAAGCGACGTTAAAATTGTCGCTTGCCTAACCCATCTCCCTGACACGGAACATAGCGAATAGTAGTCAACGGAGTTATGGTACAAAGCCTGAAAGCGACCTCAGACGAAGGGTCTGCCCGCAGGACGTGGGCTCTAATCCTCGGGGGCCTCGCCTACGTAGCACATCCCCAATAGCACTAAGAAGATGTGAACGAAACGCCGCTGTCGGATTCCAATTCTGAAATAGATAGTACCGGGTCCGAGGCGATGGAGGGTGGCGAAACCCCCATTTACGCATAGCGGTAACTTGGTCCCGGACTATTTATCAGTTGGTACCCTCGGCCCTGGTGGATGTGTTTTACTGGAATGATGTATAGACATCGCCTAC'
seq2 = 'CGCCCTACCTTTCGTGCCTATCAAAGTCGGTAGGCGCATTCCGAGCTCCGACTAGCGAGAAAACCAGCACAATGAAGTGCGCCTCCATATGTTAGTATACGGGTTGCCTCAGCTTTGGGCGGGCCTAAGGGCGGGATGAACGGCTTATTCCTGTGCAGTAGGTCGGTCCCGCCGATATGAGTTCTGGAGGTATTTAACAGCAGGGGCACGTGCACTGGGCCGCCTCAGGGGCTCACCATGACCAGGACATCGAGCGGGATCTCAACTTTGGGGACTTTCCCGGGTGATGTCTGCCACTGATGTTCGGTACCAGCGACCATACTGTGAAAAGGGGCTGCAAACATTATCTTCCCGAGTTGTCATTACACAGCATATATTACGCCGAAACTACGTCTCCGCATGTTTAACTCGCGTTTTGCCAGAGCCTAGCTAGGTTATAAAGGCGGCACAGTTTGTTCTAAAACCCGCCCTCGCCGAAGCGTTACGTACAGCGGCAAAGCGACGTTGAAATTGTCGATTGCCTATGGTTGACGCGGAACATAGCGTGACATAGTAGTGATACTCCAAGCGTGAAAGTGACCTAGCGGGTCTGAGCAGGACGTTGGCTCGAGCACGTGCAAGTTAGGGTTGTTCATGCCTGGGGAAAATAGTTTACGATTATTCGGCCATGCGTTAGCGTGCTGCTGTCCAGGCGCAGGGGCTTTCGAAAGTTTACTCCGTGATTGCGATTACCCTGAAGCCAGAGCTACACTTCCACGGACCGAACGCCGAATATCCGGATCCTGCCTGGCTTTCCGGCTTCGGAACTGAGAGAGGATCCGAGAGATAGCTGGTAA'
print(seq1.index('GCGGCTAAGCGACGTTAAAATTGTCGCTTGCCTA'))
print(seq1[:527])
print(seq1[527:])

print(seq2.index('GCGGCAAAGCGACGTTGAAATTGTCGATTGCCTA'))
print(seq2[:490])
print(seq2[490:])
'''


'''path1 = ['0-0', '512-604', '513-605', '514-606', '515-607', '516-608', '517-609', '518-610', '519-611', '520-612', '521-613', '522-614', '523-615', '524-616', '525-617', '526-618', '527-618', '528-619', '529-620', '530-621', '531-622', '532-623', '533-624', '534-625', '535-626', '536-627', '537-628', '538-629', '539-630', '540-631', '541-632', '542-633', '543-634', '544-635', '545-636', '546-637', '547-638', '548-639', '549-640', '550-641', '551-642', '552-643', '553-644', '554-645', '555-646', '556-647', '557-648', '558-649', '559-650', '560-651', '561-652']
seq1 = 'CAAGGAACTCAATCGTCGAGTCCACGGGGGGCAGAACACGCTATATTTAATCTTGATGAGGAACGCAAATAACCATGGTTGCACGTGAGGATTTTCTTTAGTGAGTTGGGTTGCTTGGTAACTTATCCACTGCTATCTTAAGGGGGTTACTTCGGGATGAACGGCTTATGACAATCACAGTGAGGTCCGTCCCGGCCGATATGAGTTCTATGTTTTAACAGCGTCACCAGTGTCACGTACGGGGCCACCTCAGGCCCTGACCAGGGAATAGAGCGATTTGGGGACTTTCCCGGGTGATGTCTACCAGGAAGTTCGGTACCACTGACTTTGAATAATACTGTCAAAGGGGCTGCACCTTCCCGAGTTCGTCGTCATTACACAGCGCATATATTACACGTTAAGCCGTTTATCCGCATGTTATGCCAATTCGCGTCTTGCCAGGTGCCAACGAGCCTGATAAAGCAGTGGGTAGCGCCGGCACAGTATGTAGCAAGTTCCCCGCCGCGCGTTGAAAGCGTTACGTACAGGCGGCTAAGCGACGTTAAAATTGTCGCTTGCCTAACCCATCTCCCTGACACGGAACATAGCGAATAGTAGTCAACGGAGTTATGGTACAAAGCCTGAAAGCGACCTCAGACGAAGGGTCTGCCCGCAGGACGTGGGCTCTAATCCTCGGGGGCCTCGCCTACGTAGCACATCCCCAATAGCACTAAGAAGATGTGAACGAAACGCCGCTGTCGGATTCCAATTCTGAAATAGATAGTACCGGGTCCGAGGCGATGGAGGGTGGCGAAACCCCCATTTACGCATAGCGGTAACTTGGTCCCGGACTATTTATCAGTTGGTACCCTCGGCCCTGGTGGATGTGTTTTACGATGCTATAGCGCGTATCGATTAGCTATGCTATCTATATTGCGCGCATATGCTAGGCTATGCTAGCTCTAGAGCAGCACACATATTCGATCGTATACGTACGTACGTACGTACGTCGTCGATGCTAGCTATCGATCGACTAGCTGGAATGATGTATAGACATCGCCTAC'
seq2 = 'CGCCCTACCTTTCGTGCCTATCAAAGATCTAGCTAGCTAGCGACGTAGCTATTATCTATCTCGAGCTACTATCGTACGAGCCGCGAGCTCTTCATCGTATATCGGCTATCGACTAGCAGCTAGCTAGCAGCCAGATCACTATTTCAGTTACGCATCGGTAGGCGCATTCCGAGCTCCGACTAGCGAGAAAACCAGCACAATGAAGTGCGCCTCCATATGTTAGTATACGGGTTGCCTCAGCTTTGGGCGGGCCTAAGGGCGGGATGAACGGCTTATTCCTGTGCAGTAGGTCGGTCCCGCCGATATGAGTTCTGGAGGTATTTAACAGCAGGGGCACGTGCACTGGGCCGCCTCAGGGGCTCACCATGACCAGGACATCGAGCGGGATCTCAACTTTGGGGACTTTCCCGGGTGATGTCTGCCACTGATGTTCGGTACCAGCGACCATACTGTGAAAAGGGGCTGCAAACATTATCTTCCCGAGTTGTCATTACACAGCATATATTACGCCGAAACTACGTCTCCGCATGTTTAACTCGCGTTTTGCCAGAGCCTAGCTAGGTTATAAAGGCGGCACAGTTTGTTCTAAAACCCGCCCTCGCCGAAGCGTTACGTACAGCGGCAAAGCGACGTTGAAATTGTCGATTGCCTATGGTTGACGCGGAACATAGCGTGACATAGTAGTGATACTCCAAGCGTGAAAGTGACCTAGCGGGTCTGAGCAGGACGTTGGCTCGAGCACGTGCAAGTTAGGGTTGTTCATGCCTGGGGAAAATAGTTTACGATTATTCGGCCATGCGTTAGCGTGCTGCTGTCCAGGCGCAGGGGCTTTCGAAAGTTTACTCCGTGATTGCGATTACCCTGAAGCCAGAGCTACACTTCCACGGACCGAACGCCGAATATCCGGATCCTGCCTGGCTTTCCGGCTTCGGAACTGAGAGAGGATCCGAGAGATAGCTGGTAA'
alignment3, max_score1 = scored_pairwise_local_alignment_affine_gap_penalty(seq1, seq2, DNA_MATRIX, DNA_KEY, -10, -5)
print(max_score1)
print(alignment3[0])
print(alignment3[1])
print(hamming_distance(alignment3[0], alignment3[1]))'''



'''seq3 = 'ENLHDMCNDHMIHGAVFKKQYKSGEPKEQIFWYPMIQFLISVCRSIVRNQKDADWPPHKHWERCYHNSTPTCCLFYGQHHHFDVAGMHKIAPMHLLVDGKTQMTILLRLQDHGAHWPYCIVTTRRAMGTEAEKSDYFNKTCKLRAVSWLYDQYLEPFNGYYSVCKNCGADRGECLMRYIAHPKMFCYYWNPNLRPPGHLMGNPRIWTMPAELFSYTFCTYYLDKLRWIQCNAMDWIHSYWETVDTNCMYITMFGLMQDNTVNTACCIPFIMGFIFSDMSNQGNSVEGMWPMMTNLQRARCNFHRSPKGEKYFTWGYDKAHGIGQAYPAKMMNQSMDQQTCRGAWYPLFGIWAVASFHSCLWESKKPHLAVDCAENMCAYHHMWLYSPGDGHVDARLDYTLTGKSKMPCWPKKQWADACDVRVCYAIFIKWLICNAEFGALMYQMCDKPVCFLMACKIQYTWMARCWADMWWNFRIEAQLDYDWVWRMLTNMMTPLCVHANYNQKYDEDSWCDSANSQFAWAWYVQWFNEWAQLTYMSEMQVWSFLSTGPIEIWMPTYEMGDASTNTDNLRWNYMFVLKVTHYVHFNDVHPYSRHNDPSLHIMWDMGSHCKFNVDMDYMIWMNHCRFCKLGTINYNNAPNEDTFNTRDWMVFSPYHSAVCEMYNWNGRPDGLVIDEAAPQSSMGINSDNRSSWNHRCNRYQCFDWYYMKDEYHERHNSSTFFYWGPIDIFQFKMPYCPRDKLGMGKWKQDNPDVVGCWGIGGVFFKCYKADLMVRSPAYHYGGCPHGVIQEPEVHRCRKCFLITHAFRGVINDWIFSMTQLDSSERHINIIDKKSHATCVYARAPLLPHVQFAVKTCWKWRDTVDVVAWQQYLITAAKGLAKIIMICWIYHLHEIHDQCPVYCGPIMPGNWQFNIEGQCCIMACKDTNQATAWHQLGSWGKHTEVRLAISNPHIHFSGSIVRLSVHCHHVSYHSIRALFETNMSCMSELQVTNLFCQFGRAANYKQSQEGKKCHNPVQQWHVLMWAMFTNEHPDEMKCAVMDFFPYPILNQNFADERHMQHFINVVTDYHQGEQYGYIVIEETEARIIVGCKSHDRLDDDLSSHVNHLYYLDLIRMDVMGDGMCMPYCPYCSQRSPPQISPMVHMGCFVWITFICKVQYIVCWIKTGPSLRYETYPWDSAKWLSWFAQINAWLVVIPCDFLMDTNNSCDHTIQLFALINMRFLSYLEWRMPAAIWAARFNCPETPFIPMYACFEMSCWLGNFPFVPHQDTTIYNTFNFSHHTFMKRKGHGVDCTKWTTREQAKCATASEENQTEEMNFDNYSFRIMRANFSLRFIFSYIERKLTSMNKTWWIGWSCRNPDQHSWLQMYFTHILCAEFEAFTAPPHYLHPMGEVMHIIHVHPRPPGMMDYLHHHMPQPFHQGDGDSIQSVIKSAAPYCFCQYGKIHYWACHIYDHDFAAMNTYITTEKAFHNYKYSQSPIMDIRGGPASTCSHRSNMAWNVRFGSNRKIYGYACMTCFSVIMMKMHRQKWKHGNCMQGAPRVCRLVAAGVFECRYHNDYQDVCTCNNRWTWEWHMMRHKDALLMSMPTTAIQPERSVISGEGENEVFAAKNQMVPQCPPTWETYRRNVNVQRSSCPRWTKVQTPHLSILVYNNRKQFMTNSSWMFAGTPTIDQMDDEPNYYKEWGNCTYAAMSQQKKVQIEFACCINHDNRFQLVAFHRAQAFNNIVYYVWFAKKRVQCPALDEHGAKMMMFYIIASYLFPNQIDFEWMKSPMEPCNMMDRKGRIRAQSKHFWMGRFYQPKAETARCSKWVNWMTIVAYWGSDHQVAKTPVRHWHSIRTVPMYVKFYCNIQPLSPEVAAPWNNIGMLFGNEHEMVARNALDKQDAPAFVAKFGHDYHAKEVFEKWAGTQPKTTNAGFISYWCDCYKRMAGCMSIWWGPVHRVPKFDWTCNLLCEWCRSMKQVDQTAIAYMTIMWCIRHFHTSKGWTHRMPAWEAIMCYVCELCDNLRERMRPFSLRLPCYEATCSDPCNEIWYFIAHPAAYEDQKTSVLKRGEDHVCSNVTHNKMKVQLSGVRHFEHYRTMQYKTRSCFKMCMRKMHFEMEQMDMIMWGVPDTYNMATWQVPHQMMSEHGHNDWRNSQHTPYHWIGRVFSGEWDQHNHAIDRRMQPNGRPDIPLSQRTEPDSLGSPAIMTLFKLISNKTNWTSQSGVNESYQFDHKNEDFNVEKSYPDIYNMQFFMQVLIRKIHKISFLEDYAVSPDWQMRKWSSHQQNNLRLPVRDRAADPSSSDFGFLDMPYDENMPNYNHYLFRQWADWMQPMEDEDSYRMQMYQNIVITEFPSFYIMCFHQFWSCVMWERTYGQLHLKNMFSRFPYGLQDHGFHPMHNSGNGPGLWARKRTYLGLQACPISVFYVRPIAPRWEGRNMTTTSHMWQAGQCWKTHPEPQPPLKGNDCSKEEWRTVFISGEEHQWWCKWRFRRPSMHSPCKYKLCWRLLDTRSKQSAPNLKGLAYNGKFLWDYEYHGLHNLEDDDFSLVIMQDRFMWQWWRRAGVMYDNTLPHFYRTWSSARLWHIVWACCMCQKWLPHNGRLGMWQMCFYNMPMMDNYEDAALLTGTPHLPIHMMMEIHFGWCGRLWHFGMQLVRMLPDLNWALYNMLCVTEWTMGHSHDHGWMNDRQVPRQCGWCTDQVMTPFLYFLNARKVHVLSCCD'
seq4 = 'LNQQIVLSGHKFRMKQELLSNNNSQMETMKPHEGACQDCLALERSNSTQCPQGWAIVWYDKLKMMMFAQDLVYEDMMSRRTNLMAWFQHSVCIWNIALNVGRNDEDSMDRSEHNVPSETQIQTMQLYMNHADLFWWPATFREWHGEKMYKTYDTWLWTIEINNCHMFAWVCSECPWWCDDWNMWDYADSYNIYVASFAFRNFGKGSAVDVMSWVPFMTIIKLRVTIVLYIQYTHHMTNVISHYTWMWMWRIMVCLMETQLIINIWFQQWVCSEACVCAAHQMMLTMNVHTVLQIVFMYMPCWPRVRPPMESFQDCRSGLCPRQTSFICRLADWRCRCVCPFHRTKQANCICCYDQSAITGPYRIIVATCGCMIKLAWVYNQESCANCTRRQWPGMAIETIIGKDFDANQVQIKRLLSNAPFYVIQEYNKAIHNYYVMCMFHLSAWWCEIIRTIQPPEDEHGMHNFWFGYPFSGNFEHLNEGAQKMVVQAPHYVYNRKWHHMYCWHYDYSLCRVAFRQKIKECQKVVWWLRPAITWAMETEQPWEYNCFKNFQNWLVWCTYPTHYHATCCYDHSCSTDVNSCNLITKFDCCINWMEIKHFTINRHPWIDNFIETHFPVHPERGMMLPCRESKDHRILERQSQQERDRHLFPMSLGCPNRKGLFYWCEKIPHEPCHGFQYNRNLAHTQWGYQVTGSNMQAHCTDWNWNAWQKHVIRIGMMNVHWMLATSVGRQRVDCHARLATTWYKTQWKWKEKEMWCMACHKTHNFYYDNQEFASPLWHCPWWPIESWWWNPKHGCNTCRGVGFRKEKQMKGRPHTGCHFQCMSMVSSERHNNIIDKKSHATCVYARAPLLSMHHFGYKHVQFAVKTCAPYTDRHQKWQCYDVDTTDVVAWQQYGITAAKGLAKIIISCPICWIYHLHYCRFVEHHDQCPVRDFLGQWGVSEGFQHIMPGNVRTSQFNIEGQMCIMACRDRNQATAWHQLGSWGKHIEVRLAVRLSVHCHNKIMFILVSYHSIRALFETNMSCMSELSVTNLFCQFGRAANPPSQEGKKCHNPVQDWHVLMWAMFTEESPDEMWMAMREDFVPYYMILERHMQHFINVVDKDYHQGQYGCIVIENIAARYCKSHDRLDDDEIYPWWYFSSHVNHLYEPGCIMLDLIRIPKLCPLDVMGDGMCMPYCPYCSQRNPPQISPMVHCFVWITFILKVQYFVCWIKTGPSLRKEKYAKWLSVFAQIIPCDFFRMDTNQALINMRQLSELTYRMPAAIWAQRFWCPETPFPMYASFEMKCWLGNFPFVPHQDTTIHHTFMVHPMRWPNYGHNNYWRSQAKCTPHRTASEENQTEAMNFDNYSKRIMRANFSLRFIFLTSMNKTWWDQHTWLQDMTFDVFQHILCAEFEAFTANTMSPPHYLHPMGEVMHIIHYPRPPGMMDVLHHHMPQPSWKKFWDRCNYCFCQYGKIPHWACHIGDHWFFAGNTYINTEKAFHNSKYSQSPIMYNIPFASTCHKMTGSHRSNMAWVVRFGKNRKIYGYTCMTCFSVITMKMHKRGCTLGFWWKHGNCMQSATAHIETWHCRLVAAFVFECRYHNYIGVMGYQDVCTCNNRWTWEQPTCHMMRHKQIHMWPDEASDFQKERMTNKVPRCPPWWPPTGCAVITYRRNVNHMQCQRSSCPKVQTPQDSILVYNNRKQLANAMFAGTPTPEEHYDRFDQMDDEPNYRHKERGNCTYPSKMGSANAMRQQGDCIRHDNEFQLVAFHRAQAFNNNVYNPWDNFQCPALDEHGAEMMMFYIIFSFQEMQYVEHAMCCNMRAQSKHFWMGRFYQPKAETRRCSKWVNWMTIVSNYWGSDHQVAKTPARTVYCDAQPLSMEECAPYVTLFYNEHEMVNALQKQDFGHDYHAKEVFDYIFLMKIGVTQLLPKDVTIISRPNNKAASDSIHKFTAKGLYVMMHTEAEGKFYEDKNQGWEIGFHIWGYQQFADSVVSVPTSISNWTVIAPQQKYYKVMRIIVSHKRDTSLILPNIMWIVWVNIFEDSMNCTLNHHNNPKAYAVPMGRIDYTIPIDIKPNVFFIMNCECEFDLEHHWCKALASCSDFLYADSLEGTGLGHRKSIVCCMYDYGECDCAEGPHTLGSYKRTAPHYTQVFERFMYPINMYMLTEAMSTGTVGHAIPDHDGECVQGSFDYNEIVMEHHIHRWRILDTKWKIRSALALKIGMYRWCCSKMCGFKPLHLGPPGSPPQYANPMEQYWLAWYYPCPYKQLLESYWMDANSQETPQVYNMGAEPTRFSIWFQNPYMENTWVHSHSERPCGYQAPYFKANCLRTCGSLWGFSTYPIYWGQSFYIQREMNQTWSNLWMTVMLLRYKGDTAECMVFCVWNHDEIFQTNKLQGNMALFMNGIWDTIYTELGHWWIHWPKGPPMCLPEPNMAWIMGSVEPSERQAGNKLTIQAYTFKSASYACHLVLYCRLCDHPASYRAFSVILMMSEDYYCQKDHKKFNEMWNIEMQNESRINVGPTPANNEFCADEMETELFCVGRDYAVNGDQQIGCVKYYHQNLTIAMVREANRECFQSIDWSRYWCRVGGQKLSPVKEWWIKVYEGQPCDMFCRTPEKEVHLFAQNPTTKWHVCSYLKDLARIVINSAMNASIAKDSSCIARQAKKPMPTYHPGEIPMGVPGDQGNYAFHKKERTHHEPWMENACFASNWCM'
alignment4, max_score2 = scored_pairwise_alignment_affine_gap_penalty(seq3,seq4, BLOSUM62, BLOSUM62_key, -10, -2)
print(max_score2)
print(alignment4[0])
print(alignment4[1])
print(hamming_distance(alignment4[0], alignment4[1]))
path2 = ['0-0', '1-119', '2-120', '3-121', '4-122', '5-123', '6-124', '7-124', '8-125', '9-126', '10-127', '11-128', '12-129', '13-130', '14-131', '15-132', '16-133', '17-134', '17-135', '17-136', '18-137', '19-138', '20-139', '21-140', '22-141', '23-142', '23-143', '23-144', '24-145', '25-146', '26-147', '27-148', '28-149', '29-150', '30-151', '31-152', '31-153', '31-154', '32-155', '33-156', '34-157', '35-158', '36-159', '37-160', '38-160', '39-160', '40-161', '41-162', '42-163', '43-164', '44-165', '45-166', '46-167', '47-168', '48-169', '49-170', '50-171', '51-172', '52-173', '53-174', '54-175', '55-176', '56-177', '57-178', '58-179', '59-179', '60-180', '61-181', '62-182', '62-183', '63-184', '64-185', '65-186', '66-187', '67-188', '68-189', '69-190', '70-191', '71-192', '72-193', '73-194', '74-195', '75-196', '76-197', '77-198', '78-198', '79-199', '80-200', '81-201', '82-202', '83-203', '84-204', '85-205', '86-206', '87-207', '88-208', '89-209', '90-210', '91-211', '92-212', '92-213', '93-214', '94-215', '95-216', '96-217', '97-218', '98-219', '99-220', '100-221', '101-222', '102-223', '103-224', '104-225', '105-226', '106-227', '107-228', '108-229', '109-230', '110-231', '111-231', '112-232', '113-233', '114-234', '115-235', '116-236', '117-237', '118-238', '119-239', '120-240', '121-241', '122-242', '122-243', '123-244', '124-245', '125-246', '126-247', '127-248', '128-248', '129-249', '130-250', '131-251', '132-251', '133-252', '134-253', '135-253', '136-254', '137-255', '138-256', '139-257', '140-258', '141-258', '142-259', '143-260', '144-261', '145-262', '146-263', '147-264', '148-265', '149-265', '150-266', '151-267', '152-268', '153-269', '154-270', '155-270', '156-271', '157-271', '158-272', '159-273', '160-273', '161-273', '162-274', '163-274', '164-275', '165-275', '166-276', '167-277', '168-278', '169-279', '170-280', '171-281', '172-282', '173-283', '174-283', '175-284', '175-285', '176-286', '177-287', '177-288', '178-289', '179-290', '180-291', '181-292', '182-293', '183-294', '184-295', '185-296', '186-297', '187-298', '187-299', '187-300', '188-301', '189-302', '190-302', '191-303', '192-304', '193-305', '194-306', '195-307', '196-308', '197-309', '198-310', '199-311', '200-312', '201-313', '202-314', '203-315', '204-316', '205-317', '206-317', '207-318', '208-319', '208-320', '209-321', '210-322', '211-323', '212-324', '213-324', '214-325', '215-326', '216-326', '217-327', '218-328', '219-329', '220-329', '221-329', '222-330', '223-331', '224-332', '225-332', '226-332', '227-333', '228-333', '229-334', '230-335', '231-336', '232-337', '233-337', '234-337', '235-337', '236-338', '237-339', '238-340', '239-341', '240-342', '241-343', '242-344', '243-345', '244-346', '245-347', '246-348', '247-349', '248-350', '248-351', '248-352', '249-353', '250-353', '251-354', '252-355', '253-356', '254-357', '255-358', '256-359', '257-360', '258-361', '258-362', '259-363', '260-364', '261-365', '262-366', '263-367', '264-368', '265-369', '265-370', '266-371', '267-372', '268-372', '269-373', '270-374', '271-375', '272-376', '273-377', '274-378', '275-379', '276-380', '277-381', '278-381', '279-382', '280-383', '281-384', '282-385', '283-386', '284-387', '285-388', '286-389', '287-390', '288-391', '289-392', '290-393', '290-394', '291-395', '292-396', '293-397', '294-398', '295-399', '296-400', '297-401', '298-402', '299-403', '300-403', '301-404', '302-405', '303-406', '304-407', '305-408', '306-409', '307-410', '308-411', '309-412', '310-413', '311-414', '312-415', '313-415', '314-416', '315-417', '316-418', '317-419', '318-420', '319-421', '320-422', '321-423', '322-424', '323-424', '324-425', '325-426', '326-427', '327-427', '328-428', '329-429', '330-430', '331-431', '332-431', '333-432', '334-433', '335-434', '336-434', '337-435', '338-436', '339-437', '340-438', '341-439', '342-439', '343-439', '344-439', '345-440', '346-441', '347-442', '348-442', '349-443', '350-444', '351-445', '352-445', '353-445', '354-445', '355-445', '356-446', '357-446', '358-446', '359-447', '360-447', '361-447', '362-448', '363-449', '364-450', '365-451', '366-451', '367-451', '368-451', '369-452', '370-453', '371-454', '372-455', '373-456', '374-457', '375-458', '376-459', '377-460', '378-461', '379-462', '380-463', '381-464', '382-465', '383-466', '384-467', '384-468', '385-469', '386-470', '386-471', '387-472', '388-473', '389-474', '389-475', '390-476', '391-477', '392-478', '393-479', '394-480', '395-481', '396-482', '397-483', '398-483', '399-484', '400-485', '401-486', '402-487', '403-488', '404-489', '405-490', '406-491', '407-492', '408-493', '409-494', '410-495', '411-496', '412-497', '413-497', '414-498', '415-499', '416-500', '417-501', '417-502', '418-503', '419-503', '420-503', '421-504', '422-505', '423-505', '424-506', '425-507', '426-507', '427-508', '428-508', '429-509', '430-509', '431-510', '432-510', '433-511', '434-512', '435-513', '436-514', '437-515', '438-516', '439-517', '440-517', '441-518', '442-519', '443-520', '444-521', '445-522', '446-523', '447-524', '448-525', '449-526', '450-527', '451-528', '452-529', '453-530', '454-531', '455-531', '456-532', '457-533', '458-533', '459-533', '460-534', '461-535', '462-536', '463-537', '464-538', '465-539', '466-539', '467-540', '468-541', '469-542', '470-542', '471-543', '472-544', '473-545', '474-546', '475-547', '476-548', '477-549', '478-550', '479-551', '480-552', '481-552', '482-553', '483-554', '483-555', '484-556', '485-557', '486-558', '487-559', '488-560', '489-561', '490-562', '491-563', '492-564', '493-565', '494-566', '495-567', '496-568', '497-568', '498-568', '499-569', '500-569', '501-570', '502-571', '503-572', '504-573', '505-574', '506-575', '507-576', '508-577', '508-578', '509-579', '510-580', '511-581', '512-582', '513-582', '514-583', '515-584', '516-585', '517-586', '518-587', '519-588', '520-588', '521-589', '522-589', '523-590', '524-591', '525-592', '526-593', '527-594', '528-594', '529-595', '530-595', '531-596', '532-597', '533-597', '534-598', '535-599', '536-600', '537-601', '538-602', '539-603', '540-604', '541-605', '542-606', '542-607', '542-608', '543-609', '544-610', '545-611', '546-612', '547-613', '548-613', '549-614', '550-615', '551-616', '552-617', '553-617', '554-618', '555-619', '556-619', '557-619', '558-620', '559-621', '560-622', '561-623', '562-624', '563-625', '564-626', '565-627', '566-628', '567-629', '568-630', '569-630', '570-631', '571-631', '572-632', '573-633', '574-634', '575-634', '576-635', '577-636', '578-637', '579-637', '580-638', '581-639', '582-639', '583-640', '584-641', '585-641', '586-642', '587-643', '588-643', '589-644', '590-644', '591-644', '592-645', '593-646', '594-647', '595-647', '596-647', '597-647', '598-647', '599-648', '600-649', '601-650', '602-651', '603-651', '604-652', '605-653', '606-654', '607-654', '608-654', '609-655', '610-656', '611-657', '612-658', '613-659', '614-660', '615-661', '616-661', '617-662', '618-662', '619-663', '620-664', '621-664', '622-664', '623-664', '624-665', '625-666', '626-666', '627-666', '628-667', '629-668', '630-669', '631-670', '632-671', '633-672', '634-673', '635-674', '636-675', '637-676', '638-677', '639-678', '640-679', '641-680', '642-681', '643-682', '644-683', '645-684', '646-685', '647-686', '648-687', '648-688', '648-689', '649-690', '650-691', '651-692', '652-693', '653-694', '654-695', '655-696', '656-697', '657-698', '658-699', '659-700', '660-701', '661-702', '662-703', '663-704', '664-705', '665-706', '666-707', '667-708', '668-709', '669-710', '670-711', '671-712', '672-713', '672-714', '673-715', '674-716', '675-717', '676-718', '677-719', '678-720', '679-721', '680-722', '681-723', '682-724', '683-725', '684-726', '685-727', '685-728', '686-729', '687-730', '688-731', '689-732', '690-733', '691-734', '692-735', '693-736', '694-737', '695-738', '695-739', '695-740', '696-741', '697-742', '698-743', '699-744', '700-745', '701-746', '702-746', '703-747', '704-748', '705-749', '706-750', '707-751', '707-752', '708-753', '709-754', '710-755', '710-756', '711-757', '712-758', '713-759', '714-760', '715-761', '716-762', '717-763', '718-764', '719-765', '720-766', '721-767', '722-768', '723-768', '724-769', '725-770', '726-771', '727-772', '728-772', '729-773', '730-774', '731-775', '732-776', '733-777', '734-778', '735-779', '736-780', '737-781', '738-782', '739-783', '740-784', '741-785', '742-786', '743-786', '744-787', '745-788', '746-789', '747-790', '748-791', '749-792', '750-793', '751-794', '752-795', '753-796', '754-797', '755-798', '756-799', '757-800', '758-801', '759-802', '760-802', '761-803', '762-803', '763-804', '764-805', '765-806', '766-806', '767-807', '768-808', '769-808', '770-809', '771-810', '772-811', '773-812', '774-813', '775-814', '775-815', '776-816', '777-817', '778-818', '779-819', '780-820', '781-820', '782-821', '783-822', '784-822', '785-823', '786-824', '787-825', '788-826', '789-826', '790-827', '791-828', '792-829', '793-829', '794-829', '795-830', '796-830', '797-831', '798-832', '799-832', '800-833', '801-834', '802-835', '803-836', '804-837', '805-838', '806-839', '807-840', '808-841', '809-842', '810-843', '811-843', '812-843', '813-843', '814-844', '815-845', '816-846', '817-847', '818-848', '819-849', '820-850', '821-850', '822-851', '823-852', '824-853', '825-854', '826-855', '827-856', '828-857', '829-858', '830-859', '831-860', '832-861', '833-862', '833-863', '834-864', '835-865', '836-866', '837-867', '838-868', '839-869', '839-870', '840-871', '841-872', '842-872', '843-873', '844-873', '845-873', '846-873', '847-874', '848-875', '849-876', '850-877', '851-878', '852-878', '853-878', '854-879', '855-879', '856-880', '857-881', '858-882', '859-882', '860-883', '861-884', '862-885', '863-886', '864-887', '865-888', '866-889', '867-890', '868-891', '869-892', '870-893', '871-894', '872-895', '873-896', '874-897', '875-898', '876-899', '877-900', '878-901', '879-902', '880-903', '881-904', '882-905', '883-906', '884-907', '884-908', '884-909', '884-910', '885-911', '886-912', '887-913', '888-914', '889-915', '890-916', '891-917', '892-918', '892-919', '892-920', '893-921', '893-922', '894-923', '894-924', '894-925', '895-926', '896-927', '897-928', '898-929', '899-930', '900-931', '900-932', '900-933', '901-934', '901-935', '901-936', '901-937', '901-938', '901-939', '901-940', '902-941', '902-942', '903-943', '903-944', '903-945', '904-946', '905-947', '906-948', '907-949', '908-950', '909-951', '909-952', '910-953', '910-954', '910-955', '911-956', '912-957', '913-958', '914-959', '915-960', '916-961', '917-962', '918-963', '919-964', '920-965', '921-966', '922-967', '923-968', '924-969', '925-970', '926-971', '927-972', '928-973', '929-974', '930-975', '931-976', '932-977', '933-978', '934-979', '935-980', '936-981', '937-982', '938-983', '939-984', '940-985', '941-986', '942-987', '943-988', '944-989', '945-990', '946-991', '947-992', '948-993', '949-994', '950-995', '951-996', '952-996', '953-997', '954-998', '955-999', '956-1000', '957-1001', '958-1002', '959-1003', '960-1004', '961-1004', '962-1005', '963-1005', '964-1006', '965-1006', '966-1006', '967-1006', '968-1007', '969-1008', '970-1009', '971-1010', '972-1011', '973-1012', '974-1013', '975-1014', '976-1015', '977-1016', '978-1017', '979-1018', '980-1019', '981-1020', '982-1021', '983-1022', '984-1023', '985-1024', '986-1025', '987-1026', '988-1027', '989-1028', '990-1029', '991-1030', '992-1031', '993-1032', '994-1033', '995-1034', '996-1035', '997-1036', '998-1037', '999-1038', '1000-1039', '1001-1040', '1002-1041', '1003-1041', '1004-1042', '1005-1043', '1006-1044', '1007-1045', '1008-1046', '1009-1047', '1010-1048', '1011-1049', '1012-1050', '1013-1051', '1014-1052', '1015-1053', '1016-1054', '1017-1055', '1018-1056', '1019-1057', '1020-1058', '1021-1059', '1022-1060', '1023-1061', '1024-1062', '1025-1063', '1026-1064', '1027-1065', '1028-1066', '1029-1067', '1030-1068', '1031-1069', '1032-1070', '1033-1071', '1034-1072', '1035-1073', '1036-1074', '1037-1075', '1038-1076', '1039-1077', '1040-1078', '1040-1079', '1041-1080', '1042-1081', '1043-1082', '1044-1083', '1045-1084', '1046-1084', '1047-1085', '1048-1086', '1049-1086', '1050-1086', '1051-1087', '1052-1088', '1053-1088', '1054-1088', '1055-1089', '1056-1090', '1057-1091', '1058-1092', '1059-1093', '1060-1094', '1061-1095', '1062-1096', '1063-1097', '1064-1098', '1065-1099', '1065-1100', '1066-1101', '1067-1102', '1068-1103', '1069-1104', '1070-1105', '1071-1106', '1072-1106', '1073-1107', '1074-1108', '1075-1109', '1076-1110', '1077-1111', '1078-1112', '1079-1113', '1080-1114', '1081-1115', '1082-1116', '1083-1117', '1084-1118', '1085-1119', '1086-1119', '1087-1120', '1088-1120', '1089-1120', '1090-1121', '1091-1122', '1092-1123', '1093-1124', '1094-1125', '1095-1126', '1096-1127', '1097-1128', '1098-1129', '1099-1130', '1099-1131', '1099-1132', '1099-1133', '1099-1134', '1099-1135', '1099-1136', '1099-1137', '1100-1138', '1101-1139', '1102-1140', '1103-1141', '1104-1142', '1105-1143', '1106-1144', '1107-1145', '1108-1146', '1108-1147', '1108-1148', '1108-1149', '1109-1150', '1109-1151', '1109-1152', '1110-1153', '1111-1154', '1112-1155', '1113-1156', '1114-1157', '1114-1158', '1114-1159', '1114-1160', '1114-1161', '1114-1162', '1114-1163', '1115-1164', '1116-1165', '1117-1166', '1118-1167', '1119-1168', '1120-1169', '1121-1170', '1122-1171', '1123-1172', '1124-1173', '1125-1174', '1126-1175', '1127-1176', '1128-1177', '1129-1178', '1130-1179', '1131-1180', '1132-1181', '1133-1182', '1134-1183', '1135-1184', '1136-1185', '1137-1186', '1138-1187', '1139-1188', '1140-1189', '1141-1190', '1142-1191', '1143-1192', '1144-1192', '1145-1192', '1146-1193', '1147-1194', '1148-1195', '1149-1196', '1150-1197', '1151-1198', '1152-1199', '1153-1200', '1154-1201', '1155-1202', '1156-1203', '1157-1204', '1158-1205', '1159-1206', '1160-1207', '1161-1208', '1162-1209', '1163-1210', '1164-1211', '1165-1212', '1166-1213', '1167-1214', '1168-1215', '1169-1216', '1170-1217', '1171-1218', '1172-1219', '1173-1220', '1174-1221', '1175-1221', '1176-1221', '1177-1221', '1178-1221', '1179-1222', '1180-1223', '1181-1224', '1182-1225', '1183-1226', '1184-1227', '1185-1228', '1186-1229', '1187-1230', '1188-1231', '1189-1231', '1190-1231', '1191-1231', '1192-1231', '1193-1231', '1194-1231', '1195-1232', '1196-1233', '1197-1234', '1198-1235', '1199-1236', '1200-1237', '1200-1238', '1201-1239', '1202-1240', '1203-1241', '1204-1241', '1205-1241', '1206-1241', '1207-1241', '1208-1241', '1209-1242', '1210-1242', '1211-1242', '1212-1243', '1213-1243', '1214-1243', '1215-1244', '1216-1245', '1217-1246', '1218-1247', '1219-1248', '1220-1249', '1221-1250', '1222-1251', '1223-1252', '1224-1253', '1225-1254', '1226-1255', '1227-1256', '1228-1257', '1229-1258', '1230-1259', '1231-1260', '1232-1261', '1233-1262', '1234-1263', '1235-1264', '1236-1265', '1237-1266', '1238-1267', '1239-1268', '1240-1269', '1241-1270', '1242-1271', '1243-1272', '1244-1273', '1245-1274', '1246-1274', '1247-1275', '1248-1276', '1249-1277', '1250-1278', '1251-1279', '1252-1280', '1253-1281', '1254-1282', '1255-1283', '1256-1284', '1257-1285', '1258-1286', '1259-1287', '1260-1288', '1261-1289', '1262-1290', '1263-1291', '1264-1292', '1265-1293', '1266-1294', '1267-1295', '1268-1296', '1269-1297', '1270-1298', '1271-1299', '1272-1300', '1273-1301', '1274-1302', '1275-1303', '1276-1303', '1277-1304', '1278-1305', '1279-1305', '1280-1306', '1281-1307', '1282-1307', '1283-1308', '1284-1309', '1285-1310', '1286-1311', '1287-1312', '1288-1313', '1289-1314', '1290-1315', '1291-1316', '1292-1316', '1293-1317', '1294-1318', '1295-1319', '1296-1319', '1297-1319', '1298-1320', '1299-1321', '1300-1322', '1301-1323', '1302-1324', '1303-1325', '1303-1326', '1304-1327', '1304-1328', '1304-1329', '1305-1330', '1306-1331', '1307-1332', '1308-1333', '1309-1334', '1310-1335', '1311-1336', '1312-1337', '1313-1338', '1314-1339', '1315-1340', '1316-1341', '1317-1342', '1318-1343', '1319-1344', '1320-1345', '1321-1346', '1322-1347', '1323-1348', '1324-1349', '1325-1350', '1326-1351', '1327-1352', '1328-1353', '1329-1354', '1330-1355', '1331-1356', '1332-1357', '1333-1358', '1334-1359', '1335-1360', '1336-1360', '1337-1360', '1338-1360', '1339-1360', '1340-1360', '1341-1360', '1342-1361', '1343-1362', '1344-1363', '1345-1364', '1346-1365', '1347-1366', '1348-1367', '1349-1367', '1350-1368', '1351-1368', '1352-1368', '1353-1369', '1354-1369', '1355-1369', '1356-1369', '1357-1369', '1358-1369', '1359-1370', '1360-1371', '1361-1372', '1362-1373', '1363-1374', '1364-1375', '1365-1376', '1365-1377', '1366-1378', '1366-1379', '1367-1380', '1367-1381', '1367-1382', '1368-1383', '1369-1384', '1370-1385', '1371-1386', '1372-1387', '1373-1388', '1374-1389', '1375-1390', '1376-1391', '1377-1392', '1378-1393', '1379-1394', '1380-1395', '1381-1396', '1381-1397', '1381-1398', '1381-1399', '1381-1400', '1382-1401', '1383-1402', '1384-1403', '1385-1404', '1386-1405', '1387-1406', '1388-1407', '1389-1408', '1390-1409', '1391-1410', '1392-1411', '1393-1412', '1394-1413', '1395-1414', '1396-1415', '1397-1416', '1398-1416', '1399-1417', '1400-1418', '1401-1419', '1402-1420', '1403-1421', '1404-1422', '1405-1423', '1406-1424', '1407-1425', '1408-1426', '1409-1427', '1410-1428', '1411-1429', '1412-1430', '1413-1431', '1414-1432', '1415-1433', '1416-1434', '1417-1435', '1418-1436', '1419-1437', '1420-1437', '1421-1437', '1422-1437', '1423-1437', '1424-1438', '1425-1439', '1426-1439', '1427-1440', '1428-1440', '1429-1441', '1430-1442', '1431-1443', '1432-1443', '1433-1443', '1434-1444', '1435-1445', '1436-1446', '1437-1447', '1438-1448', '1439-1449', '1440-1450', '1441-1451', '1442-1452', '1443-1453', '1444-1454', '1445-1455', '1446-1456', '1447-1457', '1448-1458', '1449-1459', '1450-1460', '1451-1461', '1452-1462', '1453-1463', '1454-1464', '1455-1465', '1456-1466', '1457-1467', '1458-1468', '1459-1469', '1460-1470', '1461-1471', '1462-1472', '1463-1473', '1464-1474', '1465-1475', '1466-1476', '1467-1477', '1468-1478', '1469-1479', '1470-1480', '1471-1481', '1472-1482', '1473-1483', '1474-1484', '1475-1485', '1476-1486', '1477-1487', '1478-1488', '1479-1489', '1479-1490', '1480-1491', '1481-1492', '1482-1493', '1482-1494', '1483-1495', '1484-1496', '1485-1497', '1485-1498', '1486-1499', '1487-1500', '1487-1501', '1488-1502', '1489-1503', '1490-1504', '1491-1505', '1492-1506', '1493-1507', '1494-1508', '1495-1509', '1496-1510', '1497-1511', '1498-1512', '1499-1513', '1500-1514', '1501-1515', '1502-1516', '1503-1517', '1504-1518', '1505-1519', '1506-1520', '1507-1521', '1508-1522', '1509-1523', '1510-1524', '1511-1525', '1512-1526', '1513-1527', '1514-1528', '1515-1529', '1516-1530', '1517-1531', '1518-1532', '1519-1533', '1520-1534', '1521-1535', '1522-1536', '1523-1537', '1524-1538', '1524-1539', '1525-1540', '1526-1541', '1526-1542', '1527-1543', '1527-1544', '1527-1545', '1527-1546', '1527-1547', '1528-1548', '1529-1549', '1530-1550', '1531-1551', '1532-1552', '1533-1553', '1534-1554', '1535-1555', '1536-1556', '1537-1557', '1537-1558', '1538-1559', '1539-1560', '1540-1561', '1540-1562', '1540-1563', '1540-1564', '1540-1565', '1541-1566', '1542-1567', '1543-1568', '1544-1569', '1545-1570', '1546-1571', '1547-1572', '1548-1573', '1549-1574', '1550-1575', '1551-1576', '1552-1577', '1553-1578', '1554-1579', '1555-1580', '1555-1581', '1555-1582', '1555-1583', '1555-1584', '1555-1585', '1556-1586', '1557-1587', '1558-1588', '1559-1589', '1560-1590', '1561-1591', '1562-1592', '1563-1593', '1564-1594', '1565-1595', '1566-1596', '1567-1597', '1568-1598', '1569-1599', '1570-1600', '1570-1601', '1570-1602', '1571-1603', '1571-1604', '1572-1605', '1573-1606', '1574-1607', '1575-1608', '1576-1609', '1577-1610', '1578-1611', '1579-1611', '1580-1612', '1581-1613', '1582-1614', '1583-1615', '1584-1615', '1585-1616', '1586-1617', '1587-1618', '1588-1619', '1589-1620', '1590-1621', '1591-1622', '1592-1623', '1593-1624', '1594-1625', '1595-1626', '1596-1627', '1597-1628', '1598-1629', '1599-1630', '1600-1631', '1601-1631', '1602-1632', '1603-1633', '1604-1634', '1605-1634', '1606-1635', '1607-1636', '1608-1637', '1609-1638', '1610-1639', '1611-1639', '1612-1640', '1613-1641', '1614-1642', '1615-1643', '1616-1644', '1617-1645', '1618-1646', '1619-1646', '1620-1646', '1621-1647', '1622-1648', '1623-1649', '1624-1650', '1625-1651', '1626-1652', '1626-1653', '1627-1654', '1628-1655', '1628-1656', '1628-1657', '1629-1658', '1630-1659', '1631-1660', '1632-1661', '1633-1662', '1634-1663', '1635-1663', '1636-1663', '1637-1663', '1638-1664', '1639-1665', '1640-1666', '1641-1667', '1642-1668', '1643-1669', '1644-1670', '1645-1671', '1646-1672', '1647-1673', '1648-1674', '1649-1675', '1650-1676', '1651-1677', '1652-1678', '1653-1679', '1654-1680', '1655-1680', '1656-1681', '1657-1682', '1658-1683', '1659-1683', '1660-1684', '1661-1684', '1662-1685', '1663-1686', '1664-1687', '1665-1688', '1666-1689', '1667-1690', '1668-1691', '1668-1692', '1668-1693', '1668-1694', '1668-1695', '1668-1696', '1668-1697', '1668-1698', '1669-1699', '1670-1700', '1671-1701', '1672-1702', '1673-1703', '1674-1704', '1675-1705', '1676-1706', '1677-1707', '1678-1708', '1678-1709', '1679-1710', '1680-1711', '1681-1712', '1682-1713', '1683-1714', '1684-1715', '1685-1716', '1686-1717', '1687-1718', '1688-1719', '1689-1720', '1689-1721', '1690-1722', '1691-1723', '1692-1724', '1693-1725', '1694-1726', '1695-1727', '1696-1728', '1697-1729', '1698-1730', '1699-1731', '1700-1731', '1701-1732', '1702-1733', '1703-1734', '1704-1735', '1705-1736', '1706-1737', '1707-1738', '1708-1739', '1709-1740', '1710-1741', '1711-1742', '1712-1743', '1713-1744', '1714-1745', '1715-1746', '1716-1747', '1717-1748', '1718-1749', '1719-1750', '1720-1751', '1721-1752', '1722-1753', '1723-1754', '1724-1755', '1725-1756', '1726-1757', '1727-1758', '1728-1759', '1729-1760', '1730-1760', '1731-1760', '1732-1761', '1733-1762', '1734-1762', '1735-1763', '1736-1764', '1737-1765', '1738-1766', '1739-1767', '1740-1768', '1741-1769', '1742-1770', '1743-1771', '1744-1772', '1745-1773', '1746-1774', '1747-1775', '1748-1776', '1749-1777', '1750-1778', '1751-1779', '1752-1780', '1753-1781', '1754-1781', '1755-1781', '1756-1782', '1757-1783', '1758-1784', '1759-1784', '1760-1785', '1761-1786', '1762-1787', '1763-1788', '1764-1789', '1765-1789', '1766-1789', '1767-1790', '1768-1791', '1769-1792', '1770-1793', '1771-1794', '1772-1794', '1773-1795', '1774-1796', '1775-1797', '1776-1797', '1777-1798', '1778-1798', '1779-1798', '1780-1798', '1781-1798', '1782-1798', '1783-1798', '1784-1799', '1785-1800', '1786-1801', '1787-1802', '1788-1803', '1789-1804', '1790-1805', '1791-1806', '1792-1807', '1793-1808', '1794-1809', '1795-1810', '1796-1811', '1797-1812', '1798-1813', '1799-1814', '1800-1815', '1801-1816', '1802-1817', '1803-1818', '1804-1819', '1805-1820', '1806-1821', '1807-1822', '1808-1823', '1809-1824', '1810-1825', '1811-1826', '1812-1827', '1813-1828', '1814-1829', '1815-1830', '1816-1831', '1816-1832', '1817-1833', '1818-1834', '1819-1835', '1820-1836', '1821-1837', '1822-1838', '1823-1839', '1824-1840', '1825-1841', '1826-1842', '1827-1843', '1828-1844', '1829-1844', '1830-1844', '1831-1844', '1832-1844', '1833-1844', '1834-1845', '1835-1845', '1836-1846', '1837-1847', '1838-1847', '1839-1847', '1840-1847', '1841-1847', '1842-1848', '1843-1848', '1844-1848', '1845-1849', '1846-1850', '1847-1851', '1848-1852', '1849-1853', '1850-1854', '1851-1855', '1852-1856', '1853-1857', '1854-1858', '1855-1859', '1856-1860', '1857-1861', '1858-1862', '1859-1863', '1860-1863', '1861-1863', '1862-1864', '1863-1865', '1864-1866', '1865-1867', '1866-1868', '1867-1868', '1868-1869', '1869-1870', '1870-1871', '1871-1872', '1872-1873', '1873-1874', '1874-1874', '1875-1874', '1876-1875', '1877-1876', '1878-1877', '1879-1878', '1880-1879', '1881-1880', '1882-1881', '1883-1881', '1884-1881', '1885-1881', '1886-1881', '1887-1881', '1888-1881', '1889-1881', '1890-1882', '1891-1883', '1892-1884', '1893-1885', '1894-1886', '1895-1887', '1896-1888', '1897-1889', '1898-1890', '1899-1891', '1900-1892', '1901-1893', '1902-1893', '1903-1894', '1904-1895', '1905-1895', '1906-1896', '1907-1897', '1908-1898', '1909-1899', '1910-1900', '1911-1901', '1912-1902', '1913-1903', '1914-1904', '1915-1905', '1916-1906', '1917-1907', '1918-1908', '1919-1909', '1920-1910', '1921-1911', '1922-1912', '1923-1913', '1924-1914', '1925-1915', '1926-1916', '1927-1917', '1928-1918', '1929-1918', '1930-1919', '1931-1920', '1932-1921', '1933-1921', '1934-1922', '1935-1923', '1936-1924', '1937-1925', '1938-1926', '1939-1927', '1939-1928', '1940-1929', '1941-1930', '1942-1931', '1943-1931', '1944-1932', '1945-1932', '1946-1933', '1947-1934', '1948-1935', '1949-1936', '1950-1937', '1950-1938', '1951-1939', '1952-1940', '1953-1941', '1953-1942', '1954-1943', '1955-1944', '1956-1945', '1957-1946', '1958-1947', '1959-1948', '1960-1949', '1961-1950', '1962-1951', '1963-1952', '1963-1953', '1964-1954', '1965-1955', '1966-1956', '1967-1957', '1968-1957', '1969-1958', '1970-1959', '1971-1959', '1972-1960', '1973-1961', '1974-1962', '1975-1963', '1976-1964', '1977-1965', '1978-1966', '1979-1967', '1980-1968', '1981-1969', '1982-1970', '1983-1971', '1984-1972', '1985-1973', '1985-1974', '1986-1975', '1987-1976', '1988-1977', '1989-1978', '1990-1979', '1991-1979', '1992-1980', '1993-1981', '1994-1982', '1995-1982', '1996-1982', '1997-1983', '1998-1984', '1999-1985', '2000-1985', '2001-1985', '2002-1986', '2003-1987', '2004-1988', '2005-1989', '2006-1990', '2007-1991', '2008-1992', '2009-1993', '2010-1993', '2011-1994', '2012-1995', '2013-1996', '2014-1997', '2015-1998', '2016-1999', '2017-1999', '2018-2000', '2019-2001', '2020-2002', '2021-2003', '2022-2004', '2023-2005', '2024-2006', '2025-2007', '2026-2007', '2027-2008', '2028-2009', '2029-2010', '2030-2011', '2031-2011', '2032-2011', '2033-2012', '2034-2013', '2035-2014', '2036-2015', '2037-2016', '2038-2017', '2039-2018', '2040-2019', '2041-2020', '2042-2021', '2043-2022', '2044-2023', '2045-2024', '2046-2025', '2047-2026', '2048-2027', '2049-2028', '2050-2028', '2051-2029', '2052-2030', '2053-2031', '2054-2031', '2055-2031', '2056-2032', '2057-2033', '2058-2033', '2059-2034', '2060-2035', '2061-2036', '2062-2036', '2063-2036', '2064-2036', '2065-2037', '2066-2038', '2067-2039', '2068-2039', '2069-2040', '2070-2040', '2071-2041', '2072-2041', '2073-2042', '2074-2043', '2075-2043', '2076-2044', '2077-2045', '2077-2046', '2078-2047', '2079-2048', '2080-2049', '2081-2050', '2082-2051', '2083-2052', '2084-2053', '2085-2054', '2086-2055', '2087-2056', '2088-2057', '2089-2058', '2089-2059', '2090-2060', '2091-2060', '2092-2060', '2093-2061', '2094-2062', '2095-2063', '2096-2064', '2097-2065', '2098-2066', '2099-2067', '2100-2068', '2101-2068', '2102-2069', '2103-2069', '2104-2069', '2105-2069', '2106-2070', '2107-2071', '2108-2072', '2109-2073', '2110-2074', '2111-2075', '2112-2075', '2113-2076', '2114-2077', '2115-2078', '2116-2079', '2117-2080', '2118-2081', '2119-2082', '2120-2083', '2121-2083', '2122-2084', '2123-2085', '2124-2086', '2125-2087', '2126-2088', '2127-2089', '2128-2090', '2129-2091', '2130-2092', '2131-2093', '2132-2094', '2133-2095', '2134-2096', '2135-2097', '2136-2098', '2137-2099', '2138-2100', '2139-2101', '2140-2102', '2141-2103', '2142-2104', '2143-2104', '2144-2105', '2145-2106', '2146-2107', '2147-2108', '2148-2109', '2149-2110', '2150-2111', '2151-2111', '2152-2112', '2153-2113', '2154-2114', '2155-2115', '2156-2116', '2157-2117', '2158-2118', '2159-2119', '2160-2120', '2161-2121', '2162-2122', '2163-2123', '2164-2124', '2165-2125', '2166-2125', '2167-2126', '2168-2127', '2169-2128', '2170-2129', '2171-2129', '2172-2130', '2173-2130', '2174-2131', '2175-2132', '2176-2132', '2177-2132', '2178-2132', '2179-2132', '2180-2133', '2181-2134', '2182-2134', '2183-2135', '2184-2136', '2185-2136', '2186-2137', '2187-2137', '2188-2138', '2189-2139', '2190-2140', '2191-2140', '2192-2141', '2193-2142', '2194-2143', '2195-2144', '2196-2145', '2197-2146', '2198-2147', '2199-2147', '2200-2148', '2201-2149', '2202-2149', '2203-2150', '2204-2150', '2205-2151', '2206-2152', '2207-2153', '2208-2154', '2209-2155', '2210-2156', '2211-2156', '2212-2157', '2213-2158', '2214-2158', '2215-2159', '2216-2160', '2217-2161', '2218-2161', '2219-2162', '2220-2163', '2221-2164', '2222-2165', '2223-2166', '2224-2167', '2225-2167', '2226-2168', '2227-2168', '2228-2169', '2229-2170', '2230-2170', '2231-2171', '2232-2171', '2233-2171', '2234-2172', '2235-2172', '2236-2173', '2237-2174', '2238-2175', '2239-2176', '2240-2177', '2241-2178', '2242-2179', '2243-2180', '2244-2181', '2245-2182', '2246-2183', '2247-2184', '2248-2184', '2249-2185', '2250-2185', '2251-2185', '2252-2185', '2253-2186', '2254-2186', '2255-2187', '2256-2188', '2257-2189', '2258-2190', '2259-2191', '2260-2191', '2261-2191', '2262-2192', '2263-2193', '2264-2193', '2265-2193', '2266-2194', '2267-2194', '2268-2195', '2269-2196', '2270-2197', '2271-2198', '2272-2199', '2273-2200', '2273-2201', '2274-2202', '2275-2202', '2276-2203', '2277-2204', '2278-2205', '2279-2205', '2280-2206', '2281-2207', '2282-2208', '2283-2209', '2284-2210', '2285-2211', '2286-2212', '2287-2213', '2288-2214', '2289-2215', '2290-2216', '2290-2217', '2291-2218', '2292-2218', '2293-2219', '2294-2220', '2295-2221', '2296-2222', '2297-2223', '2298-2224', '2299-2225', '2300-2226', '2301-2227', '2302-2228', '2303-2229', '2304-2229', '2305-2230', '2306-2231', '2306-2232', '2307-2233', '2308-2234', '2309-2235', '2310-2236', '2310-2237', '2311-2238', '2312-2239', '2312-2240', '2313-2241', '2314-2242', '2315-2243', '2316-2244', '2316-2245', '2317-2246', '2318-2247', '2319-2248', '2320-2249', '2321-2250', '2322-2251', '2323-2252', '2324-2252', '2325-2252', '2326-2253', '2327-2254', '2328-2254', '2329-2254', '2330-2255', '2331-2256', '2332-2257', '2333-2258', '2334-2259', '2335-2260', '2336-2261', '2337-2262', '2338-2263', '2339-2264', '2340-2265', '2341-2266', '2342-2267', '2343-2268', '2344-2269', '2345-2270', '2346-2271', '2347-2272', '2348-2273', '2349-2273', '2350-2274', '2350-2275', '2351-2276', '2352-2277', '2353-2278', '2354-2279', '2354-2280', '2355-2281', '2356-2282', '2357-2283', '2358-2284', '2359-2285', '2360-2286', '2361-2287', '2362-2288', '2363-2289', '2364-2290', '2365-2291', '2366-2291', '2367-2292', '2368-2293', '2369-2294', '2370-2295', '2371-2296', '2372-2297', '2373-2298', '2374-2299', '2375-2300', '2376-2301', '2377-2302', '2378-2303', '2379-2304', '2380-2305', '2381-2306', '2382-2307', '2383-2308', '2384-2309', '2385-2309', '2386-2310', '2387-2311', '2388-2312', '2389-2313', '2390-2314', '2391-2314', '2392-2315', '2393-2316', '2394-2317', '2395-2317', '2396-2318', '2397-2319', '2398-2320', '2399-2321', '2400-2321', '2401-2322', '2402-2323', '2403-2324', '2404-2324', '2405-2325', '2406-2326', '2407-2327', '2407-2328', '2408-2329', '2409-2330', '2410-2331', '2411-2332', '2412-2333', '2413-2334', '2414-2335', '2415-2336', '2416-2337', '2417-2338', '2418-2339', '2419-2340', '2420-2341', '2421-2342', '2422-2343', '2423-2344', '2424-2345', '2425-2346', '2426-2347', '2426-2348', '2426-2349', '2427-2350', '2428-2351', '2429-2352', '2430-2353', '2430-2354', '2430-2355', '2430-2356', '2430-2357', '2431-2358', '2431-2359', '2432-2360', '2433-2361', '2434-2361', '2435-2362', '2436-2363', '2437-2364', '2438-2365', '2438-2366', '2439-2367', '2440-2368', '2441-2369', '2441-2370', '2442-2371', '2443-2372', '2444-2373', '2445-2374', '2445-2375', '2446-2376', '2446-2377', '2447-2378', '2448-2379', '2449-2380', '2450-2381', '2451-2382', '2452-2383', '2453-2384', '2454-2385', '2455-2386', '2456-2387', '2457-2387', '2458-2388', '2459-2389', '2460-2390', '2461-2391', '2462-2392', '2463-2392', '2464-2393', '2465-2394', '2466-2395', '2467-2396', '2468-2397', '2469-2398', '2470-2398', '2471-2399', '2472-2400', '2473-2401', '2474-2402', '2475-2403', '2476-2404', '2477-2405', '2478-2406', '2479-2406', '2480-2407', '2481-2408', '2482-2409', '2483-2410', '2484-2411', '2485-2412', '2486-2413', '2487-2414', '2487-2415', '2487-2416', '2488-2417', '2489-2418', '2490-2419', '2491-2420', '2492-2421', '2493-2422', '2494-2423', '2495-2424', '2496-2425', '2497-2426', '2498-2427', '2499-2428', '2500-2429', '2501-2430', '2502-2431', '2503-2432', '2504-2433', '2505-2434', '2506-2435', '2507-2436', '2508-2437', '2509-2438', '2510-2438', '2511-2439', '2512-2440', '2513-2441', '2514-2442', '2515-2443', '2516-2443', '2517-2444', '2518-2444', '2519-2445', '2520-2446', '2521-2446', '2522-2446', '2523-2446', '2524-2446', '2525-2447', '2525-2448', '2526-2449', '2527-2450', '2528-2451', '2529-2452', '2530-2453', '2531-2454', '2532-2455', '2533-2456', '2534-2457', '2535-2457', '2536-2458', '2537-2459', '2538-2459', '2539-2460', '2540-2460', '2541-2460', '2542-2461', '2543-2462', '2544-2463', '2545-2464', '2546-2465', '2547-2465', '2548-2466', '2549-2467', '2550-2467', '2551-2468', '2552-2469', '2553-2470', '2554-2471', '2555-2472', '2556-2473', '2557-2474', '2558-2475', '2559-2476', '2560-2477', '2561-2478', '2562-2479', '2563-2479', '2564-2480', '2565-2481', '2566-2482', '2567-2483', '2568-2484', '2569-2485', '2570-2486', '2571-2487', '2571-2488', '2572-2489', '2573-2490', '2574-2491', '2575-2492', '2576-2493', '2577-2493', '2578-2494', '2579-2495', '2580-2496', '2581-2497', '2582-2498', '2583-2499', '2584-2500', '2585-2501', '2586-2502', '2587-2503', '2588-2503', '2589-2504', '2590-2505', '2591-2506', '2592-2506', '2593-2507', '2594-2508', '2595-2509', '2596-2510', '2597-2511', '2598-2511', '2599-2512', '2600-2512', '2601-2512', '2602-2512', '2603-2513', '2604-2514', '2605-2515', '2606-2516', '2607-2517', '2608-2518', '2609-2519', '2610-2520', '2611-2521', '2612-2522', '2613-2523', '2614-2523', '2615-2523', '2616-2524', '2617-2525', '2618-2526', '2619-2527', '2620-2528', '2621-2529', '2622-2529', '2623-2530', '2624-2531', '2625-2532', '2626-2533', '2627-2534', '2628-2535', '2628-2536', '2628-2537', '2629-2538', '2630-2539', '2631-2540', '2632-2541', '2633-2542', '2634-2543', '2635-2544', '2635-2545', '2635-2546', '2635-2547', '2636-2548', '2636-2549', '2636-2550', '2637-2551', '2638-2552', '2639-2553', '2640-2554', '2641-2555', '2642-2556', '2643-2557', '2644-2558', '2645-2559', '2646-2559', '2647-2559', '2648-2560', '2649-2561', '2649-2562', '2650-2563', '2651-2564', '2652-2565', '2653-2566', '2654-2567', '2655-2568', '2656-2569', '2657-2570', '2658-2571', '2659-2572', '2660-2573', '2661-2574', '2662-2574', '2663-2575', '2664-2576', '2665-2577', '2666-2578', '2667-2579', '2668-2580', '2669-2580', '2670-2581', '2671-2582', '2672-2583', '2673-2584', '2674-2585', '2674-2586', '2675-2587', '2676-2588', '2677-2589', '2678-2590', '2679-2591', '2680-2592', '2681-2593', '2682-2594', '2683-2595', '2684-2596', '2685-2597', '2686-2598', '2687-2599', '2688-2600', '2689-2601', '2689-2602', '2690-2603', '2691-2604', '2692-2605', '2693-2606', '2694-2607', '2695-2608', '2696-2608', '2697-2608', '2698-2609']
'''

'''
container_sequence1 = 'KQHFYGGYSNNRVEKQNTPGWSTYTCENCQMCGVLMTNEQLLSMCPMQQYMGMRIKQVWFTWVRMWPYHESCMRHQGYQFVNHVNWVPYSKTIHEVNKFWDFMSFDCIKAYHFLYHSWRMHLWCFIGSPYIRCVICYTCPGWANWQHHWYNHEYSEWTMCQSGYNYVSCMMSDDMAHDKFIFHPHFGSGFIYHGAQEDSMRLQHMIRWSSFWEEIPVNVIIIEHYTIGVEHMGFCYCVNTCFDMICHWRNAGRLELCEQWLRGFKWPQHMKDEVSQSKEMIDLMSEMGWMHHETTQYHSKGFREFYIIQRWCSCQEVEPQWFPKMYPGQRPFQTPVCPLNYDPKRVNEPTPLHFCCQEFPVANWATMERKQTYVWDKNQSAVWMFGNLETIFTPIITPISETSAYNGYGSDMCSVIRAVAWCEVNVTNRMQDLRASWDMMLCGMEHSRIGHNTVDYLKWMQSYNTWPAFFNRPFLKYSPWKHMNHHSQWEYHHDFEYYYLGRANHYHRIESTNQEVWYGDEMMKGGVREWKDLPCGFKHGNVQDNLAKKHKWSSMVYSGEEEMLFTWQFGGSEKHKKWFYSAVWPKPWIDWRQLHGKEAFRHWRFYYYWVNPDFIYFQCVAGTPVYTMFNPGHICCNEFVCVYPICLQAGLYVETKGCIMQQQDRVKMSQNRDTYCPMGHHKWVPDDLNMNLAQWREANEQAVIDWINEHDEMDFKISAMYWWSKIFASQKMRIERKMQCLLFWPEPRVRVLFTQMYLMETHNHVGLKWHHYPLFSVKAKKLFVISYQMQHQLRLQNFAHAVGFVTSPLINAVDNNRNTWMTHMWLHGVKTKIRSGFVLIRASDFCVNYHLHDPMTMPLPDEWIKVTQLKDTMLQNHTAPYNYRKMACLNDDACPGDNIRAKPHPGETAHAHSYHSELREMQMIQGAWDFAANGMCSVYHVVTKSPRHTLPFVHCWQMKNGEWQGFWRWWSGITISMTHIMVPCIGLWLTTLIQDGNGDY'
fitting_sequence1 = 'FGFKHGNVQDSSMVYSEMLFTWQFGGNEKHKKWFGSWVWPKPWIDGPWRQLHGKEAFRHWYYWVNPDFIYHAMRVQVFRFQCVAGTPVYTMFNPGHIC'
alignment5, max_score_3 = scored_pairwise_fitting_alignment(container_sequence1, fitting_sequence1, BLOSUM62, BLOSUM62_KEY, -10)
print(max_score_3)
print(alignment5[0])
print(alignment5[1])
475
-GFKHGNVQDNLAKKHKWSSMVYSGEEEMLFTWQFGGSEKHKKWF-YS-AVWPKPWID--WRQLHGKEAFRHWRFYYYWVNPDFIY-------F--QCVAGTPVYTMFNPGHIC
FGFKHGNVQD--------SSMVYS-E--MLFTWQFGGNEKHKKWFG-SW-VWPKPWIDGPWRQLHGKEAFRHW--YY-WVNPDFIYHAMRVQVFRFQCVAGTPVYTMFNPGHIC'''


sequence5 = 'TCTCGGTGTCTTCCTTCCGGTGCTCGTTAGTCGGGTGAGCAAGCACTTCGCCTCCAGAACGCGAGGCGTGTTCTACTTTGGCGTGTCGGTACCCAGCATTGCAGACTCAGATCTCAGAACTGCCCCGAAACGATACCGTTTGAACTACCTTAAGACTGTGAATGGACCCAGCTCCGCGCCCTAGTGGCAACCTTACCTATCGGTATTGCTAGCTGCCCCCAGCGGCTATCCAATCAACCCCTGTGGGAATATCGAAGAGGGAAGTATTACACCCTCCGCCTGGCATTGCTGCAGACCTAATCGTGTCGCCAGGGTAAACCCTTATATCCTGGACCGAGTACTCGCCGGTTCGCAGTTCAGACATGCTGCGACTGATCGCTAATCTTTTGTGTCTGCCAGATTAAGTTTACGGTGTGAAAGCGCGCCGTACCCACGGAAATGAGCGTACGTTTACCGACACGCCCCCACGGAGCCAGCTAGCTCACGGGAGACTCGTCCAAAGTTCTAACACCGGATAGTGTCTGCACCCAGAGTACGTAGCTTCTCAAGCTGTCTTCGCTGATTACAGCGACTTTATAAAATCTGGCGAC-TGTTGTGGGCGGGCATGCGGTCGTTT-TTAACA-GCTGAATGGCGGGCCGACGCGGCTCT-ACGAGCAACGGTGTCGTAACTTAAAACACTTGGTTAATGGACTTC-CGAAGGGAG-GGAAATTTAGACCAACT-ATCAAGACATC-TAA-AGCAGACGCATTAGA-CGGTATGAATATGTTACGGCGTACACTGATGAA--GGACTAGCATCAATGAAACCCCAGGAACTG-AACCGGCTGTG-ATAT-AACCTCACAGTCGGCGA'
sequence6 = '                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     TCCGACAATATGTCGTGGGTGAGCACGCTGGCGGTGCTTGACACGCCAAATGGCGGGC-GACGGGACTCTCACGAG-A-CGGTC-C-TAACTTAACACACTTGGTTCGTGGAATTCGCGA-GGGAGAGGAGAATGAGGGTTAATCATCCATGCATCCTAAGAGCAGGTA-ATTAAAACGGTATGA-TAAATTACGGTGTT-A--GA-GAATCGGACTACGATAGATGAA-CCC-AGGCACAGCAACCTGCTGTGGACATGAATC-CGCAG-CG-CGACTAAAAACCTGTGTGATTCTTAGTGATCCCATGCCGAATGATCTACGCTTACAACACGTTTATCGCTTGACCTGTTCTTCGTACACTGCGCTTAGTCCGATTATGCTTGTGACGACATCCTGCGTTAGTAAGCCTGGCTATTTCGCATCAGTCAGTCTGTCAAATTAGCATACGAAACCCCATAAAACGGGAGAGAGGCGTGGACATGTACCTATGTTCTGGTCCTCGCGCTCCAAAACACAACATTGGCAGTAGCATGTAGACTATGCGCCGTCCGATAGGGCTCTCTATTCCAACAGGACGGTCCGCTAGGTAATATGATAGGAGGGCTAGGGCACGTGTCGCTAAATTAGTTCCAGAGCTTATCAGCGGCAGATCAGCGGTTATGAGATTTTGGGACCGTCACCTATAGTGTGTCAAGAGCACATTTCCCAGTAATTTGAAGCTGACGGGGATATCATGAAATCCTAATAGCGCATACGCAGTGCATCGGCGACTCCGTACTAATAGAGGATCTCGTGCGGATGAGTCCATTGTCCACGTACTACCCGAAATCCCTTGCTTGACTAATGCATAGAGGACCGGGTGCCGGCAGGCCGATTGAACGTTGTAGAGGTAGGAATCACCAAGTTGTATCGTGTATGACGCAAGGTGGTGGACTCTATC'
alignment6, max_score_4 = scored_pairwise_overlap_alignment(sequence5, sequence6, DNA_MATRIX, DNA_KEY, -2)
print(max_score_4)
print(alignment6[0])
print(alignment6[1])

'''200
GGAACAGAGAA-GACAT--CGGAAACT--A-GCGTGGGCCCCCAAAACTCACTGAGCGTGAACACACACGCATCGATCACGCAGTTGGTTGGCTGAATACACTGGGCGCAATTGGTC-AG-GTACACCTAGTCCTTTTGTTCAG--TG-TTTGACGAACCGTCCTCCCTTTTCAAATAGAACACCGCTAACTGTGTGCTCCTTCGGCTTAACTTGACAGCTTGCAAGCTAAGCCTGGAACAGCAACAAAGCGGCGACCCCTCGGGGGGTAC-AGGAGGGTC---TA--TA-TCTA-CCT----CGTGTT-CTGTTGTCTATTACGCCCCCAAATACCTAGCAGACTACAGT-CC-TGCAAGACGCCATATCAGCGAGTTGATATTTTGTGCCCGACCGCCCCACTACCAGAGATGGCGCCCCTCGCCA---CA-GACAGGGTGGAACTAGT-CACA--T-GTACA--ACAAATAATGGCGGGTTCCTGCGATGC-ATGCCCGACCTTTAGGCAAAAACCATATATCGGGCCCTTTGCCGTGGGCCATATTAACTCTGAATCCACCGTTCTATGTACGATATAGGCTAATTACAATTTGCTGCACCC-----ACATAGGTCTCGT--CC--CC----CAAGCCGTAATCAACATCTCA-ACGTACACCTGCCACGGCGCC-A-GAGCTAGACGACTG-TACATTTACAACTGGCTAACAGCTTG---A-G---CTCGGCCCGGAATAAAGGAGCTATCTTAGAGTGCCGCGCAGTCACATATTTTGGCGAGGTGTATACAGGTGAGCTGGCTGCGT--CCCACACA-AACCA--CGAGAATG-TTGGTGAGAGGC---A-G-GTCCTAGAGC-TAAC-TCCTAGAC---GGTGGAAGGATT-ATATAAGCTTGCCTCTAGTCATTTCGCGCCCAATTGTAGAATGATTTTTACACGATG
GGCACCTCCAATGACATAATGGTAAGTAGAGGCGTGAGCCCCAAAAACTCACTGAGCGTAAACACACAC-CAGCTTGCA-TC-GAT-GTTGG--GCATACACT-GGCCAAATTTGTCTCGCGT-CAGGT-GCCATTTTGGTCAGTCTGTTTTGACGAACCGTCCT--ATGGGAGACCTCTTGACCGCTAA-----TCCT-ATTAGTGTT---TTGACA---TGCAAGCTAA-----GAACCGCAACAAAGCGGC-A---AGAGGGGGGTACAAAGAGGATCTTTTACGGATTCAAGCATGAGACATGTTGGCGGCCTCCACTACACCCCCAAATA---AGTA-AC--CAGTCCCATGCAAGACGCCAT-T-TGC--G-TGACCTTGGTTTCCCGATCGCCCCACTACCAGAAATGCCGCCCCT-GATAGTTCACGACAGGGTGGAA-TAGTACCCAGTTAATATAGGAGATCTCATGGCGGGTTCCTGCGATGCAACGTCCGACC--T---------CGATATATCGGGCGCTGTGCCGTGGGCCAAATAAACTCT-TGTCCA-CGGTCCATGTAGGCT-T---C--GTTACAATTTGCCGAACCCTTGAGAGAGAGGTCTCGTACCCATCCTTTACAAGCCG-ATTCCGTACGTAAGACCTGCAATTGACATCGCGCCAATGCGCTAGA-AATAGATTC-CAAACAACTGGCTTACAGCTTGACTAGGTACCTCGGCCCGGAATGAAGGAGCTATC--AGAGTGCCGCGCAGGCACATATTTTGGCGAGGTGTATACAGGTGAGCTGG-TACTTAAACCACCCAGAAGAATGCGGGTGCGACAGATGAGAGGCAAAATGCGTCCTAGAGCAAAACGTTACAGACTATGGTGGAAGGATTAAGATATGCTTGCCTCTAGGCATTTCGCGCCCAGTTGTAG-ATAAATCTTA-TGGTTT
'''



#with open('dna_scores.txt', 'w') as scores_txt:
#    for i in range(len(seq1)):
#        for j in range(len(seq2)):
#            match_layer_score = str(scores1[f'{i}-{j}'])
#            delete_layer_score = str(scores1[f'D-{i}-{j}'])
#            insert_layer_score = str(scores1[f'I-{i}-{j}'])
#            scores_txt.write(f'  {i}-{j}|{match_layer_score}|D-{delete_layer_score}|I-{insert_layer_score}  ')
#        scores_txt.write('\n')
