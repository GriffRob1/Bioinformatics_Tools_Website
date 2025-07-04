from graph import Graph
from motif_finder import hamming_distance

def create_basic_alignment_graph(sequence1, sequence2):
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



def longest_path_to_alignment(sequence1, sequence2, longest_path):
    longest_path = id_to_tuple(longest_path)
    aligned_sequence1 = ''
    sequence1_index = 0
    aligned_sequence2 = ''
    sequence2_index = 0

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
        ids_as_tuples.append((int(indexes[0]), int(indexes[1])))
    return ids_as_tuples



def basic_pairwise_alignment(sequence1, sequence2):
    graph = create_basic_alignment_graph(sequence1, sequence2)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(sequence1, sequence2, longest_path), score



def create_scored_alignment_graph(sequence1, sequence2, scoring_matrix, scoring_key, indel_score):
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



def scored_pairwise_alignment(sequence1, sequence2, scoring_matrix, scoring_key, indel_score):
    graph = create_scored_alignment_graph(sequence1, sequence2, scoring_matrix, scoring_key, indel_score)
    longest_path, score = graph.longest_path_DAG('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(sequence1, sequence2, longest_path), score



def create_scored_local_alignment_graph(sequence1, sequence2, scoring_matrix, scoring_key, indel_score):
    graph = Graph()
    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            graph.add_node(f'{i}-{j}')

    for i in range(len(sequence1) + 1):
        for j in range(len(sequence2) + 1):
            if i != 0  and j != 0:
                graph.add_edge('0-0', f'{i}-{j}', weight=0)
            if i != len(sequence1)  and j != len(sequence2):
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

    return graph



def scored_pairwise_local_alignment(sequence1, sequence2, scoring_matrix, scoring_key, indel_score):
    graph = create_scored_alignment_graph(sequence1, sequence2, scoring_matrix, scoring_key, indel_score)
    longest_path, score = graph.longest_path_DAG_local_alignment('0-0', f'{len(sequence1)}-{len(sequence2)}')
    return longest_path_to_alignment(sequence1, sequence2, longest_path), score
'''
        g     a     g     g     a
    0-0 - 0-1   0-2   0-3   0-4   0-5
  a           \t      
    1-0   1-1   1-2   1-3   1-4   1-5
  g                 \t           
    2-0   2-1   2-2   2-3   2-4   2-5
  g                       \t
    3-0   3-1   3-2   3-3   3-4   3-5
  a                          |   
    4-0   4-1   4-2   4-3   4-4   4-5
  a                             \t 
    5-0   5-1   5-2   5-3   5-4   5-5
    
    - insertion
    | deletion
    \t match or mismatch
    
    -aggaa
    gagg-a
    
'''
dna_matrix = [[1, -1, -1, -1],
              [-1, 1, -1, -1],
              [-1, -1, 1, -1],
              [-1, -1, -1, 1]]

dna_list = ['A', 'C', 'G', 'T']

BLOSUM62_key = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G','H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X']

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

PAM250_key = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

PAM250 = [[2, -2, 0, 0, -3, 1, -1, -1, -1, -2, -1, 0, 1, 0, -2, 1, 1, 0, -6, -3],
          [-2, 12, -5, -5, -4, -3, -3, -2, -5, -6, -5, -4, -3, -5, -4, 0, -2, -2, -8, 0],
          [0, -5, 4, 3, -6, 1, 1, -2, 0, -4, -3, 2, -1, 2, -1, 0, 0, -2, -7, -4],
          [0, -5, 3, 4, -5, 0, 1, -2, 0, -3, -2, 1, -1, 2, -1, 0, 0, -2, -7, -4],
          [-3, -4, -6, -5, 9, -5, -2, 1, -5, 2, 0, -3, -5, -5, -4, -3, -3, -1, 0, 7],
          [1, -3, 1, 0, -5, 5, -2, -3, -2, -4, -3, 0, 0, -1, -3, 1, 0, -1, -7, -5],
          [-1, -3, 1, 1, -2, -2, 6, -2, 0, -2, -2, 2, 0, 3, 2, -1, -1, -2, -3, 0],
          [-1, -2, -2, -2, 1, -3, -2, 5, -2, 2, 2, -2, -2, -2, -2, -1, 0, 4, -5, -1],
          [-1, -5, 0, 0, -5, -2, 0, -2, 5, -3, 0, 1, -1, 1, 3, 0, 0, -2, -3, -4],
          [-2, -6, -4, -3, 2, -4, -2, 2, -3, 6, 4, -3, -3, -2, -3, -3, -2, 2, -2, -1],
          [-1, -5, -3, -2, 0, -3, -2, 2, 0, 4, 6, -2, -2, -1, 0, -2, -1, 2, -4, -2],
          [0, -4, 2, 1, -3, 0, 2, -2, 1, -3, -2, 2, 0, 1, 0, 1, 0, -2, -4, -2],
          [1, -3, -1, -1, -5, 0, 0, -2, -1, -3, -2, 0, 6, 0, 0, 1, 0, -1, -6, -5],
          [0, -5, 2, 2, -5, -1, 3, -2, 1, -2, -1, 1, 0, 4, 1, -1, -1, -2, -5, -4],
          [-2, -4, -1, -1, -4, -3, 2, -2, 3, -3, 0, 0, 0, 1, 6, 0, -1, -2, 2, -4],
          [1, 0, 0, 0, -3, 1, -1, -1, 0, -3, -2, 1, 1, -1, 0, 2, 1, -1, -2, -3],
          [1, -2, 0, 0, -3, 0, -1, 0, 0, -2, -1, 0, 0, -1, -1, 1, 3, 0, -5, -3],
          [0, -2, -2, -2, -1, -1, -2, 4, -2, 2, 2, -2, -1, -2, -2, -1, 0, 4, -6, -2],
          [-6, -8, -7, -7, 0, -7, -3, -5, -3, -2, -4, -4, -6, -5, 2, -2, -5, -6, 17, 0],
          [-3, 0, -4, -4, 7, -5, 0, -1, -4, -1, -2, -2, -5, -4, -4, -3, -3, -2, 0, 10]]

#import numpy as np
#data = np.loadtxt('PAM250.txt', dtype=int).tolist()
#print(data)
'''
alignment3 = scored_pairwise_local_alignment('CAAGGAACTCAATCGTCGAGTCCACGGGGGGCAGAACACGCTATATTTAATCTTGATGAGGAACGCAAATAACCATGGTTGCACGTGAGGATTTTCTTTAGTGAGTTGGGTTGCTTGGTAACTTATCCACTGCTATCTTAAGGGGGTTACTTCGGGATGAACGGCTTATGACAATCACAGTGAGGTCCGTCCCGGCCGATATGAGTTCTATGTTTTAACAGCGTCACCAGTGTCACGTACGGGGCCACCTCAGGCCCTGACCAGGGAATAGAGCGATTTGGGGACTTTCCCGGGTGATGTCTACCAGGAAGTTCGGTACCACTGACTTTGAATAATACTGTCAAAGGGGCTGCACCTTCCCGAGTTCGTCGTCATTACACAGCGCATATATTACACGTTAAGCCGTTTATCCGCATGTTATGCCAATTCGCGTCTTGCCAGGTGCCAACGAGCCTGATAAAGCAGTGGGTAGCGCCGGCACAGTATGTAGCAAGTTCCCCGCCGCGCGTTGAAAGCGTTACGTACAGGCGGCTAAGCGACGTTAAAATTGTCGCTTGCCTAACCCATCTCCCTGACACGGAACATAGCGAATAGTAGTCAACGGAGTTATGGTACAAAGCCTGAAAGCGACCTCAGACGAAGGGTCTGCCCGCAGGACGTGGGCTCTAATCCTCGGGGGCCTCGCCTACGTAGCACATCCCCAATAGCACTAAGAAGATGTGAACGAAACGCCGCTGTCGGATTCCAATTCTGAAATAGATAGTACCGGGTCCGAGGCGATGGAGGGTGGCGAAACCCCCATTTACGCATAGCGGTAACTTGGTCCCGGACTATTTATCAGTTGGTACCCTCGGCCCTGGTGGATGTGTTTTACTGGAATGATGTATAGACATCGCCTAC', 'CGCCCTACCTTTCGTGCCTATCAAAGTCGGTAGGCGCATTCCGAGCTCCGACTAGCGAGAAAACCAGCACAATGAAGTGCGCCTCCATATGTTAGTATACGGGTTGCCTCAGCTTTGGGCGGGCCTAAGGGCGGGATGAACGGCTTATTCCTGTGCAGTAGGTCGGTCCCGCCGATATGAGTTCTGGAGGTATTTAACAGCAGGGGCACGTGCACTGGGCCGCCTCAGGGGCTCACCATGACCAGGACATCGAGCGGGATCTCAACTTTGGGGACTTTCCCGGGTGATGTCTGCCACTGATGTTCGGTACCAGCGACCATACTGTGAAAAGGGGCTGCAAACATTATCTTCCCGAGTTGTCATTACACAGCATATATTACGCCGAAACTACGTCTCCGCATGTTTAACTCGCGTTTTGCCAGAGCCTAGCTAGGTTATAAAGGCGGCACAGTTTGTTCTAAAACCCGCCCTCGCCGAAGCGTTACGTACAGCGGCAAAGCGACGTTGAAATTGTCGATTGCCTATGGTTGACGCGGAACATAGCGTGACATAGTAGTGATACTCCAAGCGTGAAAGTGACCTAGCGGGTCTGAGCAGGACGTTGGCTCGAGCACGTGCAAGTTAGGGTTGTTCATGCCTGGGGAAAATAGTTTACGATTATTCGGCCATGCGTTAGCGTGCTGCTGTCCAGGCGCAGGGGCTTTCGAAAGTTTACTCCGTGATTGCGATTACCCTGAAGCCAGAGCTACACTTCCACGGACCGAACGCCGAATATCCGGATCCTGCCTGGCTTTCCGGCTTCGGAACTGAGAGAGGATCCGAGAGATAGCTGGTAA', dna_matrix, dna_list, -5)
print(alignment3[1])
print(alignment3[0][0])
print(alignment3[0][1])
print(hamming_distance(alignment3[0][0], alignment3[0][1]))'''

alignment4 = scored_pairwise_local_alignment('ENLHDMCNDHMIHGAVFKKQYKSGEPKEQIFWYPMIQFLISVCRSIVRNQKDADWPPHKHWERCYHNSTPTCCLFYGQHHHFDVAGMHKIAPMHLLVDGKTQMTILLRLQDHGAHWPYCIVTTRRAMGTEAEKSDYFNKTCKLRAVSWLYDQYLEPFNGYYSVCKNCGADRGECLMRYIAHPKMFCYYWNPNLRPPGHLMGNPRIWTMPAELFSYTFCTYYLDKLRWIQCNAMDWIHSYWETVDTNCMYITMFGLMQDNTVNTACCIPFIMGFIFSDMSNQGNSVEGMWPMMTNLQRARCNFHRSPKGEKYFTWGYDKAHGIGQAYPAKMMNQSMDQQTCRGAWYPLFGIWAVASFHSCLWESKKPHLAVDCAENMCAYHHMWLYSPGDGHVDARLDYTLTGKSKMPCWPKKQWADACDVRVCYAIFIKWLICNAEFGALMYQMCDKPVCFLMACKIQYTWMARCWADMWWNFRIEAQLDYDWVWRMLTNMMTPLCVHANYNQKYDEDSWCDSANSQFAWAWYVQWFNEWAQLTYMSEMQVWSFLSTGPIEIWMPTYEMGDASTNTDNLRWNYMFVLKVTHYVHFNDVHPYSRHNDPSLHIMWDMGSHCKFNVDMDYMIWMNHCRFCKLGTINYNNAPNEDTFNTRDWMVFSPYHSAVCEMYNWNGRPDGLVIDEAAPQSSMGINSDNRSSWNHRCNRYQCFDWYYMKDEYHERHNSSTFFYWGPIDIFQFKMPYCPRDKLGMGKWKQDNPDVVGCWGIGGVFFKCYKADLMVRSPAYHYGGCPHGVIQEPEVHRCRKCFLITHAFRGVINDWIFSMTQLDSSERHINIIDKKSHATCVYARAPLLPHVQFAVKTCWKWRDTVDVVAWQQYLITAAKGLAKIIMICWIYHLHEIHDQCPVYCGPIMPGNWQFNIEGQCCIMACKDTNQATAWHQLGSWGKHTEVRLAISNPHIHFSGSIVRLSVHCHHVSYHSIRALFETNMSCMSELQVTNLFCQFGRAANYKQSQEGKKCHNPVQQWHVLMWAMFTNEHPDEMKCAVMDFFPYPILNQNFADERHMQHFINVVTDYHQGEQYGYIVIEETEARIIVGCKSHDRLDDDLSSHVNHLYYLDLIRMDVMGDGMCMPYCPYCSQRSPPQISPMVHMGCFVWITFICKVQYIVCWIKTGPSLRYETYPWDSAKWLSWFAQINAWLVVIPCDFLMDTNNSCDHTIQLFALINMRFLSYLEWRMPAAIWAARFNCPETPFIPMYACFEMSCWLGNFPFVPHQDTTIYNTFNFSHHTFMKRKGHGVDCTKWTTREQAKCATASEENQTEEMNFDNYSFRIMRANFSLRFIFSYIERKLTSMNKTWWIGWSCRNPDQHSWLQMYFTHILCAEFEAFTAPPHYLHPMGEVMHIIHVHPRPPGMMDYLHHHMPQPFHQGDGDSIQSVIKSAAPYCFCQYGKIHYWACHIYDHDFAAMNTYITTEKAFHNYKYSQSPIMDIRGGPASTCSHRSNMAWNVRFGSNRKIYGYACMTCFSVIMMKMHRQKWKHGNCMQGAPRVCRLVAAGVFECRYHNDYQDVCTCNNRWTWEWHMMRHKDALLMSMPTTAIQPERSVISGEGENEVFAAKNQMVPQCPPTWETYRRNVNVQRSSCPRWTKVQTPHLSILVYNNRKQFMTNSSWMFAGTPTIDQMDDEPNYYKEWGNCTYAAMSQQKKVQIEFACCINHDNRFQLVAFHRAQAFNNIVYYVWFAKKRVQCPALDEHGAKMMMFYIIASYLFPNQIDFEWMKSPMEPCNMMDRKGRIRAQSKHFWMGRFYQPKAETARCSKWVNWMTIVAYWGSDHQVAKTPVRHWHSIRTVPMYVKFYCNIQPLSPEVAAPWNNIGMLFGNEHEMVARNALDKQDAPAFVAKFGHDYHAKEVFEKWAGTQPKTTNAGFISYWCDCYKRMAGCMSIWWGPVHRVPKFDWTCNLLCEWCRSMKQVDQTAIAYMTIMWCIRHFHTSKGWTHRMPAWEAIMCYVCELCDNLRERMRPFSLRLPCYEATCSDPCNEIWYFIAHPAAYEDQKTSVLKRGEDHVCSNVTHNKMKVQLSGVRHFEHYRTMQYKTRSCFKMCMRKMHFEMEQMDMIMWGVPDTYNMATWQVPHQMMSEHGHNDWRNSQHTPYHWIGRVFSGEWDQHNHAIDRRMQPNGRPDIPLSQRTEPDSLGSPAIMTLFKLISNKTNWTSQSGVNESYQFDHKNEDFNVEKSYPDIYNMQFFMQVLIRKIHKISFLEDYAVSPDWQMRKWSSHQQNNLRLPVRDRAADPSSSDFGFLDMPYDENMPNYNHYLFRQWADWMQPMEDEDSYRMQMYQNIVITEFPSFYIMCFHQFWSCVMWERTYGQLHLKNMFSRFPYGLQDHGFHPMHNSGNGPGLWARKRTYLGLQACPISVFYVRPIAPRWEGRNMTTTSHMWQAGQCWKTHPEPQPPLKGNDCSKEEWRTVFISGEEHQWWCKWRFRRPSMHSPCKYKLCWRLLDTRSKQSAPNLKGLAYNGKFLWDYEYHGLHNLEDDDFSLVIMQDRFMWQWWRRAGVMYDNTLPHFYRTWSSARLWHIVWACCMCQKWLPHNGRLGMWQMCFYNMPMMDNYEDAALLTGTPHLPIHMMMEIHFGWCGRLWHFGMQLVRMLPDLNWALYNMLCVTEWTMGHSHDHGWMNDRQVPRQCGWCTDQVMTPFLYFLNARKVHVLSCCD', 'LNQQIVLSGHKFRMKQELLSNNNSQMETMKPHEGACQDCLALERSNSTQCPQGWAIVWYDKLKMMMFAQDLVYEDMMSRRTNLMAWFQHSVCIWNIALNVGRNDEDSMDRSEHNVPSETQIQTMQLYMNHADLFWWPATFREWHGEKMYKTYDTWLWTIEINNCHMFAWVCSECPWWCDDWNMWDYADSYNIYVASFAFRNFGKGSAVDVMSWVPFMTIIKLRVTIVLYIQYTHHMTNVISHYTWMWMWRIMVCLMETQLIINIWFQQWVCSEACVCAAHQMMLTMNVHTVLQIVFMYMPCWPRVRPPMESFQDCRSGLCPRQTSFICRLADWRCRCVCPFHRTKQANCICCYDQSAITGPYRIIVATCGCMIKLAWVYNQESCANCTRRQWPGMAIETIIGKDFDANQVQIKRLLSNAPFYVIQEYNKAIHNYYVMCMFHLSAWWCEIIRTIQPPEDEHGMHNFWFGYPFSGNFEHLNEGAQKMVVQAPHYVYNRKWHHMYCWHYDYSLCRVAFRQKIKECQKVVWWLRPAITWAMETEQPWEYNCFKNFQNWLVWCTYPTHYHATCCYDHSCSTDVNSCNLITKFDCCINWMEIKHFTINRHPWIDNFIETHFPVHPERGMMLPCRESKDHRILERQSQQERDRHLFPMSLGCPNRKGLFYWCEKIPHEPCHGFQYNRNLAHTQWGYQVTGSNMQAHCTDWNWNAWQKHVIRIGMMNVHWMLATSVGRQRVDCHARLATTWYKTQWKWKEKEMWCMACHKTHNFYYDNQEFASPLWHCPWWPIESWWWNPKHGCNTCRGVGFRKEKQMKGRPHTGCHFQCMSMVSSERHNNIIDKKSHATCVYARAPLLSMHHFGYKHVQFAVKTCAPYTDRHQKWQCYDVDTTDVVAWQQYGITAAKGLAKIIISCPICWIYHLHYCRFVEHHDQCPVRDFLGQWGVSEGFQHIMPGNVRTSQFNIEGQMCIMACRDRNQATAWHQLGSWGKHIEVRLAVRLSVHCHNKIMFILVSYHSIRALFETNMSCMSELSVTNLFCQFGRAANPPSQEGKKCHNPVQDWHVLMWAMFTEESPDEMWMAMREDFVPYYMILERHMQHFINVVDKDYHQGQYGCIVIENIAARYCKSHDRLDDDEIYPWWYFSSHVNHLYEPGCIMLDLIRIPKLCPLDVMGDGMCMPYCPYCSQRNPPQISPMVHCFVWITFILKVQYFVCWIKTGPSLRKEKYAKWLSVFAQIIPCDFFRMDTNQALINMRQLSELTYRMPAAIWAQRFWCPETPFPMYASFEMKCWLGNFPFVPHQDTTIHHTFMVHPMRWPNYGHNNYWRSQAKCTPHRTASEENQTEAMNFDNYSKRIMRANFSLRFIFLTSMNKTWWDQHTWLQDMTFDVFQHILCAEFEAFTANTMSPPHYLHPMGEVMHIIHYPRPPGMMDVLHHHMPQPSWKKFWDRCNYCFCQYGKIPHWACHIGDHWFFAGNTYINTEKAFHNSKYSQSPIMYNIPFASTCHKMTGSHRSNMAWVVRFGKNRKIYGYTCMTCFSVITMKMHKRGCTLGFWWKHGNCMQSATAHIETWHCRLVAAFVFECRYHNYIGVMGYQDVCTCNNRWTWEQPTCHMMRHKQIHMWPDEASDFQKERMTNKVPRCPPWWPPTGCAVITYRRNVNHMQCQRSSCPKVQTPQDSILVYNNRKQLANAMFAGTPTPEEHYDRFDQMDDEPNYRHKERGNCTYPSKMGSANAMRQQGDCIRHDNEFQLVAFHRAQAFNNNVYNPWDNFQCPALDEHGAEMMMFYIIFSFQEMQYVEHAMCCNMRAQSKHFWMGRFYQPKAETRRCSKWVNWMTIVSNYWGSDHQVAKTPARTVYCDAQPLSMEECAPYVTLFYNEHEMVNALQKQDFGHDYHAKEVFDYIFLMKIGVTQLLPKDVTIISRPNNKAASDSIHKFTAKGLYVMMHTEAEGKFYEDKNQGWEIGFHIWGYQQFADSVVSVPTSISNWTVIAPQQKYYKVMRIIVSHKRDTSLILPNIMWIVWVNIFEDSMNCTLNHHNNPKAYAVPMGRIDYTIPIDIKPNVFFIMNCECEFDLEHHWCKALASCSDFLYADSLEGTGLGHRKSIVCCMYDYGECDCAEGPHTLGSYKRTAPHYTQVFERFMYPINMYMLTEAMSTGTVGHAIPDHDGECVQGSFDYNEIVMEHHIHRWRILDTKWKIRSALALKIGMYRWCCSKMCGFKPLHLGPPGSPPQYANPMEQYWLAWYYPCPYKQLLESYWMDANSQETPQVYNMGAEPTRFSIWFQNPYMENTWVHSHSERPCGYQAPYFKANCLRTCGSLWGFSTYPIYWGQSFYIQREMNQTWSNLWMTVMLLRYKGDTAECMVFCVWNHDEIFQTNKLQGNMALFMNGIWDTIYTELGHWWIHWPKGPPMCLPEPNMAWIMGSVEPSERQAGNKLTIQAYTFKSASYACHLVLYCRLCDHPASYRAFSVILMMSEDYYCQKDHKKFNEMWNIEMQNESRINVGPTPANNEFCADEMETELFCVGRDYAVNGDQQIGCVKYYHQNLTIAMVREANRECFQSIDWSRYWCRVGGQKLSPVKEWWIKVYEGQPCDMFCRTPEKEVHLFAQNPTTKWHVCSYLKDLARIVINSAMNASIAKDSSCIARQAKKPMPTYHPGEIPMGVPGDQGNYAFHKKERTHHEPWMENACFASNWCM', PAM250, PAM250_key, -5)
print(alignment4)
print(alignment4[1])
print(alignment4[0][0])
print(alignment4[0][1])
print(hamming_distance(alignment4[0][0], alignment4[0][1]))