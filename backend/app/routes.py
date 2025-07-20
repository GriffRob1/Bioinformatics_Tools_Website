import sys
import traceback

from flask import Blueprint, jsonify, request, current_app
from .tools_list import tools_list
from bson import ObjectId

from .algorithms.DnaAbox_finder import approximate_pattern_occurrences
from .algorithms.DnaAbox_finder import frequent_words_mismatch_reverse_complement
from .algorithms.DnaAbox_finder import reverse_complement
from .algorithms.DnaAbox_finder import find_clumps
from .algorithms.DnaAbox_finder import find_skew
from .algorithms.DnaAbox_finder import hamming_distance
from .algorithms.DnaAbox_finder import neighbors
from .algorithms.DnaAbox_finder import find_potential_dnaa_boxes

from .algorithms.motif_finder import count_matrix
from .algorithms.motif_finder import profile_no_psuedocounts
from .algorithms.motif_finder import consensus
from .algorithms.motif_finder import score
from .algorithms.motif_finder import median_string_all_kmers
from .algorithms.motif_finder import greedy_motif_search
from .algorithms.motif_finder import iterate_randomized_motif_search
from .algorithms.motif_finder import iterate_gibbs_sampler

from .algorithms.debruijn_graph_genome_sequencer import kmer_composition
from .algorithms.debruijn_graph_genome_sequencer import path_to_genome
from .algorithms.debruijn_graph_genome_sequencer import read_pairs_path_to_genome
from .algorithms.debruijn_graph_genome_sequencer import debruijn_graph_from_kmers
from .algorithms.debruijn_graph_genome_sequencer import paired_debruijn_graph
from .algorithms.debruijn_graph_genome_sequencer import debruijn_assembler
from .algorithms.debruijn_graph_genome_sequencer import debruijn_assembler_read_pairs

from .algorithms.cyclopeptide_sequencer import transcribe
from .algorithms.cyclopeptide_sequencer import translate
from .algorithms.cyclopeptide_sequencer import linear_and_cyclo_peptide_spectrum
from .algorithms.cyclopeptide_sequencer import brute_force_cyclopeptide_sequencer
from .algorithms.cyclopeptide_sequencer import cyclopeptide_sequencing_ideal_spectrum
from .algorithms.cyclopeptide_sequencer import leaderboard_cyclopeptide_sequencer
from .algorithms.cyclopeptide_sequencer import convolution_cyclopeptide_sequencing
from .algorithms.cyclopeptide_sequencer import amino_acid_masses_to_letters

from .algorithms.sequence_alignment import scored_pairwise_global_alignment
from .algorithms.sequence_alignment import scored_pairwise_local_alignment
from .algorithms.sequence_alignment import scored_pairwise_fitting_alignment
from .algorithms.sequence_alignment import scored_pairwise_overlap_alignment
from .algorithms.sequence_alignment import scored_pairwise_global_alignment_affine_gap_penalty
from .algorithms.sequence_alignment import scored_pairwise_local_alignment_affine_gap_penalty
from .algorithms.sequence_alignment import scored_pairwise_fitting_alignment_affine_gap_penalty
from .algorithms.sequence_alignment import scored_pairwise_overlap_alignment_affine_gap_penalty
from .algorithms.sequence_alignment import greedy_multiple_alignment

api = Blueprint('api', __name__)



@api.route('/tools-list')
def get_tools_list():
    mongo = current_app.extensions['pymongo']
    mongo.db.tools.delete_many({})
    if mongo.db.tools.count_documents({}) == 0:
        mongo.db.tools.insert_many(tools_list)
    return jsonify(mongo.db.tools.find({}))






@api.route('/update-tool-popularities', methods=['POST'])
def update_tool_popularities():
    id_and_is_favorited = request.get_json()
    mongo = current_app.extensions['pymongo']
    oid = id_and_is_favorited['oid']
    is_favorited = id_and_is_favorited['isFavorited']

    if is_favorited: #value of is_favorited is from before toggle, since it only gets updated after react render
        mongo.db.tools.update_many({'_id': ObjectId(oid)}, {'$inc': {'popularity': -1}})
    else:
        mongo.db.tools.update_many({'_id': ObjectId(oid)}, {'$inc': {'popularity': 1}})
    return jsonify({'message': 'successfully updated tool popularity'})





@api.route('/tool-page/submit/<algorithm_name>', methods=['POST'])
def call_algorithm(algorithm_name):
    mongo = current_app.extensions['pymongo']
    current_tool = mongo.db.tools.find_one({'URL': algorithm_name})
    all_input_arguments = request.get_json()['inputText'].split()

    # handles tools with a list argument
    if current_tool['listArgumentIndex'] >= 0:
        index = current_tool['listArgumentIndex']
        arguments = all_input_arguments[:index] + [all_input_arguments[index:]]
    else:
        arguments = all_input_arguments

    # algorithm call and exception handling
    try:
        algorithm = globals()[algorithm_name]
        result = algorithm(*arguments)
    except Exception as e:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        result = f'An error occurred in calling the algorithm: {e}\n{traceback.format_exception(exc_type, exc_value, exc_traceback)}'


    # formatting as a string for display

    if type(result) is set:
        result = "\n".join(result)

    elif type(result) is tuple:
        temp_result = []
        for item in result:
            temp_result.append(str(item))
        result = "\n".join(temp_result)

    elif type(result) is dict:
        temp_result = []
        for key in result:
            temp_result.append(f'{key}: {', '.join(result[key])}')
        result = "\n".join(temp_result)

    elif type(result) is list:
        for i in range(0,len(result)):
            result[i] = str(result[i])
        result = "\n".join(result)

    return jsonify({'result': result})

