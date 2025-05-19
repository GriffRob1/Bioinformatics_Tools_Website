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

    if current_tool['hasListArgument']:
        index = current_tool['listArgumentIndex']
        arguments = all_input_arguments[:index] + [all_input_arguments[index:]]
    else:
        arguments = all_input_arguments

    try:
        algorithm = globals()[algorithm_name]
        result = algorithm(*arguments)
    except Exception as e:
        result = 'An error occurred in calling the algorithm'

    return jsonify({'result': f'Result: {str(result)}'})

