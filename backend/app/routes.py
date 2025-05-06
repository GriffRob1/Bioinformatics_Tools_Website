from flask import Blueprint, jsonify, request, current_app
from .tools_list import tools_list
from bson import ObjectId


api = Blueprint('api', __name__)

@api.route('/tools-list')
def get_tools_list():
    mongo = current_app.extensions['pymongo']
    if (mongo.db.tools.count_documents({}) == 0):
        mongo.db.tools.insert_many(tools_list)
    return jsonify(mongo.db.tools.find({}))

@api.route('/update-tool-popularities', methods=['POST'])
def update_tool_popularities():
    idAndIsFavorited = request.get_json()
    mongo = current_app.extensions['pymongo']
    oid = idAndIsFavorited['oid']
    isFavorited = idAndIsFavorited['isFavorited']
    if isFavorited:
        mongo.db.tools.update_many({'_id': ObjectId(oid)}, {'$inc': {'popularity': -1}})
    else:
        mongo.db.tools.update_many({'_id': ObjectId(oid)}, {'$inc': {'popularity': 1}})
    return jsonify({'message': 'successfully updated tool popularity'})