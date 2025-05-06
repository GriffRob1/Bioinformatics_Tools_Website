from flask import Flask
from .routes import api
from flask_cors import CORS
from flask_pymongo import PyMongo

mongo = PyMongo()

def create_app():
    app = Flask(__name__)
    app.config.from_pyfile('../config.py')
    mongo.init_app(app)
    app.extensions['pymongo'] = mongo
    CORS(app)
    app.register_blueprint(api)
    return app