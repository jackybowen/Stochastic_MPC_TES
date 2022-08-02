from flask import Flask, jsonify, request
from flask_restful import Resource, Api, reqparse
import subprocess
import os
from tcilc_predict import *
app = Flask(__name__)
api = Api(app)


import json


class calculation(Resource):
    '''Interface to get the log information from the optimization.'''

    def __init__(self, **kwargs):

        self.model_ckpt = kwargs.get('model_ckpt')
        self.torch_config = kwargs.get('torch_config')

    def put(self):

        inputs = request.get_json(force=True)
        if self.model_ckpt is not None and self.torch_config is not None:
            outputs = get_prediction(inputs, model_ckpt=self.model_ckpt, torch_config=self.torch_config)
        else:
            outputs = get_prediction(inputs)

        return outputs


def main(config):
    # FLASK REQUIREMENTS
    # ------------------
    api.add_resource(calculation, '/calculate', resource_class_kwargs = {"config": config})
    # --------------------------------------
    app.run(debug=False, host='0.0.0.0', port=5500)
    # --------------------------------------


if __name__ == '__main__':
    import sys
    main(sys.argv[1])
