from flask import Flask, jsonify, request
from flask_restful import Resource, Api, reqparse
import subprocess
import os
import time

app = Flask(__name__)
api = Api(app)


import json


class calculation(Resource):
    '''Interface to get the log information from the optimization.'''

    def post(self):

        inputs = request.get_json(force=True)

        input_file = 'input.json'

        with open(input_file, 'w', encoding='utf-8') as f:

                json.dump(inputs, f, indent=4)

        out_file = 'output.json'
        if os.path.exists(out_file):
            os.remove(out_file)
            time.sleep(1)

        p=subprocess.Popen(["julia", "interface.jl", input_file, out_file])
        p.wait()
        # while not os.path.exists(out_file):
        #     print("Waiting for Julia")
        #     time.sleep(10)

        for x in range(5):
            try:
                with open(out_file) as json_file:
                    outputs = json.load(json_file)
                    break
            except Exception as ex:
                print("ERROR passing default value: {}".format(ex))
                outputs = {"demand_target": 25000}

        return outputs

api.add_resource(calculation, '/calculate')

if __name__ == '__main__':
    app.run(host='0.0.0.0',debug=False)
