import os
import json
from flask import Flask, request, jsonify
from main import run_simulation
from utils.UpdateCatalogue import update_catalogue_celestrak, update_catalogue_jsr

app = Flask(__name__)

@app.route('/')
def hello():
    return 'Hello, welcome to the UCL Future Space Populations API!'

# Use the schema when the final json version is ready
# schema = {
#     "type": "object",
#     "properties": {
#         "name": {"type": "string"},
#         "age": {"type": "number"}
#     },
#     "required": ["name", "age"]
# }

@app.route('/newsim', methods=['POST'])
def post_json():

    # validate JSON data against schema
    # try:
    #     jsonschema.validate(data, schema)
    # except jsonschema.exceptions.ValidationError as e:
    #     return {'message': f'Invalid JSON data: {e}'}, 400

    policy_data = request.json

    # Create a json file with GUID as the name, this will be used as an audit record.
    id = policy_data.get('id')
    if id:
        filename = os.path.join(os.getcwd(), 'src', 'data', 'policy', f'{id}.json')
        file_created = False
        if os.path.isfile(filename):
            with open(filename, 'w') as f:
                json.dump(policy_data, f)
                file_created = True
        else:
            with open(filename, 'w') as f:
                json.dump(policy_data, f)
                file_created = True
        
        if (file_created):
            run_simulation(policy_data)
        else:
            return jsonify({'error': 'File not created and simulation not found'}), 400
        
        return jsonify({'success': True}), 201
    else:
        return jsonify({'error': 'Missing id parameter'}), 400
    

# Get endpoint that will return the json file with the GUID as the name
@app.route('/getfile', methods=['GET'])
def get_json():
    fileType = request.args.get('type')
    id = request.args.get('id')
    filename = os.path.join(os.getcwd(), 'src', 'data', fileType, f'{id}.json')
    if os.path.isfile(filename):
        with open(filename, 'r') as f:
            data = json.load(f)
            return jsonify(data), 200
    else:
        return jsonify({'error': 'File not found'}), 400

    
# Get endpoint that will update the external catalogue
@app.route('/updatecatalogue', methods=['GET'])
def update_external_cat():
    if (request.args.get('type') == 'jsr'):
        update_catalogue_jsr()
        return jsonify({'success': True}), 201
    
    elif (request.args.get('type') == 'celestrak'):
        update_catalogue_celestrak()
        return jsonify({'success': True}), 201
     
    else:
        # return an old version of the celestrak
        return jsonify({'error': 'Invalid type parameter'}), 400

if __name__ == '__main__':
    app.run(debug=True) # use this for development, in production use app.run()