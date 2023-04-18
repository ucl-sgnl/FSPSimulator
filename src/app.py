# from flask import Flask, request
# import os
# import json

# app = Flask(__name__)

# @app.route('/policy', methods=['POST'])
# def policy():
#     # Get the JSON data from the request
#     data = request.get_json()

#     # Get the name of the file from the JSON data
#     filename = data.get('filename')

#     # If the filename is not specified, return an error
#     if not filename:
#         return {'error': 'No filename specified'}

#     # Build the path to the directory to save the file
#     directory = os.path.join('data', 'policy', 'subdirectory')

#     # If the directory does not exist, create it
#     os.makedirs(directory, exist_ok=True)

#     # Build the path to the file to save
#     filepath = os.path.join(directory, filename)

#     # Write the JSON data to the file
#     with open(filepath, 'w') as f:
#         json.dump(data, f)

#     return {'success': True}

import os
import json
from flask import Flask, request, jsonify

app = Flask(__name__)

@app.route('/')
def hello():
    return 'Hello, welcome to the future space populations!'

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

    data = request.json
    id = data.get('id')

    if id:
        filename = os.path.join(os.getcwd(), 'src', 'data', 'policy', f'{id}.json')
        print(filename)
        if os.path.isfile(filename):
            with open(filename, 'w') as f:
                json.dump(data, f)
        else:
            # Create the file
            with open(filename, 'w') as f:
                json.dump(data, f)
        
        
        return jsonify({'success': True}), 201
    
    
    else:
        return jsonify({'error': 'Missing id parameter'}), 400

if __name__ == '__main__':
    app.run()