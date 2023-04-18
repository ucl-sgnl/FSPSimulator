from utils.SpaceObject import Satellite

def run_simulation(policy):
    id = policy.get('id')
    print(f'Running simulation for policy {id}')

    # Create a Space Object
    space_object = Satellite()
    print(space_object.GUID)

    space_object.cospar_id = '2020-001A'

    print(space_object.cospar_id)
