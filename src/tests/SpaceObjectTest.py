def test_spaceobject_creation():
    """
    Test that when an object is instantiated, the correct attributes are set.
    In particular, the characteristic area and mass are imputed to some value if they are not specified or 0. 
    """
    object_1 = SpaceObject(  object_type='+', 
                                    payload_operational_status='Active',
                                    application="Unknown", 
                                    operator='Starink', 
                                    mass=0, 
                                    eccentricity=0.1, 
                                    inc=2*np.pi, 
                                    argp=12, 
                                    raan=252, 
                                    source = "Guatemala",
                                    launch_date='2012-10-12', 
                                    decay_date='2013-10-12',
                                    rso_name='bigboisat',
                                    perigee='1200',
                                    apogee='1000',
                                    tle="1 53544U 22101T   23122.20221856  .00001510  00000-0  11293-3 0  9999\n2 53544  53.2176  64.0292 0001100  79.8127 280.2989 15.08842383 38928")

