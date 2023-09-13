Welcome to FSPSim's documentation!
==================================

Motivation and Significance
---------------------------

Human space endeavors have expanded dramatically over the past decades, with a significant increase in the number of satellites currently in orbit. 
The rapid growth and potential impacts of satellite launches have led to the development of various models and the need for a more accessible and adaptable modeling framework.
The FSPSim tool addresses this need by offering a transparent, configurable, and open-source solution. This tool empowers stakeholders to comprehend and address the potential impacts of satellite mega-constellations.
FSPSim aims to make source-sink evolutionary models more accessible to a broader audience. It allows stakeholders to utilize modeling insights crafted by the community, promoting a more inclusive approach to understanding the effects of various decisions on the space environment.

.. image:: FSPSim_system_diag_v1.png
   :alt: Overview of FSPSim's architecture
   :caption: Overview of FSPSim's architecture, showcasing how user inputs influence simulation parameters, the classes, and data sources they connect with, and the propagation mechanism.

Software Description
--------------------

FSPSim consists of various classes and modules designed to simulate the behavior and impact of space objects.

- **The SpaceCatalogue Class:** Manages all SpaceObject instances and fetches data from major space object catalogs.
- **The SpaceObject Class:** Represents an individual object in space with detailed attributes.
- **The LaunchModel Module:** Transforms user predictions into SpaceObjects.
- **The Propagator Module:** Drives the propagation of all SpaceObjects within a given SpaceCatalog.

Software Functionalities
------------------------

Users specify simulation parameters using a JSON file. The results are saved in a pickle file and can be further analyzed using popular data science libraries.

Illustrative Examples
---------------------

Based on provided input parameters, users can generate various visualizations to analyze the simulation data. 

.. image:: fspsim_example.png
   :alt: A sample altitude time-series plot
   :caption: A sample altitude time-series plot generated from a simulation.

Current and Future Developments
-------------------------------

Key areas of focus include refining the Launch Model and expanding the preliminary launch cost estimator.

.. toctree::
   :maxdepth: 6
   :caption: Contents:

   fspsim
   fspsim.utils

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
