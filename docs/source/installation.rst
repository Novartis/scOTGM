Setup and Installation
======================

This document outlines the steps required to set up and install the necessary environment and dependencies for the project. It is highly recommended to use a separate conda virtual environment for easy package management and to ensure there are no conflicts with existing packages.

Creating a Virtual Environment
------------------------------

A virtual environment is an isolated environment for Python projects. This ensures that each project has its own dependencies, irrespective of what dependencies every other project has. 

To create a conda virtual environment:

1. Create the Conda Environment:

   Create a new conda environment using the `requirements.txt` file. Replace `<myenv>` with your desired environment name:

   .. code-block:: bash

      conda create -n <myenv> --file requirements.txt python=3.10

2. Activate the Environment:

   Activate the newly created conda environment:

   .. code-block:: bash

      conda activate <myenv>

Installing the Project
----------------------

Once the environment is set up, you can proceed to install the project:

1. Install Project Dependencies:

   While the conda environment is active, install the project using pip:

   .. code-block:: bash

      pip install -e .

2. Setting up IPython Kernel:

   For projects involving IPython notebooks, it’s beneficial to set up an IPython kernel:

   .. code-block:: bash

      python -m ipykernel install --user --name=<myenv>

Dependency Version Check
------------------------

Inconsistencies in the environment can lead to unexpected behaviors or compatibility issues. Please ensure that the versions in your environment match those specified in the `requirements.txt` file.

Further Assistance
-------------------

For any issues or questions regarding the setup and installation process, please contact the project maintainers or refer to the project documentation.


