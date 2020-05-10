import os
import pickle

# TODO: Add doc for class
class Simulation:

    # Constructor for Simulation
    def __init__(self, nodes, faces, cells):
        # Save the domain information
        self.nodes = nodes
        self.faces = faces
        self.cells = cells

        # Initialize the time of the simulation
        self.t = 0.

        # Initialize the restart number
        self.rest_num = 0
