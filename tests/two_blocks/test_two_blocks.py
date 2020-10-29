import os
path = os.path.dirname(__file__)
import sys
sys.path.append(path+'/../../src')

import numpy as np
from pyNetUQ import network
from aria_component import aria_component

def test_solve_network():
    os.chdir(path)
    components = []
    components.append(aria_component("comp1.i", "comp1_output.e", np=1))
    components[0].add_QoI("centerpoint")
    components[0].add_required_files(["mesh.g", "comp1.i", "comp1_output.e"])
    components.append(aria_component("comp2.i", "comp2_output.e", np=1))
    components[1].add_QoI("centerpoint")
    components[1].add_required_files(["mesh.g", "comp2.i", "comp2_output.e"])
    components[0].add_endogenous_port("surface_1", "solution->temperature")
    components[1].add_endogenous_port("surface_2", "solution->temperature")
    components[0].add_exogenous_port("TL", 1000, 40)
    components[1].add_exogenous_port("TR", 500, 20)

    # Initialize NetUQ network
    my_network = network(components, "decomposed_problem")

    # Perform Simulations
    my_network.run_simulation(max_iter=5, tolerance=0.0001)

    # Compare solutioni
    centerpoint = np.loadtxt('centerpoint.gold.out')
    assert(np.all(centerpoint == np.around(my_network.QoI_pce_coeffs[0], decimals=3)))
    assert(np.all(centerpoint == np.around(my_network.QoI_pce_coeffs[1], decimals=3)))
    
