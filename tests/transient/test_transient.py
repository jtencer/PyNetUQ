import os
path = os.path.dirname(__file__)
import sys
sys.path.append(path+'/../../src')

import numpy as np
from pyNetUQ import network
from aria_component import aria_component

def test_solve_transient_network():
    os.chdir(path)
    components=[]
    components.append(aria_component("comp1.i", "comp1_output.e", np=1))
    components.append(aria_component("comp2.i", "comp2_output.e", np=1))
    components[0].add_QoI("centerpoint")
    components[0].add_required_files(["mesh.g", "comp1.i", "comp1_output.e"])
    components[1].add_QoI("centerpoint")
    components[1].add_required_files(["mesh.g", "comp2.i", "comp2_output.e"])

    components[0].add_endogenous_port("surface_1", "solution->temperature")
    components[1].add_endogenous_port("surface_2", "solution->temperature")
    components[0].add_exogenous_port("TL", 1000, 40)
    components[1].add_exogenous_port("TR", 500, 20)

    for comp in components:
        comp.generate_pc_model(2, nord=1)

    # Initialize NetUQ network
    my_network = network(components, "decomposed_problem", merge_timelines=True, restart=False)

    # Perform Simulations
    my_network.run_simulation(max_iter=10, tolerance=0.0001)

    #np.savetxt(path+'/centerpoint1.gold.out', my_network.QoI_pce_coeffs[0], fmt="%.3f")
    #np.savetxt(path+'/centerpoint2.gold.out', my_network.QoI_pce_coeffs[1], fmt="%.3f")

    # Compare solutioni
    centerpoint = np.loadtxt('centerpoint1.gold.out')
    assert(np.all(centerpoint == np.around(my_network.QoI_pce_coeffs[0], decimals=3)))
    centerpoint = np.loadtxt('centerpoint2.gold.out')
    assert(np.all(centerpoint == np.around(my_network.QoI_pce_coeffs[1], decimals=3)))
