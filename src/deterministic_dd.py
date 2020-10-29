import os
import sys
import shutil
import numpy as np
import multiprocessing

sys.path.append("/gpfs1/jtencer/UQTk_v3.0.4-install")
import PyUQTk.pce as uqtkpce
import PyUQTk.uqtkarray as uqtkarray
import PyUQTk.quad as uqtkquad
import PyUQTk.tools as uqtktools

from aria_component import aria_component

def run_ddd(components, working_directory, max_iter=30, tolerance=0.01, pce_order=3, exo_links=[], restart=False):
    """
    Perform deterministic domain decomposition wrapped in NISP
    """

    working_directory=os.path.abspath(working_directory)

    # Calculate pce dimension (total number of exogenous ports in network)
    pce_dim = 0
    for comp in components:
        pce_dim += comp.get_num_exogenous_ports()
    pce_dim -= len(exo_links)
    print("Network has %d exogenous ports." % pce_dim)

    # Generate PCE model
    pc_type="HG"
    pc_alpha=0.0
    pc_beta=1.0
    param = pce_order+1 # Parameter for quadrature point generation. Equal to number of quad points per dimension for full quadrature
    pc_model = uqtkpce.PCSet("NISP", pce_order,pce_dim,pc_type, pc_alpha,pc_beta)
    pc_model.SetQuadRule(pc_type, 'full', param)

    for comp in components:
        comp.pc_model = pc_model

    # Get Quadrature Points
    qdpts_uqtk = uqtkarray.dblArray2D()
    pc_model.GetQuadPoints(qdpts_uqtk)
    totquat = pc_model.GetNQuadPoints()
    qdpts = np.zeros((totquat,pce_dim))
    qdpts_uqtk.getnpdblArray(qdpts)

    # Create and populate working directories for each simulation
    if not restart:
        if os.path.isdir(working_directory):
            shutil.rmtree(working_directory)
        os.makedirs(working_directory)

    for quad_iter in range(totquat):
        # Create subdirectory for each quadrature point
        #   and copy required files into it
        newfolder=os.path.join(working_directory, "qdpt_%d" % quad_iter)
        if not os.path.isdir(newfolder):
            os.makedirs(newfolder)
            for c in components:
                for f in c.get_required_files():
                    shutil.copy2(f,newfolder)

    # Perform iterations
    for jacobi_iter in range(max_iter):
        print("Performing iteration %d..." % jacobi_iter)
        old_endo_data = []
        new_endo_data = []
        exo_port_idx = []
        exogenous_idx = 0
        for compidx,comp in enumerate(components):
            old_endo_data.append(comp.get_endogenous_data())
            running_jobs = []
            for quad_iter in range(totquat):
                targetfolder=os.path.join(working_directory, "qdpt_%d" % quad_iter)
                os.chdir(targetfolder)
                unique_ports=0
                my_ports = []
                for portidx,port in enumerate(comp.exogenous_ports):
                    found=False
                    for link in exo_links:
                        if compidx == link[2] and portidx ==link[3]:
                            port.value = port.mean + port.std * qdpts[quad_iter,exo_port_idx[link[0]][link[1]]]
                            my_ports.append(exo_port_idx[link[0]][link[1]])
                            found=True
                    if not found:
                        port.value = port.mean + port.std * qdpts[quad_iter,exogenous_idx + unique_ports]
                        my_ports.append(exogenous_idx + unique_ports)
                        unique_ports+=1

                job = comp.execute()
                running_jobs.append(job)

            exo_port_idx.append(my_ports)
            exogenous_idx += len(my_ports)

            # Wait for all simulations to complete
            print("Waiting for aria jobs to complete...")
            for job in running_jobs:
                job.communicate()

            new_endo_data.append(comp.get_endogenous_data())

        # Check for convergence
        maxdiff = 0
        for comp_idx in range(len(components)):
            for timestep, nodal_data in enumerate(old_endo_data[comp_idx]):
                for nid in nodal_data:
                    diff = abs(new_endo_data[comp_idx][timestep][nid] - old_endo_data[comp_idx][timestep][nid])
                    maxdiff = max(maxdiff, diff)
        print(maxdiff)
        if maxdiff < tolerance:
            break

    # Collect update PCE coefficients
    print("Collect Data")
    QoI_pce_coeffs = []
    for comp in components:
        QoI_data=np.zeros(( totquat, comp.get_num_QoI()*comp.get_num_timesteps() ))
        for quad_iter in range(totquat):
            targetfolder=os.path.join(working_directory, "qdpt_%d" % quad_iter)
            os.chdir(targetfolder)
            QoI_at_qdpt = comp.get_QoI_data()
            for i,QoI in enumerate(comp.QoIs):
                QoI_data[[quad_iter],(i*comp.get_num_timesteps()):((i+1)*comp.get_num_timesteps())] = QoI_at_qdpt[QoI]

        QoI_ck = np.zeros(( QoI_data.shape[1], pc_model.GetNumberPCTerms() ))
        for i in range(QoI_data.shape[1]):
            QoI_ck[i,...] = GalerkinProjection(pc_model, QoI_data[:,i])
        QoI_pce_coeffs.append(QoI_ck)

    return QoI_pce_coeffs

def GalerkinProjection(pc_model, f_evaluations):
    """
    Obtain PC coefficients by Galerkin Projection
    Input:
        f_evaluations: 1D numpy array (vector) with function to be projected,
                       evaluated at the quadrature points
    Output:
        Numpy array with PC coefficients
    """

    # Get parameters
    if len(f_evaluations.shape) > 1:
        print("This function can only project single variables for now")
        exit(1)

    npce = pc_model.GetNumberPCTerms()
    nqp = f_evaluations.shape[0]        # Number of quadrature points

    # UQTk array for PC coefficients for one variable
    c_k_1d_uqtk = uqtkarray.dblArray1D(npce,0.0)

    # UQTk array for function evaluations at quadrature points for that variable
    f_uqtk = uqtkarray.dblArray1D(nqp,0.0)
    for ipt in range(nqp):
        f_uqtk[ipt]=f_evaluations[ipt]

    # Galerkin Projection
    pc_model.GalerkProjection(f_uqtk,c_k_1d_uqtk)

    # Put PC coefficients in numpy array
    c_k = np.zeros(npce)
    for ip in range(npce):
        c_k[ip] = c_k_1d_uqtk[ip]

    # Return numpy array of PC coefficients
    return c_k


