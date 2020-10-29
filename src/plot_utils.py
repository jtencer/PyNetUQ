import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')
from scipy import stats
from aria_component import aria_component


def plot_QoI_pdfs(components, pce_coeffs):
    # Calculate pce dimension (total number of exogenous ports in network)
    pce_dim = 0
    for comp in components:
        pce_dim += comp.get_num_exogenous_ports()
    # Generate germ samples
    germ_samples=np.random.normal(0,1, (10000,pce_dim))
    for i,comp in enumerate(components):
        my_pce_coeffs = pce_coeffs[i]
        coeffs_per_QoI = comp.get_num_pc_terms()
        times = comp.get_solution_times()
        num_times = len(times)
        if num_times == 1:
            for QoI_idx,QoI in enumerate(comp.QoIs):
                c_k = my_pce_coeffs[QoI_idx, ...]
                # Evaluate the PCE at the germ samples
                pce_evals=comp.evaluate_pce(c_k,germ_samples)
                #Peform kernel density estimation
                xpts_pce, PDF_data_pce= KDE(pce_evals)

                plt.figure(figsize=(10,10))
                plt.plot(xpts_pce, PDF_data_pce, linewidth=2, color='r', label='NISP full quadrature method')
                # Label Axes
                plt.xlabel("Temperature [K]", size=16)
                plt.ylabel("PDF", size=16)
                # Add title
                plt.suptitle(QoI, size=20)
                # Change tick size
                plt.tick_params(axis='both', labelsize=14)
                # Pad tick labels
                plt.gca().tick_params(pad=6)
                # Save figure
                fig_name="Component_%d_%s.pdf" % (i, QoI)
                plt.savefig(fig_name)
                plt.close('all')
        else:
            for QoI_idx,QoI in enumerate(comp.QoIs):
                mean = np.zeros((num_times))
                std = np.zeros((num_times))
                for timestep in range(num_times):
                    c_k = my_pce_coeffs[QoI_idx*num_times+timestep, ...]
                    mean[timestep] = c_k[0]
                    std[timestep] = c_k[1]

                plt.figure(figsize=(10,10))
                plt.plot(times, mean, linewidth=2, color='k', label='NISP full quadrature method')
                plt.fill_between(times, mean-std, mean+std, alpha=0.3, edgecolor='#000000', facecolor='#929591')
                # Label Axes
                plt.xlabel("Time [s]", size=16)
                plt.ylabel("Temperature [K]", size=16)
                # Add title
                plt.suptitle(QoI, size=20)
                # Change tick size
                plt.tick_params(axis='both', labelsize=14)
                # Pad tick labels
                plt.gca().tick_params(pad=6)
                # Save figure
                fig_name="Component_%d_%s.pdf" % (i, QoI)
                plt.savefig(fig_name)
                plt.close('all')

def get_plot_data(components, pce_coeffs):
    # Calculate pce dimension (total number of exogenous ports in network)
    pce_dim = 0
    for comp in components:
        pce_dim += comp.get_num_exogenous_ports()
    # Generate germ samples
    germ_samples=np.random.normal(0,1, (10000,pce_dim))
    returndata = []
    for i,comp in enumerate(components):
        my_pce_coeffs = pce_coeffs[i]
        coeffs_per_QoI = comp.get_num_pc_terms()
        times = comp.get_solution_times()
        num_times = len(times)
        if num_times == 1:
            xydata={}
            for QoI_idx,QoI in enumerate(comp.QoIs):
                c_k = my_pce_coeffs[QoI_idx, ...]
                # Evaluate the PCE at the germ samples
                pce_evals=comp.evaluate_pce(c_k,germ_samples)
                #Peform kernel density estimation
                xpts_pce, PDF_data_pce= KDE(pce_evals)
                xydata[QoI] = [xpts_pce, PDF_data_pce]
            returndata.append(xydata)
        else:
            xydata={}
            for QoI_idx,QoI in enumerate(comp.QoIs):
                mean = np.zeros((num_times))
                std = np.zeros((num_times))
                for timestep in range(num_times):
                    c_k = my_pce_coeffs[QoI_idx*num_times+timestep, ...]
                    mean[timestep] = c_k[0]
                    std[timestep] = 0
                    for i in range(1,coeffs_per_QoI):
                        std[timestep] += c_k[i] * c_k[i]
                    std[timestep] = math.sqrt(std[timestep])
                xydata[QoI] = [times, mean, std]
            returndata.append(xydata)
    return returndata

def likelihood_to_exceed_from_file(comp, threshold, filename):
    varnames = []
    for coeff_idx in range(comp.get_num_pc_terms()):
        varnames.append('PCE_%d' % coeff_idx)

    e = exodus.exodus(filename, mode='r')
    times = e.get_times()
    num_times = len(times)
    num_nodes = e.num_nodes()

    vals=[]
    for timestep in range(self.num_timesteps):
        tvals=[]
        for var in varnames:
            nodal_values = e.get_node_variable_values(var,timestep+1)
            tvals.append(nodal_values)
        val.append(tvals)
    e.close()

    pce_dim=1
    germ_samples=np.random.normal(0,1, (10000,pce_dim))

    c_k = np.zeros((len(varnames),))
    e = exodus.exodus(filename, mode='a')
    for timestep in range(num_times):
        for nid in range(num_nodes):
            for k in range(len(varnames)):
                c_k[k] = vals[timestep][k][nid]
            # Evaluate the PCE at the germ samples
            pce_evals=comp.evaluate_pce(c_k,germ_samples)

            #Peform kernel density estimation
            xpts_pce, PDF_data_pce= KDE(pce_evals)


def KDE(fcn_evals):
    """
    Performs kernel density estimation
    Input:
        fcn_evals: numpy array of evaluations of the forward model
    Output:
        xpts_pce: numpy array of points at which the PDF is estimated.
        PDF_data_pce: numpy array of estimated PDF values.
    """
    # Perform KDE on fcn_evals
    kern_pce=stats.kde.gaussian_kde(fcn_evals)
    # Generate points at which to evaluate the PDF
    xpts=np.linspace(fcn_evals.min(),fcn_evals.max(),200)
    # Evaluate the estimated PDF at these points
    PDF_data=kern_pce(xpts)
    return xpts, PDF_data


