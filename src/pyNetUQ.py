import os
import shutil
import numpy as np
import multiprocessing

from aria_component import aria_component
import sync_times

import time
import logging

def run_point(network, idx, qdpt, workdir):
    return_code = network.execute_quadrature_point(idx, qdpt, workdir)
    assert(return_code == 0)

class network():

    def __init__(self, component_list, workdir=".", restart=False, merge_timelines=False):
        self.components = component_list
        self.working_directory = os.path.abspath(workdir)
        self.linked_endo_ports = []
        self.linked_exo_ports = []
        self.pce_dim = self.get_num_exo_ports()
        if restart:
            for comp in self.components:
                if comp.pc_model is None:
                    comp.generate_pc_model(self.pce_dim)
            print("Collect Data")
            self.endo_pce_coeffs, self.QoI_pce_coeffs = self.collect_data()
        else:
            print("Initialize PCE")
            self.endo_pce_coeffs = []
            self.QoI_pce_coeffs = []
            self.initialize_PCE()
            print("Create subfolders")
            self.create_folders()
        self.historyfile = None
        self.merge_timelines = merge_timelines

    def get_num_exo_ports(self):
        """
        Calculate total number of exogenous ports in network
        """
        count = 0
        for comp in self.components:
            count += comp.get_num_exogenous_ports()
        return count - len(self.linked_exo_ports)

    def get_num_endo_ports(self):
        """
        Calculate total number of endogenous ports in network
        """
        count = 0
        for comp in self.components:
            count += comp.get_num_endogenous_ports()
        return count

    def get_num_endo_pc_coefs(self):
        """
        Count total number of PCE coefficients in network
        """
        count = 0
        for comp in self.components:
            count += comp.get_num_pc_coeffs()
        return count

    def get_num_QoIs(self):
        """
        Count total number of QoIs in network
        """
        count = 0
        for comp in self.components:
            count += comp.get_num_QoI()
        return count

    def link_endogenous_ports(self, componentA, portA, componentB, portB):
        """
        Link endogenous ports between 2 components
        Useful for only evaluating subset of output BCs on components
        """
        self.linked_endo_ports.append([componentA, portA, componentB, portB])

    def link_exogenous_ports(self, componentA, portA, componentB, portB):
        """
        Link exogenous ports between 2 components
        Useful for shared material properties/BCs between components
        """
        if componentA > componentB:
            self.linked_exo_ports.append([componentB, portB, componentA, portA])
        else:
            self.linked_exo_ports.append([componentA, portA, componentB, portB])
        self.pce_dim -= 1

    def execute_quadrature_point(self, current_comp_idx, qdpt, workdir):
        """
        Execute forward problem for single component, quadrature point
        Evaluates endogenous outputs for other components and exogenous inputs
        """
        launchdir = os.getcwd()
        os.chdir(workdir)

        logging.basicConfig(filename='network.log', level=logging.DEBUG)

        start=time.time()
        # Evaluate PCEs for other components
        if len(self.linked_endo_ports) == 0:
            for idx,othercomp in enumerate(self.components):
                logging.info('Evaluating PCEs for component %d', idx)
                if idx!=current_comp_idx:
                    othercomp.set_all_endogenous_values(self.endo_pce_coeffs[idx], qdpt)
        else:
            for link in self.linked_endo_ports:
                compA = link[0]
                compB = link[2]
                if compA == current_comp_idx:
                    self.components[compB].set_port_endogenous_values(self.endo_pce_coeffs[compB], qdpt, link[3])
                elif compB == current_comp_idx:
                    self.components[compA].set_port_endogenous_values(self.endo_pce_coeffs[compA], qdpt, link[1])
        end=time.time()
        logging.info('Evaluating PCEs for other components took %f s.' % (end-start))

        start=time.time()
        # Convert Quadrature point in \xi_i to equivalent samples of input parameters
        # (taking advantage of the fact that inputs are assumed to be Gaussian)
        # This is equivalent to evaluating 1st order PC expansions for the input parameters.
        exogenous_idx = 0
        for compidx, comp in enumerate(self.components):
            for portidx,port in enumerate(comp.get_exogenous_ports()):
                found = False
                for link in self.linked_exo_ports:
                    if compidx == link[2] and portidx == link[3]:
                        port.value = self.components[link[0]].exogenous_ports[link[1]].value
                        found = True
                if not found:
                    port.value = port.mean + port.std * qdpt[0,exogenous_idx]
                    exogenous_idx += 1
        end=time.time()
        logging.info('Evaluating PCEs for exogenous ports took %f s.' % (end-start))

        start=time.time()
        # Solve forward problem on component
        running_sim = self.components[current_comp_idx].execute()
        end=time.time()
        logging.info('Starting simulation took %f s.' % (end-start))

        start=time.time()
        # Wait for simulation to complete
        running_sim.communicate()
        end=time.time()
        logging.info('Simulation execution took %f s.' % (end-start))

        # Check return code
        if not running_sim.returncode == 0:
            print("simulation failed in %s" % workdir)
            
        os.chdir(launchdir)
        return running_sim.returncode

    def create_folders(self):
        """
        Creates directory structure for working directory
        allows for simultaneous execution across components/quadrature points
        """
        # TODO: utilize symbolic links where possible to reduce overhead

        # Note that this will clear everything from working_directory.  Be careful to avoid data loss
        if os.path.isdir(self.working_directory):
            shutil.rmtree(self.working_directory)
        os.makedirs(self.working_directory)

        for compidx,comp in enumerate(self.components):
            compfolder=os.path.join(self.working_directory, "component_%d" % compidx)
            os.makedirs(compfolder)
            qdpts, totquat = comp.get_quadpts(self.pce_dim)
            for quad_iter in range(totquat):
                # Create subdirectory for each quadrature point
                #   and copy required files into it
                newfolder=os.path.join(compfolder, "qdpt_%d" % quad_iter)
                if not os.path.isdir(newfolder):
                    os.makedirs(newfolder)
                    for c in self.components:
                        for f in c.get_required_files():
                            shutil.copy2(f,newfolder)

    def consolidate_timesteps(self):
        """
        Navigate directory structure and generate composite timeline for data
        """
        for myidx,comp in enumerate(self.components):
            qdpts, totquat = comp.get_quadpts(self.pce_dim)
            folder_list = []
            
            for quad_iter in range(totquat):
                targetfolder=os.path.join(self.working_directory, "component_%d" % myidx)
                targetfolder=os.path.join(targetfolder, "qdpt_%d" % quad_iter)
                folder_list.append(targetfolder)

            filename = comp.get_output_filename()
            file_list = []
            for folder in folder_list:
                file_list.append(folder + "/" + filename)
            
            times = sync_times.get_combined_timeline(file_list)
            comp.set_solution_times(times)

            if self.historyfile:
                os.remove(self.historyfile)

    def collect_data(self):
        """
        Navigate directory structure and extract nodal-temporal values
        Perform projection to generate nodal-temporal PCE coefficients
        """
        launchdir = os.getcwd()

        endo_pce_coeffs = []
        QoI_pce_coeffs = []
        for myidx,comp in enumerate(self.components):
            qdpts, totquat = comp.get_quadpts(self.pce_dim)

            endo_data=np.zeros(( totquat, comp.get_num_endogenous_ports() ))
            QoI_data=np.zeros(( totquat, comp.get_num_QoI()*comp.get_num_timesteps() ))

            endo_pce_coeffs.append(np.zeros(( comp.get_num_endogenous_ports(), comp.get_num_pc_terms() )))
            QoI_pce_coeffs.append(np.zeros(( comp.get_num_QoI()*comp.get_num_timesteps(), comp.get_num_pc_terms() )))

            for quad_iter in range(totquat):
                targetfolder=os.path.join(self.working_directory, "component_%d" % myidx)
                targetfolder=os.path.join(targetfolder, "qdpt_%d" % quad_iter)
                os.chdir(targetfolder)

                logging.basicConfig(filename='network.log', level=logging.DEBUG)
                logging.info('Component %d, QuadPt %d' % (myidx, quad_iter))
                logging.info('endo_data.size: %d,%d' % (endo_data.shape[0], endo_data.shape[1]))
                logging.info('QoI_data.size: %d,%d' % (QoI_data.shape[0], QoI_data.shape[1]))

                # Collect endogenous data
                logging.info('Collect Endogenous Data')
                endo_at_qdpt = comp.get_endogenous_data()
                index = 0
                for nodal_data in endo_at_qdpt:
                    for nid in nodal_data:
                        endo_data[quad_iter, index] = nodal_data[nid]
                        index += 1

                # Collect QoI data
                QoI_at_qdpt = comp.get_QoI_data()

                for i,QoI in enumerate(comp.get_QoIs()):
                    QoI_data[[quad_iter],(i*comp.get_num_timesteps()):((i+1)*comp.get_num_timesteps())] = QoI_at_qdpt[QoI]

            # Perform PC projection via quadrature integration
            endo_ck = np.zeros(( endo_data.shape[1], comp.get_num_pc_terms() ))
            for i in range(endo_data.shape[1]):
                endo_ck[i,...] = comp.GalerkinProjection(endo_data[:,i])
            endo_pce_coeffs[myidx] = endo_ck

            QoI_ck = np.zeros(( QoI_data.shape[1], comp.get_num_pc_terms() ))
            for i in range(QoI_data.shape[1]):
                QoI_ck[i,...] = comp.GalerkinProjection(QoI_data[:,i])
            QoI_pce_coeffs[myidx] = QoI_ck

        os.chdir(launchdir)
        return endo_pce_coeffs, QoI_pce_coeffs

    def save_all_data(self, variable_name, comps=None):
        """
        Output nodal data to exodus
        """
        launchdir = os.getcwd()

        if comps is None:
            comps = range(len(self.components))

        for compidx in comps:
            exodus_filename = self.components[compidx].get_output_filename()
            pce_filename = "%s_pce.e" % exodus_filename.split('.')[:-1][0]

            qdpts, totquat = self.components[compidx].get_quadpts(self.pce_dim)
            all_data = None
            for quad_iter in range(totquat):
                targetfolder=os.path.join(self.working_directory, "component_%d" % compidx)
                targetfolder=os.path.join(targetfolder, "qdpt_%d" % quad_iter)
                os.chdir(targetfolder)

                # Collect nodal data
                data = self.components[compidx].get_all_data(variable_name, exodus_filename)

                if quad_iter == 0:
                    # size numpy array to hold data
                    datasize=0
                    for timestep, nodal_data in enumerate(data):
                        datasize += len(nodal_data)
                    all_data = np.zeros(( totquat, datasize ))

                for timestep, nodal_data in enumerate(data):
                    numnodes = len(nodal_data)
                    for i,node_val in enumerate(nodal_data):
                        index = timestep*numnodes + i
                        all_data[quad_iter, index] = node_val

            # Perform PC projection via quadrature integration
            nodal_ck = np.zeros(( all_data.shape[1], self.components[compidx].get_num_pc_terms() ))
            for i in range(all_data.shape[1]):
                nodal_ck[i,...] = self.components[compidx].GalerkinProjection(all_data[:,i])

            # Save nodal PCE
            os.chdir(launchdir)
            self.components[compidx].save_nodal_pce(nodal_ck, exodus_filename, pce_filename)

    def initialize_PCE(self):
        """
        Initialize all PCE coefficients
        """
        self.endo_pce_coeffs = []
        for comp in self.components:
            if comp.get_pc_model() is None:
                comp.generate_pc_model(self.pce_dim)

            my_endo_pce_coeffs = comp.initialize_PCE()

            self.endo_pce_coeffs.append(my_endo_pce_coeffs)
            self.QoI_pce_coeffs.append(np.zeros(( comp.get_num_QoI()*comp.get_num_timesteps(), comp.get_num_pc_terms() )))
            
    def compute_relative_residual(self, endo_pce_coeffs_new):
        """
        Compute relative stepsize ||x_k - x_(k-1)|| / ||x_k|| where ||.|| is the infinity norm
        """
        delta = []
        norm_factor = []
        for myidx in range(len(self.endo_pce_coeffs)):
            delta.append(np.linalg.norm(endo_pce_coeffs_new[myidx] - self.endo_pce_coeffs[myidx], np.inf))
            norm_factor.append(np.linalg.norm(endo_pce_coeffs_new[myidx], np.inf))

        return max(delta)/max(norm_factor)

    def anderson_acceleration(self, endo_pce_coeffs_new, mk=5):
        """
        Perform anderson accleration to modify step taken in DDUQ iteration
        """
        if self.historyfile is None:
            self.historyfile = os.path.join(self.working_directory, "PCE_history.csv")

        # serialize data
        compstrings = []
        for comp_coeffs in self.endo_pce_coeffs:
            flat = comp_coeffs.flatten()
            compstrings.append(','.join(['%f' % num for num in flat]))
        line = ','.join(compstrings)
        line += '\n'
        compstrings = []
        for comp_coeffs in endo_pce_coeffs_new:
            flat = comp_coeffs.flatten()
            compstrings.append(','.join(['%f' % num for num in flat]))
        line += ','.join(compstrings) + '\n'

        # Append to history file
        with open(self.historyfile, "a+") as f:
            f.write(line)

        # read in data
        data = np.genfromtxt(self.historyfile, delimiter=',', dtype=float)
        xi = data[::2].transpose()

        # Perform acceleration step
        if xi.shape[1] > 1:
            gi = data[1::2].transpose()
            fi = gi - xi
            dfi = np.diff(fi)
            dgi = np.diff(gi)

            # Solve linear least squares problem
            # TODO: use QR factorization, save off Q,R and update with scipy.linalg.qr_insert
            gamma = np.linalg.lstsq(dfi, fi[:,-1], rcond=None)[0]

            # Take step
            next_x = gi[:,-1] - np.matmul(dgi, gamma)

            # Unpack step
            endo_pce_coeffs_accel = []
            offset = 0
            for comp_coeffs in self.endo_pce_coeffs:
                myshape = comp_coeffs.shape
                mydata = next_x[offset:offset+np.product(myshape)].reshape(myshape)
                endo_pce_coeffs_accel.append(mydata)
                offset += np.product(myshape)

            # Delete first 2 lines from file if too many steps are recorded
            if xi.shape[1] > mk:
                os.system("tail -n +3 %s >> tmp.csv" % self.historyfile)
                os.system("mv tmp.csv %s" % self.historyfile)

            return endo_pce_coeffs_accel
        else:
            return endo_pce_coeffs_new
        


    def run_simulation(self, max_iter=30, tolerance=0.001):
        """
        Perform DDUQ simulation
        """

        # Perform iterations
        for jacobi_iter in range(max_iter):
            print("Performing iteration %d..." % jacobi_iter)
            args=[]
            for myidx,comp in enumerate(self.components):
                qdpts, totquat = comp.get_quadpts(self.pce_dim)
                for quad_iter in range(totquat):
                    targetfolder=os.path.join(self.working_directory, "component_%d" % myidx)
                    targetfolder=os.path.join(targetfolder, "qdpt_%d" % quad_iter)
                    args.append((self, myidx, qdpts[[quad_iter],:], targetfolder))

            # Begin simulations
            print("Starting aria jobs")
            tasks=[]
            for myargs in args:
                tasks.append(multiprocessing.Process(target=run_point, args=myargs))
                tasks[-1].start()

            # Wait for all simulations to complete
            print("Waiting for aria jobs to complete...")
            unfinished = len(tasks)
            print("%d/%d tasks to go..." % (unfinished,unfinished))
            while True:
                incomplete = [task.is_alive() for task in tasks]
                if incomplete.count(True) < unfinished:
                    unfinished = incomplete.count(True)
                    print("%d/%d tasks to go..." % (unfinished,len(incomplete)))
                if unfinished > 0:
                    time.sleep(2)
                else:
                    break
            
            for task in tasks:
                task.join()
                assert(task.exitcode==0)

            if jacobi_iter==0 and self.merge_timelines:
                self.consolidate_timesteps()

            # Collect update PCE coefficients
            print("Collect Data")
            endo_pce_coeffs_new, self.QoI_pce_coeffs = self.collect_data()
            
            if len(self.components)>1:
                print("Perform acceleration step")
                # Perform acceleration step
                if jacobi_iter==0:
                    next_step = endo_pce_coeffs_new
                    self.endo_pce_coeffs = endo_pce_coeffs_new
                else:
                    next_step = self.anderson_acceleration(endo_pce_coeffs_new)

                    print("Check for convergence")
                    # Check for convergence
                    stepsize = self.compute_relative_residual(endo_pce_coeffs_new)
                    self.endo_pce_coeffs = next_step
                    print(stepsize)
                    if stepsize < tolerance:
                        break

        return self.QoI_pce_coeffs, self.endo_pce_coeffs



