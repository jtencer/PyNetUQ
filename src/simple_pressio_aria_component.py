import abstract_component

import os
import subprocess
import numpy as np
import exodus
from output_suppression import Suppressor
from output_suppression import suppress_stdout_stderr
import sync_times

import sys
sys.path.append("/gpfs1/jtencer/UQTk_v3.0.4-install")
import PyUQTk.pce as uqtkpce
import PyUQTk.uqtkarray as uqtkarray



class exogenous_port:
    """
    Uncertain parameters set using aprepro
    """
    def __init__(self, varname, nominal_value, uncertainty):
        self.varname = varname
        self.mean = nominal_value
        self.std = uncertainty
        self.value = nominal_value

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

    def __repr__(self):
        return "exogenous_port(%s, %s, %s)" % (self.varname, self.mean, self.std)

class endogenous_port:
    """
    Connections to other components
    """
    def __init__(self, sideset_name, field_name):
        self.ssname = sideset_name
        self.varname = field_name

    def __str__(self):
        return str(self.__class__) + ":" + str(self.__dict__)

    def __repr__(self):
        return "endogenous_port(%s, %s)" % (self.ssname, self.varname)

class simple_pressio_aria_component(abstract_component.Component):
    """
    Single component in DDUQ network with foward problem implemented as an aria model
    """

    def __init__(self, inputdeckfilename, outputfilename, np=1, output_times=None):
        """
        Initialize component object.  
        Run initial simulation to make sure output exodus file exists.
        """
        self.inputdeck = inputdeckfilename
        self.outfile   = outputfilename
        self.pce_file  = "%s.pce" % outputfilename.split('.')[:-1][0]
        self.endogenous_output_ports = []
        self.exogenous_ports  = []
        self.QoIs = []
        self.num_procs = np
        self.required_files = []
        self.pc_model = None
        self.num_endogenous_nodes = 0
        self.num_timesteps = 0
        self.output_times = output_times

        if not os.path.isfile(self.outfile):
            print("initializing component: %s" % inputdeckfilename)
            aprfile = "%s.apr" % self.inputdeck.split('.')[:-1][0]
            with open(os.devnull, 'w') as devnull:
                subprocess.check_output(["aprepro", "initialize=1", self.inputdeck, aprfile], stderr=devnull)
            subprocess.check_output(["sierra", "--pre", "-n", str(self.num_procs), "aria", "-i", aprfile])
            subprocess.check_output(["sierra", "--run", "-n", str(self.num_procs), "aria", "-i", aprfile])
            subprocess.check_output(["rm", aprfile])

        self.required_files.append(inputdeckfilename)
        self.required_files.append(outputfilename)

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

    def __call__(self):
        return self.execute()

    def execute(self):
        """
        Perform forward problem
        """
        aprfile = "%s.apr" % self.inputdeck.split('.')[:-1][0]
        apr_command=["aprepro"]
        for port in self.exogenous_ports:
            apr_command.append("%s=%f" % (port.varname, port.value))
        apr_command.append(self.inputdeck)
        apr_command.append(aprfile)
        with open(os.devnull, 'w') as devnull:
            subprocess.check_output(apr_command, stderr=devnull)

        if self.num_procs > 1:
            subprocess.check_output(["sierra", "--pre", "-n", str(self.num_procs), "aria", "-i", aprfile])
        with open(os.devnull, 'w') as devnull:
            p = subprocess.Popen(["sierra", "--run", "-n", str(self.num_procs), "aria", "-i", aprfile], stdout=devnull)
        return p

    def add_endogenous_port(self, ssname, varname):
        """
        Add an endogenous port between two aria_component instances
        Port is specified on the *sending* component
        """
        with Suppressor():
            e = exodus.exodus(self.outfile, mode='r')
            ssnames = e.get_side_set_names()
            e.close()
        assert ssname in ssnames, "%s not a sideset in %s." % (ssname, self.outfile)
        my_port = endogenous_port(ssname,varname)
        self.endogenous_output_ports.append(my_port)

    def add_exogenous_port(self, varname, nominal_value, uncertainty):
        """
        Specify aprepro variable that will be used as an exogenous input
        """
        self.exogenous_ports.append(exogenous_port(varname, nominal_value, uncertainty))

    def add_QoI(self, varname):
        """
        Specify global variables from exodus output to be treated at QoIs
        """
        with Suppressor():
            e = exodus.exodus(self.outfile, mode='r')
            gnames = e.get_global_variable_names()
            e.close()
        assert varname in gnames, "%s not a global variable in %s." % (varname, self.outfile)
        self.QoIs.append(varname)

    def get_num_endogenous_ports(self):
        nodes = self.get_num_endogenous_nodes()
        steps = self.get_num_timesteps()
        return nodes * steps

    def get_num_endogenous_nodes(self):
        if (self.num_endogenous_nodes == 0) and (len(self.endogenous_output_ports) > 0):
            self.get_endogenous_data()
        return self.num_endogenous_nodes

    def get_num_timesteps(self):
        if self.num_timesteps == 0:
            self.num_timesteps = len(self.get_solution_times())
        return self.num_timesteps

    def set_solution_times(self, times):
        self.output_times = times
        self.num_timesteps = len(times)

    def get_solution_times(self):
        if self.output_times:
            print("using self.output_times")
            times = self.output_times
        else:
            print("no output times found, reading from exodus")
            with Suppressor():
                e = exodus.exodus(self.outfile, mode='r')
                times = e.get_times()
                e.close()
        return times

    def get_num_exogenous_ports(self):
        return len(self.exogenous_ports)

    def get_output_filename(self):
        return self.outfile

    def get_exogenous_ports(self):
        return self.exogenous_ports

    def get_num_QoI(self):
        return len(self.QoIs)

    def get_QoIs(self):
        return self.QoIs

    def get_num_pc_terms(self):
        return self.pc_model.GetNumberPCTerms()

    def get_pc_model(self):
        return self.pc_model

    def get_num_pc_coeffs(self):
        return self.get_num_pc_terms()*self.get_num_endogenous_ports()

    def get_endogenous_data(self):
        """
        Retreive my output data at endogenous nodes
        """
        if self.output_times:
            sync_times.interpolate_to_timeline(self.outfile, self.outfile+".tmp", self.output_times)
            os.rename(self.outfile+".tmp", self.outfile)

        with Suppressor():
            e = exodus.exodus(self.outfile, mode='r')
        ss_ids = e.get_side_set_ids()
        ss_names = e.get_side_set_names()
        dictionary = dict(zip(ss_names, ss_ids))

        # Get list of time steps for which to provide data
        times = e.get_times()
        self.num_timesteps = len(times)

        vals=[]
        for timestep in range(self.num_timesteps):
            self.num_endogenous_nodes = 0
            for port in self.endogenous_output_ports:
                endogenous_vals = {}
                ssid = e.get_side_set_node_list(dictionary[port.ssname])[1]
                nodal_values = e.get_node_variable_values(port.varname,timestep+1)
                side_set_unique_node_ids = set(ssid)
                for nid in side_set_unique_node_ids:
                    endogenous_vals[nid] = nodal_values[nid-1]
                vals.append(endogenous_vals)
                self.num_endogenous_nodes += len(endogenous_vals)

        with Suppressor():
            e.close()
        return vals

    def get_all_data(self, varname, filename=None):
        """
        Retreive my output data at all nodes
        """
        if filename==None:
            filename = self.outfile
        with Suppressor():
            e = exodus.exodus(filename, mode='r')

            # Get list of time steps for which to provide data
            times = e.get_times()
            self.num_timesteps = len(times)

            vals=[]
            for timestep in range(self.num_timesteps):
                nodal_values = e.get_node_variable_values(varname,timestep+1)
                vals.append(nodal_values)

            e.close()
        return vals

    def get_QoI_data(self):
        with Suppressor():
            e = exodus.exodus(self.outfile, mode='r')
        QoI_vals = {}
        for QoI in self.QoIs:
            QoI_vals[QoI] = e.get_global_variable_values(QoI)
        with Suppressor():
            e.close()
        return QoI_vals

    def get_required_files(self):
        return self.required_files

    def add_required_files(self, files):
        for f in files:
            if f not in self.required_files:
                self.required_files.append(f)

    def generate_pc_model(self, pce_dim, nord=3, pc_type="HG", pc_alpha=0.0, pc_beta=1.0, quadtype='full'):
        """
        Wrapper for uqtk PCSet with default values
        """
        param = nord+1 # Parameter for quadrature point generation. Equal to number of quad points per dimension for full quadrature
        self.pc_model = uqtkpce.PCSet("NISP", nord,pce_dim,pc_type, pc_alpha,pc_beta)
        self.pc_model.SetQuadRule(pc_type, quadtype, param)

    def initialize_PCE(self):
        if os.path.isfile(self.pce_file):
            # Read initial PCE values from exodus file
            my_endo_pce_coeffs = np.zeros(( self.get_num_endogenous_ports(), self.get_num_pc_terms() ))

            varnames = []
            for coeff_idx in range(self.get_num_pc_terms()):
                varnames.append('PCE_%d' % coeff_idx)

            e = exodus.exodus(self.pce_file, mode='r')
            ss_ids = e.get_side_set_ids()
            ss_names = e.get_side_set_names()
            dictionary = dict(zip(ss_names, ss_ids))

            # Get list of nodes for which to provide data
            #TODO: This likely broken from port change
            all_side_set_node_ids = []
            for port in self.endogenous_output_ports:
                side_set_node_ids = e.get_side_set_node_list(dictionary[port.ssname])[1]
                all_side_set_node_ids.append(side_set_node_ids)

            endo_map = self.get_endogenous_data()
            for timestep, node_map in enumerate(endo_map):
               print("timestep: %d" % timestep)
               for coeff_idx in range(self.get_num_pc_terms()):
                   varname = varnames[coeff_idx]
                   nodal_values = e.get_node_variable_values(varname,1)
                   for ssid in all_side_set_node_ids:
                       side_set_unique_node_ids = set(ssid)
                       for nid in side_set_unique_node_ids:
                           index = timestep*self.num_endogenous_nodes + node_map.keys().index(nid)
                           my_endo_pce_coeffs[index,coeff_idx] = nodal_values[nid-1]

            e.close()

        else:
            endo_init = self.get_endogenous_data()
            my_endo_pce_coeffs = np.zeros(( self.get_num_endogenous_ports(), self.get_num_pc_terms() ))

            index = 0
            for timestep in range(self.num_timesteps):
                for portid,port in enumerate(self.endogenous_output_ports):
                    nodal_data = endo_init[timestep*len(self.endogenous_output_ports) + portid]
                    for nid in nodal_data:
                        my_endo_pce_coeffs[index,0] = nodal_data[nid]
                        index += 1

        return my_endo_pce_coeffs



    def GalerkinProjection(self,f_evaluations):
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

        npce = self.pc_model.GetNumberPCTerms()
        nqp = f_evaluations.shape[0]        # Number of quadrature points

        # UQTk array for PC coefficients for one variable
        c_k_1d_uqtk = uqtkarray.dblArray1D(npce,0.0)

        # UQTk array for function evaluations at quadrature points for that variable
        f_uqtk = uqtkarray.dblArray1D(nqp,0.0)
        for ipt in range(nqp):
            f_uqtk[ipt]=f_evaluations[ipt]

        # Galerkin Projection
        self.pc_model.GalerkProjection(f_uqtk,c_k_1d_uqtk)

        # Put PC coefficients in numpy array
        c_k = np.zeros(npce)
        for ip in range(npce):
            c_k[ip] = c_k_1d_uqtk[ip]

        # Return numpy array of PC coefficients
        return c_k

    def evaluate_pce(self, pc_coeffs,germ_samples):
        """
        Evaluate PCE at a set of samples of the germ of this PCE
        Input:
            pc_coeffs: 1D numpy array with PC coefficients of the RVs to be evaluated.
                       Each column corresponds to one RV.
            germ_samples: numpy array with samples of the PCE germ at which the RVs
                          are to be evaluated. Each line is one sample. The number
                          of columns is the number of RVs.

        Output:
            Numpy array with PCE evaluations
        """

        # Get data set dimensions etc.
        n_test_samples = germ_samples.shape[0]
        ndim = germ_samples.shape[1]
        npce = self.pc_model.GetNumberPCTerms()

        # Put PC germ samples in a UQTk array
        std_samples_uqtk = uqtkarray.dblArray2D(n_test_samples, ndim)
        std_samples_uqtk.setnpdblArray(np.asfortranarray(germ_samples))

        # Numpy array to store all RVs evaluated from sampled PCEs
        rvs_sampled = np.zeros(n_test_samples)

        # Evaluate PCE for RVs in each dimension
        # Create and fill UQTk array for PC coefficients
        c_k_1d_uqtk = uqtkarray.dblArray1D(npce,0.0)
        for ip in range(npce):
            c_k_1d_uqtk[ip] = pc_coeffs[ip]

        # Create UQTk array to store outputs in
        rv_from_pce_uqtk = uqtkarray.dblArray1D(n_test_samples,0.0)

        # Evaluate the PCEs for each input RV at those random samples
        self.pc_model.EvalPCAtCustPoints(rv_from_pce_uqtk,std_samples_uqtk,c_k_1d_uqtk)

        # Put evaluated samples in numpy array
        for isamp in range(n_test_samples):
            rvs_sampled[isamp] = rv_from_pce_uqtk[isamp]

        # Return numpy array of PCE evaluations
        return rvs_sampled

    def save_nodal_pce(self, pc_coeffs, meshfilename, outputfilename):
        if os.path.isfile(outputfilename): os.remove(outputfilename)
        print("Save nodal PCE %s" % outputfilename)
        times = self.get_solution_times()

        varnames = []
        for coeff_idx in range(pc_coeffs.shape[1]):
            varnames.append('PCE_%d' % coeff_idx)

        e = exodus.copy_mesh(meshfilename, outputfilename)
        e.close()
        e = exodus.exodus(outputfilename, mode='a')
        exodus.add_variables(e, nodal_vars=varnames)

        numnodes = pc_coeffs.shape[0]/self.num_timesteps

        for timestep in range(self.num_timesteps):
            for coeff_idx in range(pc_coeffs.shape[1]):
                varname = varnames[coeff_idx]
                nodal_values = e.get_node_variable_values(varname,1)
                for nidx in range(numnodes):
                    index = timestep*numnodes + nidx
                    nodal_values[nidx] = pc_coeffs[index,coeff_idx]
                e.put_node_variable_values(varname,timestep+1,nodal_values)
                e.put_time(timestep+1,times[timestep])

        e.close()

    def set_all_endogenous_values(self, pc_coeffs, germ):
        """
        Sample polynomial chaos expansion for endogenous values at germ
        Assign those values to the nodes on the exodus mesh
        """
        endo_map = self.get_endogenous_data()

        with Suppressor():
            e = exodus.exodus(self.outfile, mode='a')
        ss_ids = e.get_side_set_ids()
        ss_names = e.get_side_set_names()
        dictionary = dict(zip(ss_names, ss_ids))

        index = 0
        for timestep in range(self.num_timesteps):
                for portid,port in enumerate(self.endogenous_output_ports):
                    node_map = endo_map[timestep*len(self.endogenous_output_ports) + portid]
                    nodal_values = e.get_node_variable_values(port.varname,timestep+1)

                    ssid = e.get_side_set_node_list(dictionary[port.ssname])[1]
                    side_set_unique_node_ids = set(ssid)
                    for nid in side_set_unique_node_ids:
                        idx = index + node_map.keys().index(nid)
                        endo_val = self.evaluate_pce(pc_coeffs[idx,...],germ)
                        nodal_values[nid-1] = endo_val
                    index += len(side_set_unique_node_ids)
                    e.put_node_variable_values(port.varname,timestep+1,nodal_values)

        with Suppressor():
            e.close()


    def set_port_endogenous_values(self, pc_coeffs, germ, portid):
        """
        Sample polynomial chaos expansion for endogenous values at germ
        Assign those values to the nodes on the exodus mesh
        """

        endo_map = self.get_endogenous_data()

        with Suppressor():
            e = exodus.exodus(self.outfile, mode='a')
        ss_ids = e.get_side_set_ids()
        ss_names = e.get_side_set_names()
        dictionary = dict(zip(ss_names, ss_ids))

        # Get list of nodes
        ssname = self.endogenous_output_ports[portid].ssname
        varname = self.endogenous_output_ports[portid].varname
        side_set_node_ids = e.get_side_set_node_list(dictionary[ssname])[1]

        for timestep, node_map in enumerate(endo_map):
            nodal_values = e.get_node_variable_values(varname,timestep+1)
            ssid = e.get_side_set_node_list(dictionary[ssname])[1]
            side_set_unique_node_ids = set(ssid)
            for nid in side_set_unique_node_ids:
                index = timestep*self.num_endogenous_nodes + node_map.keys().index(nid)
                endo_val = self.evaluate_pce(pc_coeffs[index,...],germ)
                nodal_values[nid-1] = endo_val
            e.put_node_variable_values(varname,timestep+1,nodal_values)

        with Suppressor():
            e.close()

    def set_port_endogenous_values_experimental(self, pc_coeffs, germ, portid):
        """
        Sample polynomial chaos expansion for endogenous values at germ
        Assign those values to the nodes on the exodus mesh
        """

        endo_map = self.get_endogenous_data()

        with Suppressor():
            e = exodus.exodus(self.outfile, mode='a')
        ss_ids = e.get_side_set_ids()
        ss_names = e.get_side_set_names()
        dictionary = dict(zip(ss_names, ss_ids))

        # Get list of nodes
        ssname = self.endogenous_output_ports[portid].ssname
        side_set_node_ids = e.get_side_set_node_list(dictionary[ssname])[1]

        for timestep, node_map in enumerate(endo_map):
            print(e.get_side_set_variable_names())
            #nodal_values = e.get_node_variable_values(self.output_variable,timestep+1)
            nodal_values = e.get_side_set_variable_values(dictionary[ssname], self.output_variable,timestep+1)

            side_set_unique_node_ids = set(side_set_node_ids)
            for nindix, nid in enumerate(side_set_unique_node_ids):
                index = timestep*self.num_endogenous_nodes + node_map.keys().index(nid)
                endo_val = self.evaluate_pce(pc_coeffs[index,...],germ)
                nodal_values[nindex] = endo_val
            #e.put_node_variable_values(self.output_variable,timestep+1,nodal_values)
            e.put_side_set_variable_values(dictionary[ssname], self.output_variable, timestep+1, nodal_values)

        with Suppressor():
            e.close()



    def get_quadpts(self,ndim):
        """
        Generates quadrature points
        Input:
            pc_model: PC object with info about PCE
            ndim: number of dimensions of the PCE
        Output:
            qdpts: numpy array of quadrature points
        """
        # Get the quadrature points
        qdpts_uqtk = uqtkarray.dblArray2D()
        self.pc_model.GetQuadPoints(qdpts_uqtk)
        totquat = self.pc_model.GetNQuadPoints() # Total number of quadrature points

        # Convert quad points to a numpy array
        qdpts = np.zeros((totquat,ndim))
        qdpts_uqtk.getnpdblArray(qdpts)
        return qdpts, totquat


