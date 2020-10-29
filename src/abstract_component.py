import abc

class Component():
    """Single component in DDUQ network"""
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self):
        pass

    @abc.abstractmethod
    def execute(self):
	"""Execute component forward simulation including evaluation of appropriate endogenous and exogneous inputs"""
        pass

    @abc.abstractmethod
    def add_endogenous_port(self):
	"""Define new endogenous port between 2 components"""
        pass

    @abc.abstractmethod
    def add_exogenous_port(self):
	"""Define new exogenous input for this component"""
        pass

    @abc.abstractmethod
    def add_QoI(self):
        """Specify global variables from exodus output to be treated as QoIs"""
        pass

    @abc.abstractmethod
    def get_endogenous_data(self):
        pass

    @abc.abstractmethod
    def get_QoI_data(self):
        pass

    @abc.abstractmethod
    def get_required_files(self):
        pass

    @abc.abstractmethod
    def add_required_files(self, files):
        pass

    @abc.abstractmethod
    def get_num_endogenous_ports(self):
        pass

    @abc.abstractmethod
    def get_num_exogenous_ports(self):
        pass

    @abc.abstractmethod
    def get_num_pc_coeffs(self):
        pass

    @abc.abstractmethod
    def get_num_QoI(self):
        pass

    @abc.abstractmethod
    def set_all_endogenous_values(self, pce_coeffs, qdpt):
        pass

    @abc.abstractmethod
    def set_port_endogenous_values(self, pce_coeffs, qdpt, link):
        pass

    @abc.abstractmethod
    def get_quadpts(self, pce_dim):
        pass 

    @abc.abstractmethod
    def get_num_timesteps(self):
        pass

    @abc.abstractmethod
    def get_num_pc_terms(self):
        pass

    @abc.abstractmethod
    def GalerkinProjection(self, endogenous_data):
        pass

    @abc.abstractmethod
    def get_all_data(self, varname, filename):
        pass

    @abc.abstractmethod
    def save_nodal_pce(self, nodal_ck, mesh_filename, output_filename):
        pass

    @abc.abstractmethod
    def generate_pc_model(self, pce_dim):
        pass

    @abc.abstractmethod
    def initialize_PCE(self):
        pass

    @abc.abstractmethod
    def get_exogenous_ports(self):
        pass

    @abc.abstractmethod
    def get_QoIs(self):
        pass

    @abc.abstractmethod
    def get_pc_model(self):
        pass
