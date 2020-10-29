import exodus
import os
from output_suppression import Suppressor
import numpy as np

def get_nodal_variable_values(filename, varname, step=1):
    """
    Extracts nodal field data from exodus file and returns a numpy array of nodal values
    """
    with Suppressor():
        e = exodus.exodus(filename, array_type='numpy', mode='r')
        vals=e.get_node_variable_values(varname, step)
        e.close()
    return vals

def get_nodal_variable_names(filename):
    """
    Returns list of nodal variables present in exodus file
    """
    with Suppressor():
        e = exodus.exodus(filename, array_type='numpy', mode='r')
        names = e.get_node_variable_names()
        e.close()
    return names

def get_num_elems(filename):
    """
    Returns the total number of elements in all blocks
    """
    with Suppressor():
        e=exodus.exodus(filename, array_type='numpy', mode='r')
        val=e.num_elems()
        e.close()
    return val

def save_nodal_fields_transient(meshfilename, outputfilename, fieldnames, fielddata):
    # assert len(fieldnames) == fielddata.shape[1]
    # assert get_num_nodes(meshfilename) == fielddata.shape[2]
    if os.path.isfile(outputfilename): os.remove(outputfilename)
    with Suppressor():
        e = exodus.copy_mesh(meshfilename, outputfilename)
        e.close()
        e = exodus.exodus(outputfilename, mode='a', array_type="numpy")
        exodus.add_variables(e, nodal_vars=fieldnames)
        for i,name in enumerate(fieldnames):
            for ts in range(fielddata.shape[0]):
                e.put_node_variable_values(name,ts+1,fielddata[ts,i,:])
                e.put_time(ts+1, ts)
        e.close()

def save_nodal_fields_from_structured(meshfilename, outputfilename, fieldnames, fielddata):
    if os.path.isfile(outputfilename): os.remove(outputfilename)

    with Suppressor():
        e = exodus.exodus(meshfilename, array_type='numpy', mode='r')
        x,y,z = e.get_coords()
        e.close()

    ux = np.unique(x.round(decimals=6))
    uy = np.unique(y.round(decimals=6))

    dx = ux[1] - ux[0]
    dy = uy[1] - uy[0]

    i = np.rint((x-min(x))/dx)
    j = np.rint((y-min(y))/dy)

    with Suppressor():
        e = exodus.copy_mesh(meshfilename, outputfilename)
        e.close()
        e = exodus.exodus(outputfilename, mode='a', array_type="numpy")
        exodus.add_variables(e, nodal_vars=fieldnames)
        for vidx,name in enumerate(fieldnames):
            for ts in range(fielddata.shape[0]):
                e.put_node_variable_values(name,ts+1,fielddata[ts,vidx,i.astype(int),j.astype(int)])
                e.put_time(ts+1, ts)
        e.close()

def save_nodal_fields(meshfilename, outputfilename, fieldnames, fielddata):
    if os.path.isfile(outputfilename): os.remove(outputfilename)
    with Suppressor():
        e = exodus.copy_mesh(meshfilename, outputfilename)
        e.close()
        e = exodus.exodus(outputfilename, mode='a', array_type="numpy")
        exodus.add_variables(e, nodal_vars=fieldnames)

        for i,name in enumerate(fieldnames):
            e.put_node_variable_values(name,1,fielddata[i])

        e.close()

def normalize_data(filename, outputfilename="normalized-output.e"):
    if os.path.isfile(outputfilename): os.remove(outputfilename)

    # Copy Mesh
    e = exodus.copy_mesh(filename, outputfilename)
    e.close()
    exo_out = exodus.exodus(outputfilename, mode='a', array_type="numpy")
    exo_in = exodus.exodus(filename, mode='r', array_type="numpy")

    # Add Variable Names
    var_names = exo_in.get_node_variable_names()
    gvar_names = exo_in.get_global_variable_names()
    exodus.add_variables(exo_out, nodal_vars=var_names)
    exo_out.set_global_variable_number(len(gvar_names))
    for i,gvar in enumerate(gvar_names):
        exo_out.put_global_variable_name(gvar, i+1)

    # Compute Var Min/Max
    minmax = []
    for var in var_names:
        print(var)
        vmin = float("inf")
        vmax = -vmin
        for step in range(exo_in.num_times()):
            data = exo_in.get_node_variable_values(var, step+1)
            vmin = min(vmin, min(data))
            vmax = max(vmax, max(data))

        print((vmin,vmax))
        minmax.append((vmin,vmax))

    # Add Data
    for step in range(exo_in.num_times()):
        for i,var in enumerate(var_names):
            data = exo_in.get_node_variable_values(var, step+1)
            vmin,vmax = minmax[i]
            exo_out.put_node_variable_values(var,step+1,(data-vmin)/(vmax-vmin))
        exo_out.put_time(step+1, step)

    # Add Global Data

    exo_in.close()
    exo_out.close()            

def append_exodus(filenamelist, outputfilename="joined-output.e", skip_first=0, skip_last=0):

    if os.path.isfile(outputfilename): os.remove(outputfilename)

    with Suppressor():
        # Copy Mesh
        e = exodus.copy_mesh(filenamelist[0], outputfilename)
        e.close()
        e = exodus.exodus(outputfilename, mode='a', array_type="numpy")

        # Add Variable Names
        var_names = []
        gvar_names = []
        for f in filenamelist:
            exo = exodus.exodus(f, mode='r', array_type="numpy")
            var_names.extend(exo.get_node_variable_names())
            gvar_names.extend(exo.get_global_variable_names())
            exo.close()
        var_names = list(set(var_names))
        gvar_names = list(set(gvar_names))
        exodus.add_variables(e, nodal_vars=var_names)
        e.set_global_variable_number(len(gvar_names))
        for i,gvar in enumerate(gvar_names):
            e.put_global_variable_name(gvar, i+1)


    # Add Variable Data
    ts = 1
    for f in filenamelist:
        exo = exodus.exodus(f, mode='r', array_type="numpy")
        for step in range(skip_first, exo.num_times()-skip_last):
            for var in exo.get_node_variable_names():
                e.put_node_variable_values(var, ts, exo.get_node_variable_values(var, step+1))
            if len(gvar_names)>0:
                gvar_vals = []
                for gvar in exo.get_global_variable_names():
                    gvar_vals.append(exo.get_global_variable_values(gvar)[step])
                e.put_all_global_variable_values(ts, gvar_vals)
            e.put_time(ts, ts-1)
            ts += 1
        exo.close()
    e.close()

def append_exodus_ss(filenamelist, outputfilename="joined-output.e", labels=None):
    if os.path.isfile(outputfilename): os.remove(outputfilename)

    with Suppressor():
        # Copy Mesh
        e = exodus.copy_mesh(filenamelist[0], outputfilename)
        e.close()
        e = exodus.exodus(outputfilename, mode='a', array_type="numpy")

        # Add Variable Names
        nvar_names = []
        gvar_names = []
        evar_names = []
        for f in filenamelist:
            exo = exodus.exodus(f, mode='r', array_type="numpy")
            nvar_names.extend(exo.get_node_variable_names())
            gvar_names.extend(exo.get_global_variable_names())
            evar_names.extend(exo.get_element_variable_names())
            exo.close()
        if labels:
            gvar_names.extend(labels.keys())
        nvar_names = list(set(nvar_names))
        gvar_names = list(set(gvar_names))
        evar_names = list(set(evar_names))

        exodus.add_variables(e, nodal_vars=nvar_names, element_vars=evar_names)
        e.set_global_variable_number(len(gvar_names))
        for i,gvar in enumerate(gvar_names):
            e.put_global_variable_name(gvar, i+1)

        gvar_vals = {}
        for gvar in gvar_names:
            gvar_vals[gvar] = []

        ts = 1
        for f in filenamelist:
            exo = exodus.exodus(f, mode='r', array_type="numpy")
            step = exo.num_times()
            e.put_time(ts, ts)
            for var in exo.get_node_variable_names():
                e.put_node_variable_values(var, ts, exo.get_node_variable_values(var, step))
            for evar in exo.get_element_variable_names():
                # TODO: only works for 1 block
                e.put_element_variable_values(1, evar, ts, exo.get_element_variable_values(1,evar,step))
            for gvar in exo.get_global_variable_names():
                val = exo.get_global_variable_value(gvar,step)
                gvar_vals[gvar].append(val)
            if labels:
                for key in labels:
                    val = labels[key][ts-1]
                    gvar_vals[key].append(val)
            ts += 1
            exo.close()
        for ts in range(1,e.num_times()+1):
            vals=[]
            for gvar in e.get_global_variable_names():
                vals.append(gvar_vals[gvar][ts-1])
            e.put_all_global_variable_values(ts,vals)
        e.close()

def isin_sideset(filename, ssname):
    with Suppressor():
        e = exodus.exodus(filename, mode='r')
    ss_ids = e.get_side_set_ids()
    ss_names = e.get_side_set_names()
    dictionary = dict(zip(ss_names, ss_ids))

    vals = np.zeros(e.num_nodes())

    ssid = e.get_side_set_node_list(dictionary[ssname])[1]
    side_set_unique_node_ids = set(ssid)
    for nid in side_set_unique_node_ids:
        vals[nid-1] = 1

    with Suppressor():
        e.close()

    return vals

def get_coords(filename):
    """
    Returns the spatial coordinates of nodes in all blocks
    """
    with Suppressor():
        e = exodus.exodus(filename, array_type='numpy', mode='r')
        x,y,z = e.get_coords()
        e.close()
    return x,y,z

def get_node_id_map(filename):
    """
    Returns mapping between node index and node id from exodus file
    """
    with Suppressor():
        e = exodus.exodus(filename, array_type='numpy', mode='r')
        nid=e.get_node_id_map()
        e.close()
    return nid

def add_global_variable(filename, name, vals):
    """
    Adds global variable and fills with values
    """
    with Suppressor():
        e = exodus.exodus(filename, array_type='numpy', mode='a')
        num_gvars = e.get_global_variable_number()
        e.set_global_variable_number(num_gvars+1)
        e.put_global_variable_name(name,num_gvars+1)

        for i,val in enumerate(vals):
            e.put_global_variable_value(name, i+1, val)
        e.close()
    

def get_num_globals(filename):
    """
    Returns number of global variables in the exodus file
    """
    with Suppressor():
        e = exodus.exodus(filename, array_type='numpy', mode='r')
        n=e.get_global_variable_number()
        e.close()
    return n

def get_all_global_variable_values(filename):
    """
    Returns global variable values for all times
    """
    with Suppressor():
        e = exodus.exodus(filename, array_type='numpy', mode='r')
        if e.get_global_variable_number() > 0:
            global_data = np.zeros((e.num_times(),e.get_global_variable_number()))
            for timestep in range(e.num_times()):
                global_data[timestep,:] = e.get_all_global_variable_values(timestep+1)
            return global_data
        else:
            return None

def copy_all_global_variables(filename_source, filename_dest):
    source = exodus.exodus(filename_source, array_type='numpy', mode='r')
    n = source.get_global_variable_number()
    
    if n>0:
        names = source.get_global_variable_names()

        dest = exodus.exodus(filename_dest, array_type='numpy', mode='a')
        dest.set_global_variable_number(n)
        for i,name in enumerate(source.get_global_variable_names()):
            dest.put_global_variable_name(name,i+1)
        for timestep in range(dest.num_times()):
            dest.put_all_global_variable_values(timestep+1, source.get_all_global_variable_values(timestep+1))

def get_num_nodes(filename):
    """
    Returns the total number of nodes in all blocks
    """
    with Suppressor():
        e = exodus.exodus(filename,array_type='numpy', mode='r')
        n = e.num_nodes()
        e.close()
    return n

def get_times(filename):
    """
    Returns list of times corresponding to time planes in the exodus file
    """
    with Suppressor():
        e = exodus.exodus(filename,array_type='numpy', mode='r')
        t = e.get_times()
        e.close()
    return t

def set_times(filename, times):
    with Suppressor():
        e = exodus.exodus(filename,array_type='numpy', mode='a')

    assert len(times)==e.num_times()

    for i,t in enumerate(times):
        e.put_time(i+1, t)

    with Suppressor():
        e.close()
