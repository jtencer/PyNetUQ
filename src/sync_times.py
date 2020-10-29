import exodus
import numpy as np
from output_suppression import Suppressor

def get_combined_timeline(filenames):
    times = []
    for filename in filenames:
        with Suppressor():
            e = exodus.exodus(filename, array_type='numpy', mode='r')
            times.append(e.get_times())
            e.close()

    times = [item for sublist in times for item in sublist]
    return sorted(set(times))

def interpolate_to_timeline(input_filename, output_filename, output_times):
    with Suppressor():
        e = exodus.copy_mesh(input_filename, output_filename)
        e.close()

        exo_in = exodus.exodus(input_filename, mode='r', array_type="numpy")
        exo_out = exodus.exodus(output_filename, mode='a', array_type="numpy")

    input_times = exo_in.get_times()
    #if list(input_times) == output_times:
    #    print("skip")
    #    return

    #print("interpolating %d input times to %d output times..." % (len(input_times), len(output_times)))

    node_varnames = exo_in.get_node_variable_names()
    global_varnames = exo_in.get_global_variable_names()

    exodus.add_variables(exo_out, nodal_vars=node_varnames, global_vars=global_varnames)

    for step,time in enumerate(output_times):
        idx = np.searchsorted(input_times, time)
        if idx == 0:
            left_idx = idx
        else:
            left_idx = idx-1

        if idx == len(input_times):
            right_idx = idx-1
        else:
            right_idx = idx

        for varname in node_varnames:
            left_data = exo_in.get_node_variable_values(varname, left_idx+1)
            if left_idx == right_idx:
                exo_out.put_node_variable_values(varname,step+1,left_data)
            else:
                right_data = exo_in.get_node_variable_values(varname, right_idx+1)
                alpha = (time-input_times[left_idx])/(input_times[right_idx]-input_times[left_idx])
                interp_data = left_data + alpha*(right_data-left_data)
                exo_out.put_node_variable_values(varname,step+1,interp_data)

        for varname in global_varnames:
            left_data = exo_in.get_global_variable_value(varname, left_idx+1)
            if left_idx == right_idx:
                exo_out.put_global_variable_value(varname,step+1,left_data)
            else:
                right_data = exo_in.get_global_variable_value(varname, right_idx+1)
                alpha = (time-input_times[left_idx])/(input_times[right_idx]-input_times[left_idx])
                interp_data = left_data + alpha*(right_data-left_data)
                exo_out.put_global_variable_value(varname,step+1,interp_data)            

        exo_out.put_time(step+1, time)

    with Suppressor():
        exo_in.close()
        exo_out.close()
