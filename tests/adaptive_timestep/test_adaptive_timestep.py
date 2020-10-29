import os
path = os.path.dirname(__file__)
import sys
sys.path.append(path+'/../../src')
import sync_times
from subprocess import check_output

def test_get_combined_timeline():
    path = os.path.dirname(__file__)
    t_comb = sync_times.get_combined_timeline([path+"/TL_500.e", path+"/TL_1000.e"])
    assert len(t_comb)==193

def test_interpolate_to_timeline():
    path = os.path.dirname(__file__)
    t_comb = sync_times.get_combined_timeline([path+"/TL_500.e", path+"/TL_1000.e"])

    outputfilename=path+"/TL_500_fine.e"
    if os.path.isfile(outputfilename): os.remove(outputfilename)
    sync_times.interpolate_to_timeline(path+"/TL_500.e", outputfilename, t_comb)

    goldfilename = path+"/TL_500.gold.e"
    check_output(["exodiff", outputfilename, goldfilename])

    outputfilename=path+"/TL_1000_fine.e"
    if os.path.isfile(outputfilename): os.remove(outputfilename)
    sync_times.interpolate_to_timeline(path+"/TL_1000.e", outputfilename, t_comb)

    goldfilename = path+"/TL_1000.gold.e"
    check_output(["exodiff", outputfilename, goldfilename]) 


test_get_combined_timeline()
test_interpolate_to_timeline()
