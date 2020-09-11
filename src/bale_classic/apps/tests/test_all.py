import subprocess
import os

def determine_launcher():
    ret = subprocess.run('srun --help',shell=True)
    if ret.returncode == 0: return('srun -n {1} {0}')
    ret = subprocess.run('aprun --help', shell=True)
    if ret.returncode == 0: return('aprun -n {1} {0}')
    ret = subprocess.run('oshrun --help', shell=True)
    if ret.returncode == 0: return('oshrun -n {1} {0}')
    ret = subprocess.run('upcrun --help', shell=True)
    if ret.returncode == 0: return('upcrun -n {1} {0}')
    return("{0} -n {1}")

# parameters to this script are handled in conftest.py
# --path : specify a path to executables
def test_all(path, node_range, implementation_mask):
    apps = ["histo", "ig", "topo", "randperm", "transpose_matrix", "permute_matrix", "triangles"]
    #apps = ['triangles']
    launcher = determine_launcher()
    print(launcher)
    if node_range is not None:
        print(type(node_range))
        l = [int(i) for i in node_range.split(',')]
        if len(l) == 2:
            node_range = range(l[0],l[1])
        if len(l) == 3:
            node_range = range(l[0],l[1],l[2])
        else:
            print("Error")
            return()
    else:
        node_range = range(1,4)
    print(node_range)


    for app in apps:
        runs = []
        runs.append("--help ")
        if app == 'histo' or app == 'ig':
            runs.append("-b 16 -n 1000 ")
            runs.append("-b 35 -n 813 ")
            runs.append("-b 35 -n 2344 -T 10 ")
            runs.append("-b 120 -n 19988 -T 10000 ")
        if app == 'topo' or app == 'transpose_matrix' or app == 'permute_matrix' or app == 'triangles':
            runs.append("-b 120 -n 1000 -F -z 2")
            runs.append("-b 120 -n 1042 -G -z 4 ")
            runs.append("-b 31 -n 3042 -F -z 4 ")
            runs.append("-b 31 -n 3042 -F -z 6 ")
            runs.append("-b 140 -n 4442 -F -z 30 ")
        if app == 'randperm':
            runs.append("-b 16 -n 1000  ")
            runs.append("-b 35 -n 813 ")
            runs.append("-b 35 -n 2344 ")
            runs.append("-b 120 -n 19988 ")
        if app == 'triangles':            
            runs.append("-b 244 -K 0:3x4x5 ")
            runs.append("-b 244 -K 1:3x4x5 ")
            runs.append("-b 244 -K 2:3x4x5 ")
            runs.append("-b 345 -K 0:2x4x7 ")
            runs.append("-b 345 -K 1:2x4x7 ")
            runs.append("-b 345 -K 2:2x4x7 ")
            
        for pes in node_range:
            if pes == 0: continue
            for run in runs:
                cmd = launcher.format(os.path.join(path,app),pes) +" "+run+" -M "+implementation_mask
                print(cmd)
                cp = subprocess.run(cmd, shell=True)
                assert(cp.returncode == 0)

            
        
