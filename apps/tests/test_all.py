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
    apps = ["histo", "ig", "topo", "randperm", "transpose_matrix", "permute_matrix"]
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

    runs = []
    runs.append("-b 16 -n 1000 -M "+implementation_mask)    
    runs.append("-b 35 -n 813 -T 103 -M "+implementation_mask)
    for app in apps:        
        for pes in node_range:
            if pes == 0: continue
            for run in runs:
                if run.count('-T'):
                    continue
                cmd = launcher.format(os.path.join(path,app),pes) +" "+run
                print(cmd)
                cp = subprocess.run(cmd, shell=True)
                assert(cp.returncode == 0)

            
        
