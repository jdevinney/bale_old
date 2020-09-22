import subprocess
import os

def run_command(cmd):
  ret = subprocess.run(cmd,shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
  return(ret.returncode)

def determine_launcher():
  ret = run_command('srun --help')
  if ret == 0: return('srun -n {1} {0}')
  ret = run_command('aprun --help')
  if ret == 0: return('aprun -n {1} {0}')
  ret = run_command('oshrun --help')
  if ret == 0: return('oshrun -n {1} {0}')
  ret = run_command('upcrun --help')
  if ret == 0: return('upcrun -n {1} {0}')
  return("{0} -n {1}")

# parameters to this script are handled in conftest.py
# --path : specify a path to executables
def test_all(path, node_range, implementation_mask):

  apps = []
  apps.append("histo")
  apps.append("ig")
  apps.append("randperm")
  apps.append("transpose_matrix")
  apps.append("permute_matrix")
  apps.append("topo")
  apps.append("triangles")
  #apps.append("sssp")
  apps.append("write_sparse_matrix")

  launcher = determine_launcher()

  if node_range is not None:
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

  print()
  for app in apps:
    runs = []
    runs.append("--help ")
    if app == 'histo' or app == 'ig' or app == 'randperm':
      runs.append("-b 10 -n 2 ")
      runs.append("-b 16 -n 1000 ")
      runs.append("-b 35 -n 813 ")
      runs.append("-b 24 -n 2509 ")
    if app == 'histo' or app == 'ig':
      runs.append("-b 35 -n 2344 -T 10 ")
      runs.append("-b 120 -n 1998 -T 10000 ")
    if app == 'topo' or app == 'transpose_matrix' or app == 'permute_matrix' or app == 'triangles' or app == 'sssp':
      runs.append("-b 120 -n 1000 -F -z 2 ")
      runs.append("-b 120 -n 1042 -G -z 4 ")
      runs.append("-b 31 -n 3042 -F -z 4 ")
      runs.append("-b 31 -n 3042 -F -z 6 ")
      runs.append("-b 140 -n 4442 -F -z 20 ")
      runs.append("-b 64 -n 834 -F -d ")
      runs.append("-b 64 -n 834 -F -w ")
      runs.append("-b 64 -n 834 -F -l ")
      runs.append("-b 64 -n 834 -F -d -w ")
      runs.append("-b 64 -n 834 -F -w -l ")
      runs.append("-b 64 -n 834 -F -d -l -w ")

    if app == 'topo':
      if os.path.exists("../../../example_matrices/toposort_input.mm"):
        runs.append("-f ../../../example_matrices/toposort_input.mm")
      #runs.append("-b 512 -n 1340 -G -e .01 ") # this fails for SOS on SMP
    if app == 'triangles' or app == 'transpose_matrix' or app == 'permute_matrix':
      runs.append("-b 244 -K 0:3x4x5 ")
      runs.append("-b 244 -K 1:3x4x5 ")
      runs.append("-b 244 -K 2:3x4x5 ")
      runs.append("-b 345 -K 0:2x4x7 ")
      runs.append("-b 345 -K 1:2x4x7 ")
      runs.append("-b 345 -K 2:2x4x7 ")    
      if os.path.exists("../../../example_matrices/"):
        runs.append("-b 44 -f ../../../example_matrices/undirected_flat_100.mm")
        runs.append("-b 44 -f ../../../example_matrices/undirected_geometric_100.mm")
        if app != 'triangles':
          runs.append("-b 44 -f ../../../example_matrices/directed_flat_100.mm")
          runs.append("-b 44 -f ../../../example_matrices/directed_geometric_100.mm")

    print()
    print("*************************************************************")
    for pes in node_range:
      if pes == 0: continue
      for run in runs:
        cmd = launcher.format(os.path.join(path,app),pes) +" "+run+" -M "+implementation_mask
        print(launcher.format(app, pes)+" "+run+" -M "+implementation_mask)
        ret = run_command(cmd)
        assert(ret == 0)



