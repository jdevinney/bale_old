import subprocess
import os

def test_all(path, implementation_mask):

  apps = ["histo", "ig", "topo", "randperm", "transpose_matrix", "permute_matrix", "triangles", "sssp", "unionfind"]
  for app in apps:

    assert(os.path.exists(os.path.join(path,app)))
    runs = []
    runs.append('--help ')
    if app == 'histo' or app == "ig":
      runs.append("-T 10 -n 13084 ")
      runs.append("-T 1 -n 12 ")
      runs.append("-T 5 -n 24242 ")
      runs.append("-n 1289 ")
    if app == 'topo':
      runs.append("-n 1082 -F -e 0.1 ")
      runs.append("-n 2341 -F -e 0.03")
      runs.append("-n 663 -F -e 0.3 ")
      runs.append("-n 1331 -G -e 0.08 ")
      runs.append("-n 133 -G -e 0.3 ")
      runs.append("-n 3332 -G -e 0.03 ")
    if app == 'transpose' or app == 'permute_matrix':
      runs.append("-n 1303 -e 0.1 ")
      runs.append("-n 233 -e 0.1 ")
      runs.append("-n 21 -e 0.3 ")
      runs.append("-n 2100 -e 0.03 ")
    if app == 'randperm':
      runs.append("-n 13093 ")
      runs.append("-n 133 ")
      runs.append("-n 2382 ")
      runs.append("-n 12 ")
    if app == 'triangles':
      runs.append("-n 3109 -e 0.05 ")
      runs.append("-n 601 -e 0.5 ")
      runs.append("-n 4409 -e 0.02 ")
      runs.append("-n 3109 -G -e 0.05 ")
      runs.append("-n 601 -G -e 0.5 ")
      runs.append("-n 4409 -G -e 0.02 ")
      runs.append('-K "0:3x4x5" ')
      runs.append('-K "1:3x4x5" ')
      runs.append('-K "2:3x4x5" ')
    if app == 'sssp':
      runs.append("-n 719 -F -e 0.05 ")
      runs.append("-n 601 -F -e 0.5 ")
      runs.append("-n 1213 -F -e 0.02 ")
      runs.append("-n 1109 -G -e 0.05 ")
      runs.append("-n 601 -G -e 0.5 ")
      runs.append("-n 1409 -G -e 0.02 ")
    if app == 'unionfind':
      runs.append("-n 719 -F -e 0.05 ")
      runs.append("-n 719 -G -e 0.05 ")
      runs.append("-n 601 -F -e 0.02 ")
      runs.append("-n 601 -G -e 0.02 ")
      runs.append("-n 1213 -F -e 0.01 ")
      runs.append("-n 1213 -G -e 0.01 ")
      
    for run in runs:
      cmd = "{0} {1} -M {2}".format(os.path.join(path, app), run, implementation_mask)
      print(cmd)
      cp = subprocess.run(cmd, shell=True)
      assert(cp.returncode == 0)
