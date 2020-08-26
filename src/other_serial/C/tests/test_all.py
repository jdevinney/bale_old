import subprocess
import os

def test_all(path, implementation_mask):

  apps = ["histo", "ig", "topo", "randperm", "transpose_matrix", "permute_matrix", "triangles", "sssp"]
  for app in apps:

    assert(os.path.exists(os.path.join(path,app)))
    runs = []
    runs.append('-h ')
    if app == 'histo' or app == "ig":
      runs.append("-T 10 -n 13084 ")
      runs.append("-T 1 -n 12 ")
      runs.append("-T 5 -N 24242 ")
      runs.append("-n 1289 ")
    if app == 'topo':
      runs.append("-n 1082 -f 0.1 ")
      runs.append("-n 2341 -f 0.03")
      runs.append("-n 663 -f 0.3 ")
      runs.append("-n 1331 -g 0.08 ")
      runs.append("-n 133 -g 0.3 ")
      runs.append("-n 3332 -g 0.03 ")
    if app == 'transpose' or app == 'permute_matrix':
      runs.append("-n 1303 -m 1303 -e 0.1 ")
      runs.append("-n 233 -m 303 -e 0.1 ")
      runs.append("-n 21 -m 30 -e 0.3 ")
      runs.append("-n 2100 -m 2120 -e 0.03 ")
    if app == 'randperm':
      runs.append("-n 13093 ")
      runs.append("-n 133 ")
      runs.append("-n 2382 ")
      runs.append("-n 12 ")
    if app == 'triangles':
      runs.append("-n 3109 -e 0.05 ")
      runs.append("-n 601 -e 0.5 ")
      runs.append("-n 4409 -e 0.02 ")
      runs.append("-n 3109 -g 0.05 ")
      runs.append("-n 601 -g 0.5 ")
      runs.append("-n 4409 -g 0.02 ")
      runs.append('-k "0: 3 4 5" ')
      runs.append('-k "1: 3 4 5" ')
      runs.append('-k "2: 3 4 5" ')
      runs.append('-k "0: 4 5 9" ')
      runs.append('-k "1: 4 5 9" ')
      runs.append('-k "2: 4 5 9" ')
      runs.append('-k "0: 4 5 9 13" ')
      runs.append('-k "1: 4 5 9 13" ')
      runs.append('-k "2: 4 5 9 13" ')
    if app == 'sssp':
      runs.append("-n 719 -e 0.05 ")
      runs.append("-n 601 -e 0.5 ")
      runs.append("-n 1213 -e 0.02 ")
      runs.append("-n 1109 -g 0.05 ")
      runs.append("-n 601 -g 0.5 ")
      runs.append("-n 1409 -g 0.02 ")
      
    for run in runs:
      cmd = "{0} {1} -M {2}".format(os.path.join(path, app), run, implementation_mask)
      print(cmd)
      cp = subprocess.run(cmd, shell=True)
      assert(cp.returncode == 0)
