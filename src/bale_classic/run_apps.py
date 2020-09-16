import subprocess
import os
import sys

def determine_launcher():
    ret = subprocess.run('srun --help',capture_output=True,shell=True)
    if ret.returncode == 0: return('srun -n {1} {0}')
    ret = subprocess.run('aprun --help', capture_output=True,shell=True)
    if ret.returncode == 0: return('aprun -n {1} {0}')
    ret = subprocess.run('oshrun --help', capture_output=True,shell=True)
    if ret.returncode == 0: return('oshrun -n {1} {0}')
    ret = subprocess.run('upcrun --help', capture_output=True,shell=True)
    if ret.returncode == 0: return('upcrun -n {1} {0}')
    return("{0} -n {1}")

apps = ["histo", "ig", "topo", "randperm", "transpose_matrix", "permute_matrix", "triangles", "write_sparse_matrix"]


def append_to_file(master_file, file_to_add, first):
    fin = open(file_to_add, 'r')
    newdata = fin.read()
    fin.close()
    fout = open(master_file, 'a')
    if not first:
        fout.write(",")
    fout.write(newdata)
    fout.close()        
        
def run_app(path, node_range, app_list, option_str, impl_mask, json_file):
  for app in app_list:
    runs = []
    if option_str == None:
      runs.append(" ")
    else:
      runs.append("{0} ".format(option_str))
      
    if impl_mask is not None:
      for i,run in enumerate(runs):
        runs[i] = "{0} -M {1}".format(run, impl_mask)
    first = True
    if json_file is not None:
      fout = open(json_file,'w')
      fout.write("[\n")
      fout.close()
      tmp_json_file = "tmp_{0}".format(json_file)  
      for i,run in enumerate(runs):
        runs[i] = "{0} -j {1}".format(run, tmp_json_file)

    
    for pes in node_range:
      #if pes == 0: continue      
      for run in runs:
        cmd = launcher.format(os.path.join(path,app),2**pes) +" "+run
        print(cmd)
        ret = subprocess.run(cmd, capture_output=True, shell=True)
        assert(ret.returncode == 0)
        lines = ret.stderr.decode('utf-8')
        with open('run_all.log','a') as wp:
          wp.write(lines)
        if json_file is not None:
            append_to_file(json_file, tmp_json_file, first)
            first=False
            
    if json_file is not None:
        fout = open(json_file,"a")
        fout.write("]")
        fout.close()

if __name__ == '__main__':

  if sys.version_info[0] < 3 or sys.version_info[1] < 7:
    print("This script requires at least Python version 3.7")
    sys.exit(1)
  
  from optparse import OptionParser
  parser = OptionParser(description="""Script to run any or all of bale apps.""", usage="python %prog [options]")

  parser.add_option('-a','--app_list',  action="store", dest='app_list',  help="Specify a list of apps to run. Must be of form [<app>, <app2>,...] "
                    " and apps must be in {0}".format(apps), default="All")
  parser.add_option('-j','--json_file', action="store", dest='json_file', help="Pass -j <file> to all apps", default=None)
  parser.add_option('-M','--impl_mask', action="store", dest='impl_mask', help="Pass -M <mask> to all apps", default=None)
  parser.add_option(    '--node_range', action="store", dest='node_range', help="Specify the node range to run on as a range <start>:<end>:<stride>. "
                        "The job will run each app with 2^x PEs where x is in range(node_range) PEs", default="1,4,1")
  parser.add_option('-o',"--option_str", action="store", dest='option_str', help="Specify a string to pass to all apps (must be valid for all apps!", default=None)
  parser.add_option('-P','--path',      action="store", dest='path',      help="Specify path to binaries", default=None)
  
  (options, ignore) = parser.parse_args()

  launcher = determine_launcher()
  #print(launcher)

  if options.app_list == 'All':
    app_list = apps
  else:
    #print(options.app_list)
    app_list = options.app_list
    if not app_list.startswith('[') or not app_list.endswith(']'):
      print("bad app list")
      exit(0)
    app_list = app_list[1:-1]
    app_list = app_list.split(",")
    app_list = [i.strip() for i in app_list]
    #print(app_list)
    
  #print(type(options.node_range))
  node_range = options.node_range
  l = [int(i) for i in node_range.split(',')]
  if len(l) == 2:
    node_range = range(l[0],l[1])
  if len(l) == 3:
    node_range = range(l[0],l[1],l[2])
  else:
    print("Error: illegal node_range")
    exit(1)
  
  
  run_app(options.path, node_range, app_list, options.option_str, options.impl_mask, options.json_file)
