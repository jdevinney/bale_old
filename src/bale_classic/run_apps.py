# script to run a suite of bale apps for experimentation or testing
import sys
if sys.version_info[0] < 3 or sys.version_info[1] < 7:
  print("This script requires at least Python version 3.7")
  sys.exit(1)

import subprocess
import os
import argparse
from shutil import which

def run_command(cmd):
  #ret = subprocess.run(cmd,shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
  ret = subprocess.run(cmd,shell=True, capture_output=True)
  return(ret)


def determine_launcher(launcher):
  if launcher is not None:
    return(""+launcher+" -n {0} {1} {2}")
  ret = run_command('srun --help')
  if ret.returncode == 0: return('srun -n {0} {1} {2}')
  ret = run_command('aprun --help')
  if ret.returncode== 0: return('aprun -n {0} {1} {2}')
  ret = run_command('oshrun --help')
  if ret.returncode== 0: return('oshrun -n {0} {1} {2}')
  ret = run_command('upcrun --help')
  if ret.returncode== 0: return('upcrun -n {0} {1} {2}')
  return("{2} -n {0} {1} --")

all_apps = []
all_apps.append("histo")
all_apps.append("ig")
all_apps.append("randperm")
all_apps.append("transpose_matrix")
all_apps.append("permute_matrix")
all_apps.append("topo")
all_apps.append("triangles")
#all_apps.append("sssp")
all_apps.append("write_sparse_matrix")


def append_to_file(master_file, file_to_add, first):
    fin = open(file_to_add, 'r')
    newdata = fin.read()
    fin.close()
    fout = open(master_file, 'a')
    if not first:
        fout.write(",")
    fout.write(newdata)
    fout.close()        


def run_apps(app_dict, node_range, launcher_opts, option_str, impl_mask, json_file):
  
  # add user options
  base_cmd = "{0} ".format(option_str)
  
  # add implementation mask to the base_cmd
  if impl_mask is not None:
    base_cmd = "{0} -M {1}".format(base_cmd, impl_mask)

  # if the user wants the output to go to json, add that option
  # and write the initial opening brace in the json file
  if json_file is not None:
    # write initial [ to master json file
    fout = open(json_file,'w')
    fout.write("[\n")
    fout.close()
    tmp_json_file = "tmp_{0}".format(json_file)
    # add -j option to base_cmd
    base_cmd = "{0} -j {1}".format(base_cmd, tmp_json_file)
    
  for app in app_dict:
    first = True
    for pes in node_range:
      cmd = launcher.format(pes, launcher_opts, app_dict[app]) +" "+base_cmd
      print(cmd)
      ret = run_command(cmd)
      assert(ret.returncode == 0)
      lines = ret.stderr.decode('utf-8')
      with open('run_all.log','a') as wp:
        wp.write(lines)
      if json_file is not None:
        append_to_file(json_file, tmp_json_file, first)
        first=False

  if json_file is not None:
    # write closing brace and close json file
    fout = open(json_file,"a")
    fout.write("]")
    fout.close()



############################################################
if __name__ == '__main__':
  
    
  parser = argparse.ArgumentParser(description="""Script to run any or all of bale apps.""")

  parser.add_argument('-a','--app',  action="append", dest='app_list',  help="Specify an app to run. "
                      " Adding multiple -a options will create a list of apps to run. "
                      " Apps must be in {0}".format(all_apps), default=[])
  parser.add_argument('-j','--json_file', action="store", dest='json_file',
                      help="Pass -j option to all apps. The final json file will have a record"
                      " for each run and will be written to filename which is the argument of this option",
                      default=None)
  parser.add_argument('-L','--launcher', action="store", dest="launcher",
                      help="Specify the launcher to use (for example: srun). "
                      "We do our best to automatically detect the launcher, but "
                      "there are cases where we will fail (for instance if you "
                      "have a UPC build but you have oshrun in your PATH).",default=None)
  parser.add_argument(   '--launcher_opts', action="store", dest='launcher_opts',
                      help="Pass these options onto the launcher, these could be srun options for example."
                      " The -n (num_tasks) option to launchers is handled separately"
                      " (using the --node_range option)."
                      " Do not specify the -n option for the launcher here.", default="")
  parser.add_argument('-M','--impl_mask', action="store", dest='impl_mask',
                      help="Pass \"-M IMPL_MASK\" to all apps", default=None)
  parser.add_argument(     '--node_range', action="store", dest='node_range',
                           help="Specify the node range to run on as a range <start>,<end>,<stride>. "
                           "The job will run each app with x PEs where x is in range(node_range) PEs",
                           default="1,4,1")
  parser.add_argument('-o',"--option_str", action="store", dest='option_str',
                      help="Specify a string to pass to all apps. Must be valid for all apps!", default="")
  parser.add_argument('-P','--path',      action="store", dest='path',
                      help="Specify path to binaries. If no path is specified "
                      " this script checks ./ and also the PATH environment variable.", default="")

  
  args = parser.parse_args()

  launcher = determine_launcher(args.launcher)
  
  app_dict = {}
  if len(args.app_list) == 0:
    for app in all_apps:
      app_dict[app] = None
  else:
    for app in args.app_list:
      if app not in all_apps:
        print("I don't know app {0}".format(app))
        exit(1)
      else:
        app_dict[app] = None


  # verify that all app binaries are in our given PATH
  if args.path != "":
    for app in app_dict:
      if not os.path.exists(os.path.join(args.path, app)):
        print("Can't find {0} in path {1}.".format(app, args.path))
        exit(1)
      else:
        app_dict[app] = os.path.join(args.path,app)

  else:
    # the user didn't supply a path, make sure the apps are in the PATH
    for app in app_dict:
      if os.path.exists(app):
        app_dict[app] = "./{0}".format(app)

      elif which(app) is not None:
        app_dict[app] = which(app)
        
      else:
        print("Can't find {0} in user $PATH.".format(app))
        exit(1)
  print()
  for app in app_dict:
    print("Using path {0} for {1}".format(app_dict[app], app))
  print()
    
  #print(type(args.node_range))
  node_range = args.node_range
  l = [int(i) for i in node_range.split(',')]
  if len(l) == 2:
    node_range = range(l[0],l[1])
  if len(l) == 3:
    node_range = range(l[0],l[1],l[2])
  else:
    print("Error: illegal node_range")
    exit(1)
  
  run_apps(app_dict, node_range, args.launcher_opts, args.option_str, args.impl_mask, args.json_file)
