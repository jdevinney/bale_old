
I put a RSA private key in my home directory at CCS.  The key has a password.  It is in the pwd.txt file.  You should have read permissions for the directory and the files:


```
[npmolin@uccsportal:~/for_jgdevin]$ pwd
/u01/npmolin/for_jgdevin
[npmolin@uccsportal:~/for_jgdevin]$ ll -d
drwxr-x--- 2 npmolin ccs 4096 Apr 23 19:23 .
[npmolin@uccsportal:~/for_jgdevin]$ ll *
-rw-r----- 1 npmolin ccs 1766 Apr 23 18:28 id_rsa_bale_temp
-rw-r----- 1 npmolin ccs   13 Apr 23 19:23 pwd.txt
```

you have an account on an AWS machine that has docker installed.  You have root on the machine.   you can log in with the key I just sent you.   ssh -i /path/to/your/id_rsa_bale_temp jgdevin@35.171.244.97
you login like this:

```
[npmolin@uccsportal:~]$ ssh -i /path/to/your/id_rsa_bale_temp jgdevin@35.171.244.97
Enter passphrase for key 'id_rsa_bale_temp': 
=============================================================================

       __|  __|_  )

       _|  (     /   Deep Learning AMI (Ubuntu 18.04) Version 27.0

      ___|\___|___|

=============================================================================

Welcome to Ubuntu 18.04.4 LTS (GNU/Linux 4.15.0-1060-aws x86_64v)

Please use one of the following commands to start the required environment with the framework of your choice:
for MXNet(+Keras2) with Python3 (CUDA 10.1 and Intel MKL-DNN) ____________________________________ source activate mxnet_p36
for MXNet(+Keras2) with Python2 (CUDA 10.1 and Intel MKL-DNN) ____________________________________ source activate mxnet_p27
for MXNet(+AWS Neuron) with Python3 ___________________________________________________ source activate aws_neuron_mxnet_p36
for TensorFlow(+Keras2) with Python3 (CUDA 10.0 and Intel MKL-DNN) __________________________ source activate tensorflow_p36
for TensorFlow(+Keras2) with Python2 (CUDA 10.0 and Intel MKL-DNN) __________________________ source activate tensorflow_p27
for Tensorflow(+AWS Neuron) with Python3 _________________________________________ source activate aws_neuron_tensorflow_p36
for TensorFlow 2(+Keras2) with Python3 (CUDA 10.1 and Intel MKL-DNN) _______________________ source activate tensorflow2_p36
for TensorFlow 2(+Keras2) with Python2 (CUDA 10.1 and Intel MKL-DNN) _______________________ source activate tensorflow2_p27
for PyTorch with Python3 (CUDA 10.1 and Intel MKL) _____________________________________________ source activate pytorch_p36
for PyTorch with Python2 (CUDA 10.1 and Intel MKL) _____________________________________________ source activate pytorch_p27
for PyTorch (+AWS Neuron) with Python3 ______________________________________________ source activate aws_neuron_pytorch_p36
for Chainer with Python2 (CUDA 10.0 and Intel iDeep) ___________________________________________ source activate chainer_p27
for Chainer with Python3 (CUDA 10.0 and Intel iDeep) ___________________________________________ source activate chainer_p36
for base Python2 (CUDA 10.0) _______________________________________________________________________ source activate python2
for base Python3 (CUDA 10.0) _______________________________________________________________________ source activate python3

Official Conda User Guide: https://docs.conda.io/projects/conda/en/latest/user-guide/
AWS Deep Learning AMI Homepage: https://aws.amazon.com/machine-learning/amis/
Developer Guide and Release Notes: https://docs.aws.amazon.com/dlami/latest/devguide/what-is-dlami.html
Support: https://forums.aws.amazon.com/forum.jspa?forumID=263
For a fully managed experience, check out Amazon SageMaker at https://aws.amazon.com/sagemaker
When using INF1 type instances, please update regularly using the instructions at: https://github.com/aws/aws-neuron-sdk/tree/master/release-notes

=============================================================================

 * Documentation:  https://help.ubuntu.com
 * Management:     https://landscape.canonical.com
 * Support:        https://ubuntu.com/advantage

  System information as of Thu Apr 23 22:34:27 UTC 2020

  System load:  0.06                Processes:              155
  Usage of /:   32.0% of 193.82GB   Users logged in:        1
  Memory usage: 1%                  IP address for ens5:    172.31.94.103
  Swap usage:   0%                  IP address for docker0: 172.17.0.1

 * Ubuntu 20.04 LTS is out, raising the bar on performance, security,
   and optimisation for Intel, AMD, Nvidia, ARM64 and Z15 as well as
   AWS, Azure and Google Cloud.

     https://ubuntu.com/blog/ubuntu-20-04-lts-arrives

35 packages can be updated.
0 updates are security updates.


*** System restart required ***
The programs included with the Ubuntu system are free software;
the exact distribution terms for each program are described in the
individual files in /usr/share/doc/*/copyright.

Ubuntu comes with ABSOLUTELY NO WARRANTY, to the extent permitted by
applicable law.

To run a command as administrator (user "root"), use "sudo <command>".
See "man sudo_root" for details.

jgdevin@ip-172-31-94-103:~$ 
```

You can see you have root/sudo and docker groups

```
jgdevin@ip-172-31-94-103:~$ whoami
jgdevin
jgdevin@ip-172-31-94-103:~$ groups
jgdevin sudo docker
jgdevin@ip-172-31-94-103:~$ # you should git clone bale_private here
```

For now, I'll just put a temporary one here that you can rm -rf

```
jgdevin@ip-172-31-94-103:~$ mkdir neiltemp
jgdevin@ip-172-31-94-103:~$ cd neiltemp
jgdevin@ip-172-31-94-103:~/neiltemp$ git clone https://github.com/jdevinney/bale_private.git
Cloning into 'bale_private'...
Username for 'https://github.com': npmolino
Password for 'https://npmolino@github.com': 
remote: Enumerating objects: 425, done.
remote: Counting objects: 100% (425/425), done.
remote: Compressing objects: 100% (183/183), done.
remote: Total 425 (delta 239), reused 423 (delta 238), pack-reused 0
Receiving objects: 100% (425/425), 694.65 KiB | 24.81 MiB/s, done.
Resolving deltas: 100% (239/239), done.
jgdevin@ip-172-31-94-103:~/neiltemp$ 
```

Next, let's just get the docker image from the link I sent you and see what it can do

```
jgdevin@ip-172-31-94-103:~/neiltemp$ git clone https://github.com/naughtont3/docker-builds.git
Cloning into 'docker-builds'...
remote: Enumerating objects: 1052, done.
remote: Total 1052 (delta 0), reused 0 (delta 0), pack-reused 1052
Receiving objects: 100% (1052/1052), 33.01 MiB | 56.62 MiB/s, done.
Resolving deltas: 100% (661/661), done.
jgdevin@ip-172-31-94-103:~/neiltemp$ cd docker-builds/gups-oshmem/

```

Great.  Now we build the image.  NOTE: You don't have to do this step.  It is for your reference

```
jgdevin@ip-172-31-94-103:~/neiltemp/docker-builds/gups-oshmem$ docker build -t gups-oshmeme .
Sending build context to Docker daemon  23.55kB
Step 1/16 : FROM ubuntu:14.04
14.04: Pulling from library/ubuntu
2e6e20c8e2e6: Pull complete 
30bb187ac3fc: Pull complete 
b7a5bcc4a58a: Pull complete 
Digest: sha256:ffc76f71dd8be8c9e222d420dc96901a07b61616689a44c7b3ef6a10b7213de4
Status: Downloaded newer image for ubuntu:14.04
 ---> 6e4f1fe62ff1
Step 2/16 : MAINTAINER Thomas Naughton <naughtont@ornl.gov>
 ---> Running in a1bb2d506488
Removing intermediate container a1bb2d506488
 ---> cfa1d25a2223
Step 3/16 : ARG NPROCS=4
 ---> Running in 588d7e0d8e1a
Removing intermediate container 588d7e0d8e1a
 ---> ce58000adef5
Step 4/16 : ARG GITHUB_TOKEN=thetoken

ABSURDLY VERBOSE....

Removing intermediate container 617b7260dc6c
 ---> 670400f61ced
Step 12/16 : ENV PATH=/usr/local/bin:$PATH
 ---> Running in 828308517e4c
Removing intermediate container 828308517e4c
 ---> e12d91f552be
Step 13/16 : ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
 ---> Running in bdb47d00da37
Removing intermediate container bdb47d00da37
 ---> d08866d4b32c
Step 14/16 : RUN cd ${PREFIX}/src/ &&     git clone https://github.com/openshmem-org/gups-shmem.git
 ---> Running in 1c8d4d287004
Cloning into 'gups-shmem'...
Removing intermediate container 1c8d4d287004
 ---> 35de87e56e18
Step 15/16 : RUN cd ${PREFIX}/src/gups-shmem &&     make &&     cp gups /usr/local/bin
 ---> Running in 69acfcce2c0d
oshcc -O2 -I./include   -c -o RandomAccess.o RandomAccess.c
RandomAccess.c: In function 'HPCC_SHMEMRandomAccess':
RandomAccess.c:340:3: warning: implicit declaration of function 'HPCC_Power2NodesSHMEMRandomAccessCheck' [-Wimplicit-function-declaration]
   HPCC_Power2NodesSHMEMRandomAccessCheck(logTableSize, TableSize, LocalTableSize,
   ^
oshcc -O2 -I./include   -c -o SHMEMRandomAccess.o SHMEMRandomAccess.c
oshcc -O2 -I./include   -c -o verification.o verification.c
oshcc -O2  RandomAccess.o SHMEMRandomAccess.o verification.o -o gups  -lm
Removing intermediate container 69acfcce2c0d
 ---> 8db3e3407971
Step 16/16 : CMD ["/run-gups.sh"]
 ---> Running in c68f8d176f29
Removing intermediate container c68f8d176f29
 ---> edcb1f2d6e71
Successfully built edcb1f2d6e71
Successfully tagged gups-oshmeme:latest
jgdevin@ip-172-31-94-103:~/neiltemp/docker-builds/gups-oshmem$ 
```


What can this do?  Run it as they intended.  Note you are root and have a prompt in there

```
jgdevin@ip-172-31-94-103:~/neiltemp/docker-builds/gups-oshmem$ docker run --rm -it gups-oshmeme bash
root@fc6e59b7f60e:/# ls
bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  run-gups.sh  sbin  srv  sys  tmp  usr  var
root@fc6e59b7f60e:/# ./run-gups.sh 
Running on 1 processors (PowerofTwo)
Total Main table size = 2^14 = 16384 words
PE Main table size = 2^14 = 16384 words/PE
Default number of updates (RECOMMENDED) = 65536
Real time used = 0.116008 seconds
0.000564926 Billion(10^9) Updates    per second [GUP/s]
0.000564926 Billion(10^9) Updates/PE per second [GUP/s]
Verification:  Real time used = 0.000389 seconds
Found 0 errors in 16384 locations (passed).
root@fc6e59b7f60e:/# exit
exit
jgdevin@ip-172-31-94-103:~/neiltemp/docker-builds/gups-oshmem$ 
```


I have no idea what that run-gups.sh thing is doing, but maybe that output means something to you?
Now let's get bale_private inside the container.  We'll mount part of our drive in /opt with the -v option

```
jgdevin@ip-172-31-94-103:~/neiltemp$ cd docker-builds/gups-oshmem/
jgdevin@ip-172-31-94-103:~/neiltemp/docker-builds/gups-oshmem$ ls
Dockerfile  GET-DOCKER.md  README.md  TODO.md  dockerhub-build-push.sh  files  patches
jgdevin@ip-172-31-94-103:~/neiltemp/docker-builds/gups-oshmem$ cd ../..
jgdevin@ip-172-31-94-103:~/neiltemp$ ls
bale_private  docker-builds
jgdevin@ip-172-31-94-103:~/neiltemp$ docker run --rm -it -v /home/jgdevin/neiltemp/bale_private:/opt/bale_private gups-oshmeme bash
root@713ddf42d5fb:/# ls /opt/bale_private/
Doxyfile  README     apps              convey   install.sh  mainpage.h      runall.sh  spmat
INSTALL   README.md  clang_upc_run.sh  exstack  libgetput   run_autoreconf  serial_C   uconvey.pdf
root@713ddf42d5fb:/# echo YAY WE HAVE BALE IN THE CONTAINER!
YAY WE HAVE BALE IN THE CONTAINER!
root@713ddf42d5fb:/# # let's try it
root@713ddf42d5fb:/# cd /opt/bale_private/
root@713ddf42d5fb:/opt/bale_private# ./install.sh 

*****************************************************
libgetput
*****************************************************
../../libgetput/configure --prefix=/opt/bale_private/build_unknown --with-upc
./install.sh: line 112: ../../libgetput/configure: No such file or directory
configure of libgetput failed!
root@713ddf42d5fb:/opt/bale_private#  # OOF; not sure what's wrong here, but maybe we can trouble shoot this together?
root@713ddf42d5fb:/opt/bale_private# exit
exit
jgdevin@ip-172-31-94-103:~/neiltemp$ 
```

So, I can't get bale_private to install in the container.  Maybe we can poke at this together tomorrow or next week?
