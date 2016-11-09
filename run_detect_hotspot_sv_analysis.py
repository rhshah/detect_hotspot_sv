"""
run_processed-psuedogene_analysis
~~~~~~~~~~~~~~~~

:Description: run_processed-psuedogene_analysis is a wrapper to run the processed-psuedogene on MSKCC data
:author:     Ronak H Shah
:copyright:  (c) 2015-2016 by Ronak H Shah for Memorial Sloan Kettering Cancer Center. All rights reserved.
:license:    Apache License 2.0
:contact:    rons.shah@gmail.com
:deffield    updated: Updated

"""


import argparse
import os
import sys
import time
import glob
import pandas as pd
from subprocess import Popen, PIPE
import shlex
import re
# from threading import Thread, Lock
from gridmap import Job, process_jobs


def main():
    
    parser = argparse.ArgumentParser(
        prog='run_processed-psuedogene_analysis.py',
        description='Run processed-psuedogene_analysis on selected pools/samples using MSK data',
        usage='%(prog)s [options]')
    parser.add_argument(
        "-mif",
        "--metaInformationFile",
        action="store",
        dest="mif",
        required=True,
        metavar='folders.fof',
        help="Full path folders and sample names of files.")
    parser.add_argument(
        "-qc",
        "--qcLocation",
        action="store",
        dest="qcLocation",
        required=True,
        metavar='/some/path/qcLocation',
        help="Full path qc files.")
    parser.add_argument(
        "-P",
        "--python",
        action="store",
        dest="python",
        required=False,
        default='/dmp/resources/prod/tools/system/python/production/bin/python',
        metavar='/somepath/python',
        help="Full path Pyhton executables.")
    parser.add_argument(
        "-ppg",
        "--processpsuedogene",
        action="store",
        dest="ppg",
        required=True,
        metavar='/somepath/process-pusedogene.py',
        help="Full path processed_psuedogene.py executables.")
    parser.add_argument(
        "-ias",
        "--iAnnotateSV",
        action="store",
        dest="ias",
        required=True,
        metavar='/somepath/iAnnotateSV.py',
        help="Full path iAnnotate.py executables.")
    parser.add_argument(
        "-q",
        "--queue",
        action="store",
        dest="queue",
        required=False,
        default='test.q',
        metavar='all.q or clin.q',
        help="Name of the SGE queue")
    parser.add_argument(
        "-qsub",
        "--qsubPath",
        action="store",
        dest="qsub",
        required=False,
        default='/common/sge/bin/lx-amd64/qsub',
        metavar='/somepath/qsub',
        help="Full Path to the qsub executables of SGE.")
    parser.add_argument(
        "-t",
        "--threads",
        action="store",
        dest="threads",
        required=False,
        default='5',
        metavar='5',
        help="Number of Threads to be used to run processed-psuedogene analysis on")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        help="make lots of noise [default]")
    parser.add_argument(
        "-o",
        "--outDir",
        action="store",
        dest="outdir",
        required=True,
        metavar='/somepath/output',
        help="Full Path to the output dir.")

    args = parser.parse_args()
    if(args.verbose):
        print "Starting the Process to Run processed-psuedogene."
    
    metaInfoDF = pd.read_table(args.mif,sep="\t")
    for index,rows in metaInfoDF.iterrows():
        poolName = rows.loc['Run']
        id = rows.loc['LIMS_ID']
        delvcf = SetupRun(poolName, id,args)
        pooldir = args.outdir + "/" + poolName
        if(os.path.isdir(pooldir)):
            if(args.verbose):
                print "Pool Output Dir:", pooldir, "exists!!!"
        else:
            os.mkdir(pooldir)
            os.chmod(pooldir, 0o755)
        sampledir = pooldir + "/" + id
        if(os.path.isdir(sampledir)):
            if(args.verbose):
                print "Sample Output Dir:", sampledir, "exists!!!"
        else:
            os.mkdir(sampledir)
            os.chmod(sampledir, 0o755)
        RunPerPool(delvcf, id, sampledir,index,args)
    if(args.verbose):
        print "Finished the Process to Run processed-psuedogene."


def SetupRun(poolName, id, args):
    """This will setup the run to be analyzed.


    :param str poolName: str of pool to be analyzed
    :param str id: str of sample ID to be analyzed
    :param Namespace args: Namespace of args to get other variables
    :return: str with deletion file
    :rtype: str

    """
    # Get the qc data location
    baseqclocation = args.qcLocation + "/" + poolName + "/*"
    # print bamlocation,baseqclocation
    all_qc_subdirs = getSubDirs(baseqclocation)
    # print all_qc_subdirs
    qclocation = max(all_qc_subdirs, key=os.path.getmtime)
    if(os.path.isdir(qclocation)):
        if(args.verbose):
            print "\tQC Location:", qclocation, "\n"
        sv_dir = qclocation + "/StrVarAnalysis/"
        if(os.path.isdir(sv_dir)):
            if(args.verbose):
                print "\tSV Location:", sv_dir, "\n"
            delFile = glob.glob(os.path.join(sv_dir , id + "*/*del.vcf"))
            delFile = delFile[0]
            dupFile = glob.glob(os.path.join(sv_dir , id + "*/*dup.vcf"))
            dupFile = dupFile[0]
            invFile = glob.glob(os.path.join(sv_dir , id + "*/*inv.vcf"))
            invFile = invFile[0]
            jmpFile = glob.glob(os.path.join(sv_dir , id + "*/*jmp.vcf"))
            jmpFile = jmpFile[0]   
    else:
        if(args.verbose):
            print "\tQC LOCATION", qclocation, " DOES NOT EXISTS!!!, Please Review you qcLocation INPUT\n"
    
    return(delFile)


def RunPerPool(vcfFile,id,sampledir,count,args):
    """This will run the pool to be analyzed.


    :param str vcfFile: str of vcf file name
    :param str id: str of sample id
    :param Namespace args: Namespace of args to get other variables
    :return: None
    :rtype: None

    """
    jobs = []
    if(os.path.isfile(vcfFile)):
        jobId = "run_ppg_" + str(count) + "_" + str(id)
        cmdList = []
        cmd = args.python + " " + args.ppg + " " + vcfFile + " " + sampledir + " -s " + id + " --iAnnotateSV " + args.ias + " --genome hg19"
        # cmd = str(cmd)
        threads = int(args.threads)
        threads = threads + 1
        qsub_cmd = args.qsub + " -q " + args.queue + " -N " + jobId + " -o " + jobId + ".stdout" + " -e " + jobId + ".stderr" + \
            " -V -l h_vmem=6G,virtual_free=6G -pe smp " + str(threads) + " -wd " + sampledir + " -sync y " + " -b y " + cmd
        print "qsub_cmd:", qsub_cmd, "\n"
        cmdList.append(qsub_cmd)
        job = Job(
            RunJob,
            cmdList,
            kwlist=None,
            cleanup=True,
            mem_free="2G",
            name=jobId,
            num_slots=1,
            queue=args.queue)
        jobs.append(job)
    print("sending function jobs to cluster")
    print("")

    job_outputs = process_jobs(
        jobs,
        max_processes=10,
        temp_dir='/dmp/analysis/SCRATCH/',
        white_list=None,
        quiet=False,
        local=False)

    print("results from each job")
    for (i, result) in enumerate(job_outputs):
        print("Job {0}- result: {1}".format(i, result))

    return


def RunJob(cmd):
    """Given a command run the job.


    :param str cmd: str of command to be run on the local machine
    :return: None
    :rtype: None

    """
    args = shlex.split(cmd)
    proc = Popen(args)
    proc.wait()
    retcode = proc.returncode
    if(retcode >= 0):
        print "I have finished running process using SGE"

'''
    #iterate over jobs and put each into the queue in sequence
    for job in jobs:
        print "inserting job into the queue: %s"%(job)
        jobqueue.put(job)
    #start some threads, each one will process one job from the queue#
    for i in range(mp.cpu_count()-1):
        th = Thread(target=processor,args = (i,jobqueue))
        th.setDaemon(True)
        th.start()
    #wait until all jobs are processed before quitting
    jobqueue.join()
    with jobqueue.mutex:
        jobqueue.queue.clear()
'''



def getSubDirs(dirLocation):
    """
    Get all sub directories.


    :param str dirLocation: str of directory location
    :return: list of all sub directories
    :rtype: list

    """
    dirs = []
    for d in glob.glob(dirLocation):
        if os.path.isdir(d):
            dirs.append(d)
    return(dirs)


def processor(i, jobqueue):
    """Operate on a jobqueue.


    :param int i: count of the job
    :param Namespace jobqueue: Namespace for jobqueue
    :return: None
    :rtype: None

    """
    devnull = open('/dev/null', 'w')
    if jobqueue.empty():
        print "the Queue is empty!"
        sys.exit(1)
    try:
        job = jobqueue.get()
        print "I'm operating on job item: %s\n" % i, job
        qsub_args = shlex.split(job)
        p = Popen(qsub_args, stdout=PIPE, stderr=devnull)
        jobqueue.task_done()
    except:
        print "Failed to operate on job\n"

if __name__ == "__main__":
    start_time = time.time()
    # mp.freeze_support()
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))