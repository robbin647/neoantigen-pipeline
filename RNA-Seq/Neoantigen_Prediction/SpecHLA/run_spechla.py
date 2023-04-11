import logging
import os 
import pandas as pd
import subprocess
import argparse
import sys
from subprocess import Popen, PIPE
from datetime import datetime
## The purpose of writing this manager is to run the spechla process 
##     on a set of related samples, to a designated space
##
## This tool must have the following properties:
##     1. Verbose on error. The automation process should 
##   catch exceptions/errors from its supervised child process
##   and print them to the STDOUT

##     2. Resume with ease. The successful/failed runs
##   should be registered. Thus upon rerun, the manager should 
##   pick up exactly from where it broke down last time.
##   In other words, repeating previous successful stages 
##   should be avoided.

##     3. Separate rule and configuration. Some static info
##   like directory, sample names, should be extracted from
##   the process and written in a pre-defined config file.
##   Such separation should also make sure that the process 
##   can be turned to run on another sample just by changing 
##   the configurations in the config file.

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description='SpecHLA automation script.')

parser.add_argument('--sample-list', dest='sample_list', type=str, required=True, default="sample_list.tsv", help='file containing information of samples')
parser.add_argument('--work-dir', dest='work_dir', type=str, required=True, help='working directory to use when running this script')
parser.add_argument('--threads', dest='threads', type=int, default=5, required=False, help='number of threads to use. Default: 5')
parser.add_argument('--spec-hla-root', dest='spec_hla_root', type=str, required=False, default="/data3/yangtianxia/SOFTWARE/spechla", help='root directory of the SpecHLA software')
parser.add_argument('--spec-hla', dest='spec_hla', type=str, required=False, default="script/whole/SpecHLA.sh", help='path to the SpecHLA.sh script')
parser.add_argument('--extract-hla', dest='extract_hla', type=str, required=False, default="script/ExtractHLAread.sh", help='path to the ExtractHLAread.sh script')
args = parser.parse_args()

if not os.path.exists(args.work_dir):
    logger.info(f"Working dir {args.work_dir} does not exist. It will be created. ")
    os.mkdir(os.path.abspath(args.work_dir))
if not os.path.exists(args.spec_hla_root):
    raise FileNotFoundError(f"Error! SpecHLA dir not found: {args.spec_hla_root}")
    exit(0)
if not os.path.exists(os.path.join(args.spec_hla_root, args.spec_hla)):
    raise FileNotFoundError(f"Error! SpecHLA.sh is not found in your SpecHLA directory")
    exit(0)
if not os.path.exists(os.path.join(args.spec_hla_root, args.extract_hla)):
    raise FileNotFoundError(f"Error! ExtractHLAread.sh is not found in your SpecHLA directory")
    exit(0)

CONFIG = {
            "WORK_DIR": args.work_dir,
            "SPEC_HLA_ROOT": args.spec_hla_root,
            "SAMPLE_LIST": args.sample_list,
            "SPEC_HLA": args.spec_hla,
            "EXTRACT_HLA": args.extract_hla,
            "THREADS": args.threads
        }
# CONFIG = {
#             "WORK_DIR": "/data5/yangtianxia/spechla",
#             "SPEC_HLA_ROOT": "/data3/yangtianxia/SOFTWARE/spechla",
#             "SAMPLE_LIST": "/data3/yangtianxia/SOFTWARE/spechla/automation_20230311/sample_list.tsv",
#             "SPEC_HLA": "script/whole/SpecHLA.sh",
#             "EXTRACT_HLA": "script/ExtractHLAread.sh",
#             "THREADS" : 5
#             }

# What are the steps?
# Step 1. check if spechla_env CONDA env is activated (Need This?)
# Step 2. Extract HLA-realted reads
# bash script/ExtractHLAread.sh -s <sample_id> -b <bamfile> -r <refGenome> -o <outdir>
# Step 3. Exon HLA Typing (-u 1 for exon)
# bash script/whole/SpecHLA.sh -n <sample> -1 <sample.fq.1.gz> -2 <sample.fq.2.gz> -u 1 -o outdir/

# write a program to automate the process of running a list of Bash commands and catch the errors 
# and print them to stdout
def automate_manager(bash_commands, output_file):
    with open(output_file, 'w') as f:
        for command in bash_commands:
            try:
                subprocess.check_call(command, shell=True)
            except subprocess.CalledProcessError as e:
                print(e, file=f)

def prepare_command(config):
    """
    Construct a Bash command to be used by ``run_bash``

    Interpolate the template string with variables in ``config``
    """
    command = """
echo `hostname`; \ 
echo ========Start at $(date)=========== && \\
cd {SPEC_HLA_ROOT} && \\
. /data1/apps/software/Anaconda3/2022\.05/etc/profile.d/conda.sh && \\
conda activate "{SPEC_HLA_ROOT}/spechla_env" && \\
conda env list && \\
    """.format(
        SPEC_HLA_ROOT = config['SPEC_HLA_ROOT']
    )
    sample_table = pd.read_table(config['SAMPLE_LIST']).set_index("sample_id", drop=False)
    for index, row in sample_table.iterrows():
        command += """
echo [`date`] Working on {SAMPLE_ID} && \\
mkdir -p {WORK_DIR}/{SAMPLE_ID} && \\
mkdir -p {WORK_DIR}/{SAMPLE_ID}/extract_hla && \\
bash {SPEC_HLA_ROOT}/{EXTRACT_HLA} -s {SAMPLE_ID} \\
     -b {BAM_FILE} \\
     -r hg38 \\
     -o {WORK_DIR}/{SAMPLE_ID}/extract_hla && \\
mkdir -p {WORK_DIR}/{SAMPLE_ID}/spec_hla && \\
bash {SPEC_HLA_ROOT}/{SPEC_HLA} -n {SAMPLE_ID} \\
     -1 {WORK_DIR}/{SAMPLE_ID}/extract_hla/{SAMPLE_ID}_extract_1.fq.gz -2 {WORK_DIR}/{SAMPLE_ID}/extract_hla/{SAMPLE_ID}_extract_2.fq.gz \\
     -o {WORK_DIR}/{SAMPLE_ID}/spec_hla \\
     -u 1 \\
     -j {THREADS} && \\
echo [`date`] Done with {SAMPLE_ID} && \\
        """.format(
            SAMPLE_ID=row["sample_id"],
            BAM_FILE=row["bam"],
            WORK_DIR=config["WORK_DIR"],
            SPEC_HLA_ROOT=config["SPEC_HLA_ROOT"],
            SPEC_HLA=config["SPEC_HLA"],
            EXTRACT_HLA=config["EXTRACT_HLA"],
            THREADS=config["THREADS"]
        )
    command += f"\necho ========Finish at $(date)============="
    return command

def run_bash(cmd: str) -> None:
    '''
    ``cmd`` is the path to a runnable Bash script
    '''
    if os.path.exists(cmd):
        cmd = os.path.join('.', cmd)
        try:
            logger.info("Executing bash script: %s", cmd)
            worker = Popen(["/usr/bin/sh", cmd])
            # worker = Popen(['/usr/bin/sh', cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout_data, stderr_data = worker.communicate()
            stdout_data = stdout_data.decode()
            logger.debug("Script output:\n%s", stdout_data)
        except subprocess.CalledProcessError as e:
            print("Error!\tError code: " + str(e.returncode))
            print(e)
    else:
        raise FileNotFoundError(f"Tried to execute {os.path.abspath(cmd)} but it was not found!")

if __name__ == "__main__":
    # create a directory to store all generated scripts
    # scripts will be marked distinct by timestamps
    try:
        os.mkdir(".rules")
    except FileExistsError:
        pass
    timestamp = datetime.now().isoformat(timespec="seconds")
    bash_script = os.path.abspath(f".rules/run_spechla.{timestamp}.sh")
    with open(bash_script, "wt" ) as fout:
        fout.write(prepare_command(CONFIG))
    if os.path.exists(bash_script):
       run_bash(bash_script)
    else:
       print(bash_script, "not exists")