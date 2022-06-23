from datetime import datetime
import sys 
from pathlib import Path
from distutils.spawn import find_executable
import os 
import re
import argparse
import shutil

def db_dir(args):
    # Search for directory with models
    if str(args.path) == 'db':
        if Path('db').is_dir() == False:
            warn("stderr", "Models for PSSM and HMM not found")
            sys.exit(1)
    else:
        if Path(str(args.path).split('/')[0]).is_dir() == False:
            warn("stderr", "Models for PSSM and HMM not found")
            sys.exit(1)

    return None 

def output_dir(args):
    # Create output directory, if already exists, reuse it 
    if os.path.isdir(args.out):
        if args.force:
            warn("stdout", f"Reusing output directory {args.out}")
            shutil.rmtree(args.out)
            os.mkdir(args.out)
        else:
            warn("stderr",
                f"Your choosen output folder '{args.out}' already exist!"
                + " Please change it using --out option or use --force"
                + " to reuse it. ")
            sys.exit(1)
    else:
        warn("stdout",f"Creating output directory {args.out}")
        os.mkdir(args.out)

def cpus(args):
    # Check for available cpus
    cpus = args.cpu
    available_cpus = os.cpu_count()
    
    if args.cpu == 0:
        cpus = available_cpus
    elif args.cpu > available_cpus:
        warn("stdout",
            f"Option --cpus asked for {args.cpu} cores,"
            + f" but system has only {available_cpus}.")
        cpus = available_cpus
    warn("stdout",f"Using {cpus} cores.")
    return None 

def search_tools():
    #Check for needed tools
    needed_tools = ("hmmscan", "rpsblast","blastp")
    for tool in needed_tools:
        if find_executable(tool) is None:
            warn("stderr",f"Tool Not Found: {tool}")
            sys.exit(1)
    return None 

def setup(args):
    #warn("stdout","Setup started")
    db_dir(args)
    output_dir(args)
    cpus(args)
    search_tools()
    return None 

def warn(tp,text):
    now = datetime.now().strftime("%H:%M:%S")
    line = f"[{now}] {text}"
    if tp == "stderr":
        print("[ERROR]" + line, file=sys.stderr)
    elif tp == "stdout":
        print(line, file=sys.stdout)
    return None 

def remove_temp_files(args):
    warn("stdout","Removing temporary files")
    temp_files = ["out.hmmer.domtab", "merged_outputs.tsv","out.blastp.toxprot","out.pssm","out.hmmer.domtab.parsed"]
    for f in temp_files:
        os.remove(Path(args.out, f))
    return None

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line + content)