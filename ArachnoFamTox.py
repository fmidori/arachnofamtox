import argparse
import sys 
from pathlib import Path
from distutils.spawn import find_executable
import os 
import re 
import subprocess
import warnings
import csv
from Bio import SearchIO 
from operator import itemgetter
import pandas as pd 
import operator
import shutil
from datetime import datetime

parser = argparse.ArgumentParser(prog="ArachnoFamTox.py",description= 'Prediction and Classification of Toxins and Venom Proteins for Arachnid species')

parser.add_argument('-fasta', metavar='<fasta file>' ,help='Specify fasta file')

parser.add_argument('-eHMM',metavar='<evalue>',
        type=float,help='e-value for HMMSCAN. Default=1e-1',default=1e-1) 

parser.add_argument('-ePSSM',metavar='<evalue>', 
        type=float,help='e-value for PSSM. Default=1e-5',default=1e-5)

parser.add_argument('-eBLASTP',metavar='<evalue>', 
        type=float,help='e-value for BLASTp. Default=1e-5',default=1e-5)

parser.add_argument('-qcovsfilter',metavar='<float>', 
        type=float,help='Filter BLASTp output for qcovs >= <float>. Default=Off',default=0)

parser.add_argument('-pposfilter',metavar='<float>', 
        type=float,help='Filter BLASTp output for ppos >= <float>. Default=Off',default=0)

parser.add_argument('-pidentfilter',metavar='<float>', 
        type=float,help='Filter BLASTp output for pident >= <float>. Default=Off',default=0) 

parser.add_argument('-evaluefilter',metavar='<float>', 
        type=float,help='Filter BLASTp output for evalue <= <evalue>. Default=1e-10',default=1e-10) 

parser.add_argument('-cpu', metavar='<int>',
        default=1, type=int,help='Specify number of threads. Default=1')

parser.add_argument('-out', metavar='<output folder>',type=Path,default='output',
        help="Specify directory to output")

parser.add_argument("--force",action="store_true",
        help="Force re-use output directory. Default=Off.")

parser.add_argument("--tempfiles",action="store_true",
        help="Maintain temporary files. Default=Off.")


args = parser.parse_args()

def db_dir():
    # Search for directory with models
    if Path('db').is_dir() == False:
        warn("stderr", "Models for PSSM and HMM not found")
        sys.exit(1)
    return None 

def output_dir():
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

def cpus():
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

def setup():
    warn("stdout","Setup started")
    db_dir()
    output_dir()
    cpus()
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

def run_hmmscan(fasta):
    warn("stdout",f"running hmmscan with evalue {args.eHMM}")
    subprocess.run([
                    "hmmscan",
                    "--cpu", str(args.cpu),
                    "-E" , str(args.eHMM),
                    "--domE", str(args.eHMM),
                    "--domtblout" , Path(args.out, "out.hmmer.domtab"),
                    Path("db", "Arachnida.hmm"),
                    fasta],
                    stdout=subprocess.DEVNULL)
    return str(Path(args.out, "out.hmmer.domtab")) 

def run_rpsblast(fasta):
    warn("stdout",f"running RPS-BLAST with evalue {args.ePSSM}")
    subprocess.run([
                    "rpsblast", 
                    "-db", Path("db", "ArachnidaToxinsV3.pssm"), #here
                    "-query", fasta, 
                    "-out",  Path(args.out, "out.pssm"),
                    "-evalue", str(args.ePSSM),
                    "-outfmt", "6 qseqid sseqid evalue bitscore ppos pident qlen qstart qend",
                    "-num_threads",  str(args.cpu)])
    return str(Path(args.out, "out.pssm")) 

def run_blastp(fasta):
    warn("stdout",f"running BLASTp with evalue {args.eBLASTP}")
    subprocess.run([
                    "blastp",
                    "-db", Path("db", "toxprot"),
                    "-query", fasta, 
                    "-out",  Path(args.out, "out.blastp.toxprot"),
                    "-evalue", str(args.eBLASTP),
                    "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs ppos", 
                    "-num_threads",  str(args.cpu)])
    return str(Path(args.out, "out.blastp.toxprot")) 

def parse_hmmscan(domtab):
    output = open(domtab + '.parsed', 'w')
    with open(domtab, 'r') as f:
        for record in SearchIO.parse(f, 'hmmscan3-domtab'):
            hits = record.hits
            hsps = record.hsps
            num_hits = len(hits)
            num_hsps = len(hsps)

            if num_hits > 0: 
                for i in range(0,num_hits):
                    if (hsps[i].evalue <= 0.1) or (hsps[i].evalue == 0):
                        output.write(str(record.id) + '\t' + str(hits[i].id) + '\t' + str(hsps[i].evalue)  
                            + '\t' + str(hsps[i].bitscore) + '\t' + str(hsps[i].acc_avg) +  
                            '\t' +  str(record.seq_len) + '\t' + str(hits[i].description) + '\t' + 
                            str(hsps[i].env_start) + '\t' + str(hsps[i].env_end) + '\n')
            else:
                if (hsps.evalue <= 0.1) or (hsps.evalue == 0):
                    output.write(str(record.id) + '\t' + str(hits.id) + '\t' + str(hsps.evalue) + 
                        '\t' + str(hsps.bitscore) + '\t' + str(hsps.acc_avg) + 
                        '\t' +  str(record.seq_len) + '\t' + str(hits.description) + '\t' + 
                        str(hsps.env_start) + '\t' + str(hsps.env_end) + '\n' )
    f.close()
    return domtab + ".parsed"

def merge_and_sort(tsv1, tsv2):
    # Merge two tsvs and sort by first column 
    # Generates merged output
    warn("stdout","Merging predictions")
    file1 = pd.read_csv(tsv1, sep='\t', names=['0','1','2','3','4','5','6','7','8'])
    file2 = pd.read_csv(tsv2, sep='\t', names=['0','1','2','3','4','5','6','7','8'])
    concatfiles = pd.concat([file1,file2])
    sorted_concatfiles = concatfiles.sort_values(concatfiles.columns[0], ascending =True)
    sorted_concatfiles.to_csv( Path(args.out, "merged_outputs.tsv"), sep="\t", index=False, header=None)
    return str(Path(args.out, "merged_outputs.tsv"))
    
def sp_lines(tsv):
    # Returns first query of table and last line 
    with open(tsv, 'r') as t:
        prev_query = t.readline().split('\t')[0] #first line
        if '//' not in t:
            with open(tsv, 'a') as t:
                    t.write('//') #when implemented delete if 
    t.close()
    last_line = '//' 
    return prev_query,last_line

def pred(dic,maxormin):
    # Returns min or max value of a dictionary
    res = 'results' 
    maxormin = str(maxormin)
    if maxormin == 'max': 
        res = max(dic.items(), key=operator.itemgetter(1))[0]
    elif maxormin == 'min':
        res = min(dic.items(), key=operator.itemgetter(1))[0]        
    return res

def define_final(famname):
    # Define final name of family 
    final = famname
    final_name = {'KTx':'Scorpion_KTx', 'theraphotoxin':'theraphotoxin', 
                  'NaTx':'Scorpion_NaTx', 'Metalloprotease':'Metalloprotease',
                  'Metalloproteinase':'Metalloprotease', 'ctenitoxin':'ctenitoxin',
                  'hexatoxin':'hexatoxin', 'oxotoxin':'oxotoxin', 'zodatoxin':'zodatoxin',
                  'agatoxin':'agatoxin', 'segestritoxin':'segestritoxin', 
                  'filistatoxin':'filistatoxin', 'lycotoxin':'lycotoxin', 
                  'sparatoxin':'sparatoxin' }
    
    for i in final_name.keys():
        if i in famname:
            final = final_name[i]
    return final 

def gen_final_class_and_evalue(dic,method):

    Target_toxin = { "beta-":"Sodium channels", "gamma-":"Non-specific cations HCN channels",
           "delta-":"Voltage-Gated Sodium channels", "kappa-":"Voltage-Gated Potassium channels", 
           "mu-":"Voltage-Gated Sodium channels", "tau-":"TRP channels",
           "omega-":"Calcium channels", "M-":"Membranolytic Activity"}

    Target = {"KTx":"Potassium channels", "NaTx":"Sodium channels", 
              "CAP":"CRISPs, Antigen-5 or Pathogen-related", "Lectin":"Carbohydrate-binding protein",
              "PeptidaseS1":"Serine Protease", "PhospholipaseA2":"Hydrolysis of phospholipids",
              "Cystatin":"Cysteine Protease Inhibitor", "Kunitz-type":"Serine Protease Inhibitor",
              "TickDefensin":"Antimicrobial Activity", "NDBP":"Antimicrobial Activity", 
              "AstacinLikeMetalloprotease":"Astacin-Like Metalloprotease",
              "TickMetalloprotease":"Neprilysin", "Cytolytic":"Cytolytic Activity",
              "Dermonecrotic":"Sphingomyelin Phosphodiesterase D", "GlycosylHydrolase_Hyaluronidase":"Hydrolase"}
    
    target = "-"

    if method == "hmm": #for hmm get min value 
        result_class = pred(dic,"min")
    elif method == "pssm": #for pssm get max value 
        result_class = pred(dic,"max")
    
    if "toxin" in str(result_class):
        for i in Target_toxin.keys():
            if i in str(result_class):
                target = Target_toxin[i]  #define target for neurotoxins
    else:
        for i in Target.keys():
            if i in str(result_class):
                target = Target[i]  #define target for other toxins

    result_class_final = define_final(result_class) #final name of family
    result_evalue_final = dic[result_class] #final e-value 
    return result_class,result_class_final,result_evalue_final,target


def gen_score(result_hmm_final,result_pssm_final):
    # Generate score from PSSM and HMM results
    result_final = 'None'
    score = 0 
    if result_hmm_final == result_pssm_final:
        result_final = result_hmm_final 
        score = 5
    elif result_hmm_final == 'None':
        result_final = result_pssm_final
        score = 3
    elif result_pssm_final == 'None':
        result_final = result_hmm_final
        score = 2 
    else: 
        result_final = result_pssm_final #use PSSM as final  
        score = 4
    return result_final,score

def filter_results(result_final,score):
    write = True
    if result_final == "GlycosylHydrolase_Hyaluronidase" and score == 2:
        write = False
    elif result_final == "alpha-latrotoxin" and score == 2:
        write = False
    elif result_final == "Dermonecrotic" and score == 2:
        write = False
    return write

def define_general_class(result_final):

    VP = ["Metalloprotease", "CAP", "Cystatin", "GlycosylHydrolase_Hyaluronidase", "Kunitz-type",
      "Lectin", "NDBP","PeptidaseS1","PhospholipaseA2", "Prokineticin", "TickDefensin",
      "TickLipocalin"]

    TX = ["alpha-latrotoxin", "Scorpion_KTx", "Scorpion_NaTx", "theraphotoxin", "ctenitoxin",
        "hexatoxin", "oxotoxin", "zodatoxin", "agatoxin", "segestritoxin", "filistatoxin",
        "lycotoxin", "sparatoxin"]
    
    if result_final in VP:
        gen_class = "Venom Protein"
    else:
        gen_class = "Toxin"
    return gen_class

def define_target_final(target_hmm,target_pssm):
    target_final = '-'
    if target_hmm == target_pssm:
        target_final = target_pssm
    elif target_pssm == '-' and target_hmm != '-':
        target_final = target_hmm
    elif target_hmm == '-' and target_pssm != '-':
        target_final = target_pssm
    else:
        target_final = target_pssm
    return target_final

def have_results(tsv):
    # Check if file exists
    if (os.stat(tsv).st_size) == 0:
        return False
    else: return True

def classify_and_gen_out(tsv):
    warn("stdout","Generating Classification Results")
    famHMM = {}      #dictionary of HMM families with evalues 
    famPSSM = {}     #dictionary of PSSM families with bitscores
    evaluePSSM = {}  #dictionary of PSSM families with evalues
    score = 0        #score for family
    found = 0        #init number of toxins found 

    output = open(Path(args.out, "classification_results.tsv"), "w")
    output.write("Query" + '\t' + "General Classification" + '\t' + "Family" + '\t'
                            + "Possible Target/Function"  + '\t' + "Score" 
                            + '\t' + "PSSM classification" + '\t' + "HMM classification" + '\t' + "PSSM e-value"  
                            + '\t' + "HMM e-value" + '\t' + "PSSM bistcore" + '\n') 
    prev_query,last_line = sp_lines(tsv) # first and last lines 
    with open(tsv, 'r') as t:
        for line in t:
            line = line.strip()
            query = line.split('\t')[0]
            if (query != prev_query) or (line == last_line): 
                target_pssm = '-'
                if bool(famHMM) == True: #if HMM results exist, get final HMM result  
                    result_hmm,result_hmm_final,evalue_hmm_final,target_hmm = gen_final_class_and_evalue(famHMM,"hmm")
                else: 
                    result_hmm_final = 'None'
                    evalue_hmm_final = 'None' 
                    target_hmm = '-'

                if bool(famPSSM) == True: #if PSSM results exist, get final PSSM result 
                    result_pssm,result_pssm_final,bitscore_pssm_final,target_pssm = gen_final_class_and_evalue(famPSSM,"pssm")
                    evalue_pssm_final = evaluePSSM[result_pssm]
                else:   
                    result_pssm_final = 'None' 
                    evalue_pssm_final = 'None' 
                    bitscore_pssm_final = 'None'
                    target_pssm = '-'

                result_final,score = gen_score(result_hmm_final,result_pssm_final) #generate final merged result and score
                write = filter_results(result_final,score) #filter results to write 
                general_class = define_general_class(result_final) #define general class for final family
                target_final = define_target_final(target_hmm,target_pssm) #define final target 

                if write == True:
                    found += 1 
                    output.write(str(prev_query) + '\t' + str(general_class) + '\t' + str(result_final)
                            + '\t' + str(target_final) + '\t' + str(score) + '\t' + str(result_pssm_final) 
                            + '\t' + str(result_hmm_final) + '\t' + str(evalue_pssm_final)  
                            + '\t' + str(evalue_hmm_final) + '\t' + str(bitscore_pssm_final) + '\n') 
                famHMM = {} 
                famPSSM = {} 
                evaluePSSM = {}
                prev_query = query 

            if line == last_line:
               break 
            hit   = line.split('\t')[1]
            fam   = hit.split('|')[0]

            if ("|TX" in hit) or ("|VP" in hit): #Save PSSM evalue and bitscore in dictionary 
                pssm = hit.split('|')[0]
                current_evalue_pssm = float(line.split('\t')[2])
                current_bitscore = float(line.split('\t')[3])
                if pssm not in famPSSM:
                    famPSSM[pssm] = current_bitscore
                if pssm not in evaluePSSM:
                    evaluePSSM[pssm] = current_evalue_pssm 
            else: 
                evalue_hmm = float(line.split('\t')[2]) #Save HMM evalue in dictionary 
                if fam in famHMM:
                    famHMM[fam] *= float(evalue_hmm) #multiply evalues if more than one for family  
                else: 
                    famHMM[fam] = float(evalue_hmm) 
    t.close()
    output.close()
    warn("stdout",f"Found {found} toxins classified by PSSM and HMM.")
    return str(Path(args.out, "classification_results.tsv"))

def read_ids_gen_dic(tsv):
    # Generate list of headers found by PSSM and HMM
    l = []
    with open(tsv, 'r') as t:
        for line in t:
            header = line.split('\t')[0].strip()
            if header not in l:
                l.append(header)
    t.close()
    return l 

def gen_toxprot_output(toxprot,pssmhmm):
    # BLASTp search against toxprot, and filter for headers not found by PSSM or HMM
    warn("stdout",f"Filtering toxprot results for evalue <= {args.evaluefilter}")

    if pssmhmm == "Not_found": # No results from PSSM or HMM
        found_headers = {}
    else:
        found_headers = read_ids_gen_dic(pssmhmm)
    output = open(Path(args.out, "toxprot_results.tsv"), "w")
    with open(toxprot, 'r') as tp:
        found = 0 
        prev_header = "prev_header"
        output.write('qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs\tppos\n')
        for line in tp:
            header = line.split('\t')[0].strip()
            pident = float(line.split('\t')[2].strip())
            evalue_tp = float(line.split('\t')[10].strip())
            qcovs  = float(line.split('\t')[12].strip())
            ppos   = float(line.split('\t')[13].strip())
            if header != prev_header:
                if header not in found_headers:
                    toxprot_out = "Yes"
                    prev_header = header
                    if (pident >= args.pidentfilter) and (qcovs >= args.qcovsfilter) and (ppos >= args.pposfilter) and (evalue_tp <= args.evaluefilter):
                        output.write(str(line))
                        found += 1 
    tp.close()
    warn("stdout",f"Found {found} putative toxins by BLASTp vs ToxProt Database.")
    return str(Path(args.out, "toxprot_results.tsv"))

def remove_temp_files():
    warn("stdout","Removing temporary files")
    temp_files = ["out.hmmer.domtab", "merged_outputs.tsv","out.blastp.toxprot"]
    for f in temp_files:
        os.remove(Path(args.out, f))
    return None

def main():
    setup()
    fasta = str(args.fasta)
    domtab = run_hmmscan(fasta) 
    hmmscan = parse_hmmscan(domtab)
    rps = run_rpsblast(fasta)
    toxprot = run_blastp(fasta)
    merged_results = merge_and_sort(hmmscan,rps)
    if have_results(merged_results) == True:
        final = classify_and_gen_out(merged_results)
        gen_toxprot_output(toxprot,final)
    else: 
        warn("stdout","No toxins found")
        gen_toxprot_output(toxprot,"Not_found")

    if args.tempfiles == False:
        remove_temp_files()
    warn("stdout", "Finished analysis.")
    
    #hmmscan = parse_hmmscan(str(Path(args.out, "out.hmmer.domtab")))
    #rps = str(Path(args.out, "out.pssm"))
    #parsed = str(Path(args.out, "out.hmmer.domtab.parsed"))
    #toxprot = str(Path(args.out, "out.toxprot.blastp"))
    #final = str(Path(args.out, "final_results.tsv"))
    #merged_results = str(Path(args.out, "merged_outputs.tsv"))
    #merged_results = merge_and_sort(parsed,rps)
    #classify_and_gen_out(merged_results)
    #gen_toxprot_output(toxprot,final)
    return None

if __name__ == "__main__":
    main()
    #print(have_results(str(Path(args.out, "merged_outputs.tsv"))))
    #print(define_target_final("hkasjd","asdasd"))
    #classify_and_gen_out(str(Path(args.out, "merged_outputs.tsv")))
    #results(str(Path(args.out, "classification_results.tsv")))
    #toxprot = str(Path(args.out, "out.blastp.toxprot"))
    #final = str(Path(args.out, "classification_results.tsv"))
    #gen_toxprot_output(toxprot,final)
