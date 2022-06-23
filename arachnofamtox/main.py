import argparse
import sys 
from pathlib import Path
import os 
import pandas as pd 
from .tools import run_hmmscan,run_rpsblast,run_blastp,parse_hmmscan
from .define_final_results import define_final,gen_final_class_and_evalue,gen_score,filter_results,define_general_class,define_target_final
from .util import db_dir,output_dir,cpus,search_tools,setup,warn,remove_temp_files,line_prepender


parser = argparse.ArgumentParser(prog="ArachnoFamTox",description= 'Prediction and Classification of Toxins and Venom Proteins for Arachnid species')

parser.add_argument('-fasta', metavar='<fasta file>' ,help='Specify fasta file')

parser.add_argument("-path", metavar='<string>', default="db",
        help="Specify models directory path. Default=db")

parser.add_argument("-model_name", metavar='<string>', default="ArachnoFamTox",
        help="Specify models name. Default=ArachnoFamTox")

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

def merge_and_sort(tsv1, tsv2):
    # Merge two tsvs and sort by first column 
    # Generates merged output
    warn("stdout","Merging predictions")
    file1 = pd.read_csv(tsv1, sep='\t', names=['0','1','2','3','4','5','6','7','8']) # Header 
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

def run():
    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit(1)
    
    setup(args)

    fasta = str(args.fasta)
    domtab = run_hmmscan(fasta,args) 
    hmmscan = parse_hmmscan(domtab)
    rps = run_rpsblast(fasta,args)
    toxprot = run_blastp(fasta,args)
    merged_results = merge_and_sort(hmmscan,rps)
    if have_results(merged_results) == True:
        final = classify_and_gen_out(merged_results)
        gen_toxprot_output(toxprot,final)
    else: 
        warn("stdout","No toxins found")
        gen_toxprot_output(toxprot,"Not_found")
    
    line_prepender(hmmscan,'record\thits\te-value\tbitscore\taccuracy\tlength\tdescription\tstart\tend\n')
    line_prepender(rps, 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs\tppos\n')

    if args.tempfiles == False:
        remove_temp_files(args)
    warn("stdout", "Finished analysis.")
    return None

if __name__ == "__main__":
    run()

