import operator
from operator import itemgetter

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