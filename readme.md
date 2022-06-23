## ArachnoFamTox: prediction and classification of Arachnids venom and toxin families 
 description

## Installation 

#### Download ArachnoFamTox 
```
git clone https://github.com/yutakajr/ArachnoFamTox.git ArachnoFamTox
```

#### Change to directory 
```
cd ArachnoFamTox
```

#### Create and activate a python virtual environment
```
python -m venv arachnofamtox_env
source arachnofamtox_env/bin/activate
```

#### Install packages in new environment
```
pip install -r requirements.txt
```

#### Install ArachnoFamTox
```
python setup.py install 
```

## Usage

#### Example
```
ArachnoFamTox -fasta test.pep -out out_test_dir 
```

#### Command line options
```
  -h, --help            show this help message and exit
  -fasta <fasta file>   Specify fasta file
  -path <string>        Specify models path. Default=db
  -eHMM <evalue>        e-value for HMMSCAN. Default=1e-1
  -ePSSM <evalue>       e-value for PSSM. Default=1e-5
  -eBLASTP <evalue>     e-value for BLASTp. Default=1e-5
  -qcovsfilter <float>  Filter BLASTp output for qcovs >= <float>. Default=Off
  -pposfilter <float>   Filter BLASTp output for ppos >= <float>. Default=Off
  -pidentfilter <float> Filter BLASTp output for pident >= <float>. Default=Off
  -evaluefilter <float> Filter BLASTp output for evalue <= <evalue>. Default=1e-10
  -cpu <int>            Specify number of threads. Default=1
  -out <output folder>  Specify directory to output
  --force               Force re-use output directory. Default=Off.
  --tempfiles           Maintain temporary files. Default=Off.
``` 

## Running ArachnoFamTox
### Models 

* The default PSSM and HMM models directory is "db" (default). To use ArachnoFamTox databases, make sure to run commands inside directory ArachnoFamTox.

* To test with other PSSM, HMM and BLASTp models, use flag -path to specify the directory with new models. The new models must have the same names and should be specified with -model_name. For example, files for PSSM, HMM and BLASTp are named MODEL.pssm and MODEL.hmm, inside directory "new_models": 
```
ArachnoFamTox.py -fasta test.pep -out out_test_dir -path new_models -model_name MODEL 
```

### Output

Inside output directory, results files are created:
* classification_results: file with proteins classified by PSSM and HMM merged results;
* toxprot_results: file with proteins identified only by BLASTp search against ToxProtDB; 
* out.pssm (optional): Output from RPS-BLAST search against ArachnoFamTox PSSM database; 
* out.hmmer.domtab.parsed (optional): Output from HMMScan search against ArachnoFamTox HMM database; 
* merged_outputs.tsv (optional): Temporary file with merged outputs from PSSM and HMM searches;
* out.blastp.toxprot (optional): Temporary file with toxprot search result. 


## Dependencies

* [**HMMER 3**](https://hmmer.org)  
  Used for HMM profile prediction.   
  *Eddy SR, Accelerated Profile HMM Searches. PLOS Computational Biology 2011, 10.1371/journal.pcbi.1002195*

* [**BLAST 2.11**](https://blast.ncbi.nlm.nih.gov/Blast.cgi)  
  Used for RPS-BLAST and BLASTp prediction.    
  *Altschul, Stephen F., et al. "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic acids research 25.17 (1997): 3389-3402.*


## Licence

[GPL v3](https://github.com/yutakajr/ArachnoFamTox/LICENSE)

## Authors

* [Fernanda Midori Abukawa](https://orcid.org/0000-0002-9304-7566)
* [Milton Yutaka Nishiyama Jr](https://orcid.org/0000-0002-2410-0562)