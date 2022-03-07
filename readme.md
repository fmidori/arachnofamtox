## ArachnoFamTox: prediction and classification of Arachnids venom and toxin families 

Something here 

## Installation 

#### Download ArachnoFamTox 
```
git clone https://github.com/fmidori/arachnofamtox.git ArachnoFamTox
```

#### Change to directory 
```
cd ArachnoFamTox
```

#### Create and activate a conda environment (optional)
```
conda create --name ArachnoFamTox --file requirements.txt
conda activate ArachnoFamTox
```

#### OR Download packages with pip (optional)
```
pip install -r requirements.txt
```

#### Test if ArachnoFamTox is running
```
ArachnoFamTox.py -h
```

## Usage

#### Example
```
ArachnoFamTox.py -fasta test.pep -out out_test_dir 
```

#### Command line options
```
  -h, --help            show this help message and exit
  -fasta <fasta file>   Specify fasta file
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

## Dependencies

* [**HMMER 3**](https://hmmer.org)  
  Used for HMM profile prediction.   
  *Eddy SR, Accelerated Profile HMM Searches. PLOS Computational Biology 2011, 10.1371/journal.pcbi.1002195*

* [**BLAST 2.11**](https://blast.ncbi.nlm.nih.gov/Blast.cgi)  
  Used for RPS-BLAST and BLASTp prediction.    
  *Altschul, Stephen F., et al. "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic acids research 25.17 (1997): 3389-3402.*


## Licence

[GPL v3](https://github.com/fmidori/arachnofamtox/GNU_license)

## Authors

* [Fernanda Midori Abukawa](https://orcid.org/)
* [Milton Yutaka Nishiyama Jr](https://www.researchgate.net/profile/)