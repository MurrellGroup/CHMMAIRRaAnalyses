# Chernyshev et al. 3 individuals IGM, with 8 PCR conditions each

### For data access, please submit a request to EGA for the accession EGADXXXXXXXXXXX (not available at time of writing due to EGA service issues)

### Step 1: run IgDiscover for genotyping
### Genotyping is done on the libraries created with our standard PCR protocol (labelled C1)
```
igdiscover init --database ../data/OGRDB_human_IGH_9 --reads1 ../data/libraries/PCR_conditions/D1/D1-C1-IGM/D1-C1-IGM_S2_L001_R1_001.fastq.gz ../data/igdiscover22/PCR_conditions_genotyping/D1-C1-IGM/
igdiscover init --database ../data/OGRDB_human_IGH_9 --reads1 ../data/libraries/PCR_conditions/D2/D2-C1-IGM/D2-C1-IGM_S8_L001_R1_001.fastq.gz ../data/igdiscover22/PCR_conditions_genotyping/D2-C1-IGM/
igdiscover init --database ../data/OGRDB_human_IGH_9 --reads1 ../data/libraries/PCR_conditions/D3/D3-C1-IGM/D3-C1-IGM_S5_L001_R1_001.fastq.gz ../data/igdiscover22/PCR_conditions_genotyping/D3-C1-IGM/
igdiscover batch run ../data/igdiscover22/PCR_conditions_genotyping/
```

### Step 2: run IgDiscover assignment on the genotypes from the previous step
```
igdiscover batch init --database ../data/igdiscover22/PCR_conditions_genotyping/D1-C1-IGM/final/database/ ../data/libraries/PCR_conditions/D1/ ../data/igdiscover22/PCR_conditions/D1/
igdiscover batch init --database ../data/igdiscover22/PCR_conditions_genotyping/D2-C1-IGM/final/database/ ../data/libraries/PCR_conditions/D2/ ../data/igdiscover22/PCR_conditions/D2/
igdiscover batch init --database ../data/igdiscover22/PCR_conditions_genotyping/D3-C1-IGM/final/database/ ../data/libraries/PCR_conditions/D3/ ../data/igdiscover22/PCR_conditions/D3/
igdiscover batch config --set iterations 0 ../data/igdiscover22/PCR_conditions/
igdiscover batch run ../data/igdiscover22/PCR_conditions/
```

### Step 3: Run CHMMAIRRa
### NOTE: This script runs CHMMAIRRa on all real datasets in the paper, so comment out the datasets you don't want to run.
```
julia +1.10.5 scripts/run_CHMMAIRRa.jl
```

### Step 4: Collect CHMMAIRRa output
```
igdiscover batch collect --file final/CHMMAIRRa_out.tsv.gz ../data/igdiscover22/PCR_conditions/ ../data/igdiscover22/PCR_conditions/collected/collected_CHMMAIRRa_out.tsv
pigz -p 20 ../data/igdiscover22/PCR_conditions/collected/collected_CHMMAIRRa_out.tsv
```