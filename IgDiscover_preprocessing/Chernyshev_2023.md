
## Chernyshev et al. 2023, 2 individuals IgM and IgG Rep-seq (doi.org/10.17044/scilifelab.21518142, doi.org/10.1038/s41467-023-37972-1)

### For data access, email datacenter@scilifelab.se

### Step 1: For preprocessing, see [the repository](https://gitlab.com/gkhlab/vaccination_of_sars-cov-2_convalescents)

### Step 2: Run chimera detection
### NOTE: This script runs CHMMAIRRa on all real datasets in the paper, so comment out the datasets you don't want to run.
```
julia +1.10.5 scripts/run_CHMMAIRRa.jl
```

### Step 3: Collect chimeras
```
igdiscover batch collect --file final/CHMMAIRRa_out.tsv.gz ../data/igdiscover22/IML369/ ../data/igdiscover22/IML369/collected/collected_CHMMAIRRa_out.tsv
pigz -p 20 ../data/igdiscover22/IML369/collected/collected_CHMMAIRRa_out.tsv
```
