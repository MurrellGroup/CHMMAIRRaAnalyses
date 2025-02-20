## Galson et al. 2015 (PRJNA308641, https://doi.org/10.1038/icb.2015.57)

### Step 1: After downloading the data, separate it into an IgM and IgG folder.
### SRAtoolkit's fasterq-dump can download the data. 
```
# you can get the SRA codes for the IgM libraries with the following command
grep IgM PRJNA308641_SraRunInfo.csv | cut -c 1-10
```

### Step 2: Cut VH primers off of the reference database and the libraries.
```
VH_primers=PRJNA308641_VH_primers.fasta
IgM_dataset_path=../data/libraries/PRJNA308641/IgM/
IgG_dataset_path=../data/libraries/PRJNA308641/IgG/

# run cutadapt on the reference database
mkdir ../data/OGRDB_human_IGH_9/trimmed/
cutadapt -j 8 --error-rate 0.2 -g file:${VH_primers} -o ../data/OGRDB_human_IGH_9/trimmed/cutadapt-V.fasta ../data/OGRDB_human_IGH_9/V.fasta

# cutadapt for the actual libraries was run at an error rate of 0.1
for d in ${IgM_dataset_path}/SRR*/ ;
do
    id=$(basename ${d})
    mkdir -p ${IgM_dataset_path}/trimmed/${id}/
    cutadapt -m 10 -j 12 --error-rate 0.1 -G file:${VH_primers} -g file:${VH_primers} -o ${IgM_dataset_path}/trimmed/${id}/cutadapt-${id}_1.fastq.gz -p ${IgM_dataset_path}/trimmed/${id}/cutadapt-${id}_2.fastq.gz ${d}/${id}_1.fastq.gz ${d}/${id}_2.fastq.gz;
done

for d in ${IgG_dataset_path}/SRR*/ ;
do
    id=$(basename ${d})
    mkdir -p ${IgG_dataset_path}/trimmed/${id}/
    cutadapt -m 10 -j 12 --error-rate 0.1 -G file:${VH_primers} -g file:${VH_primers} -o ${IgG_dataset_path}/trimmed/${id}/cutadapt-${id}_1.fastq.gz -p ${IgG_dataset_path}/trimmed/${id}/cutadapt-${id}_2.fastq.gz ${d}/${id}_1.fastq.gz ${d}/${id}_2.fastq.gz;
done

```

### Step 3: Run IgDiscover on IgM libraries to infer VDJ genotypes
```
# igdiscover.yaml : no barcode on 3 or 5 primer, reverse_complement true, aux file, and trimmed database with VH_primers.fasta. 
igdiscover batch run --database ../data/OGRDB_human_IGH_9/trimmed/ --igdiscover-yaml igdiscover_yamls/Galson_2015/PRJNA308641_genotyping_igdiscover.yaml ../data/libraries/PRJNA308641_IgM/trimmed/ ../data/igdiscover22/PRJNA308641/
```

### Step 4: Run assignment on IgM and IgG libraries using the genotypes from the previous step
### This is coded up in julia.
```
julia +1.10.5 PRJNA308641_final_inits.jl
```

### Step 5: Run CHMMAIRRa on the IgM and IgG libraries
### NOTE: This script runs CHMMAIRRa on all real datasets in the paper, so comment out the datasets you don't want to run.
```
julia +1.10.5 scripts/run_CHMMAIRRa.jl
```

### Step 6: Collect CHMMAIRRa output
```
igdiscover batch collect --file final/CHMMAIRRa_out.tsv.gz ../data/igdiscover22/PRJNA308641/ ../data/igdiscover22/PRJNA308641/collected/collected_CHMMAIRRa_out.tsv
pigz -p 20 ../data/igdiscover22/PRJNA308641/collected/collected_CHMMAIRRa_out.tsv
```


