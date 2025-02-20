
## Rubelt F et al. 2016, five pairs of twins TCRAB Rep-seq (PRJNA300878,  https://doi.org/10.1038/ncomms11112)

### Step 1: Preprocess reads
```
# Download from SRA
accessions=(SRR2905657  SRR2905659  SRR2905663  SRR2905665  SRR2905669  SRR2905671  SRR2905675  SRR2905677  SRR2905681  SRR2905683  SRR2905687  SRR2905689  SRR2905693  SRR2905695  SRR2905699  SRR2905701  SRR2905705  SRR2905707  SRR2905711  SRR2905713 SRR2905658  SRR2905660  SRR2905664  SRR2905666  SRR2905670  SRR2905672  SRR2905676  SRR2905678  SRR2905682  SRR2905684  SRR2905688  SRR2905690  SRR2905694  SRR2905696  SRR2905700  SRR2905702  SRR2905706  SRR2905708  SRR2905712  SRR2905714)
d=../data/libraries/PRJNA300878
mkdir -p ${d}/raw ${d}/merged ${d}/combined
for accession in "${accessions[@]}"
do
	fasterq-dump --split-3 --threads 10 -O ${d}/raw/${accession} ${accession}
	gzip ${d}/${accession}/${accession}_1.fastq
	gzip ${d}/${accession}/${accession}_2.fastq
done

# merge
conda activate igdiscover22
for accession in "${accessions[@]}"
do
	mkdir -p ${d}/merged/${accession}/
	igdiscover merge --threads 10 ${d}/raw/${accession}/${accession}_1.fastq.gz ${d}/raw/${accession}/${accession}_2.fastq.gz ${d}/merged/${accession}/${accession}.fastq.gz
done

# combine libraries per donor
twin_sets=("SRR2905657 SRR2905658 SRR2905659 SRR2905660" "SRR2905663 SRR2905664 SRR2905665 SRR2905666" "SRR2905669 SRR2905670 SRR2905671 SRR2905672" "SRR2905675 SRR2905676 SRR2905677 SRR2905678" "SRR2905681 SRR2905682 SRR2905683 SRR2905684" "SRR2905687 SRR2905688 SRR2905689 SRR2905690" "SRR2905693 SRR2905694 SRR2905695 SRR2905696" "SRR2905699 SRR2905700 SRR2905701 SRR2905702" "SRR2905705 SRR2905706 SRR2905707 SRR2905708" "SRR2905711 SRR2905712 SRR2905713 SRR2905714")
twin_set_names=(1A 1B 2A 2B 3A 3B 4A 4B 5A 5B)
for twin_index in "${!twin_sets[@]}"
do
	mkdir -p ${d}/combined/${twin_set_names[$twin_index]}
	for twin_accession in ${twin_sets[$twin_index]}
	do
		ln -fs ../../merged/${twin_accession}/${twin_accession}.fastq.gz ${d}/combined/${twin_set_names[$twin_index]}/${twin_accession}.fastq.gz
	done
	cat ${d}/combined/${twin_set_names[$twin_index]}/SRR*.fastq.gz > ${d}/combined/${twin_set_names[$twin_index]}/${twin_set_names[$twin_index]}.fastq.gz
done


# split by chain (TRA/TRB) using constant regions
mkdir -p ${d}/TRA/ ${d}/TRB/
for twin_set_name in "${twin_set_names[@]}"
do
	echo ${twin_set_name}
	mkdir -p ${d}/TRA/${twin_set_name} ${d}/TRB/${twin_set_name}
	zgrep --no-group-separator -B 1 AACCCTGACCCTG ${d}/combined/${twin_set_name}/${twin_set_name}.fastq.gz | sed '1~2s/^/>/' > ${d}/TRA/${twin_set_name}/${twin_set_name}.fasta
	zgrep --no-group-separator -B 1 AGGACCTGAAAAACGTG ${d}/combined/${twin_set_name}/${twin_set_name}.fastq.gz | sed '1~2s/^/>/' > ${d}/TRB/${twin_set_name}/${twin_set_name}.fasta
done
```

### Step 2: init, adjust barcode length in config, and run IgDiscover22
### NOTE: The KITDB database used here is v0.0.1 of the database at https://gkhlab.gitlab.io/tcr/sequences/
```
# init
igdiscover batch init --igdiscover-yaml ./igdiscover_yamls/GKH_TCR/TRA.yaml --database /home/mchernys/references/KITDB/TRA/ ../data/libraries/PRJNA300878/TRA/ ../data/igdiscover22/PRJNA300878/TRA/
igdiscover batch init --igdiscover-yaml ./igdiscover_yamls/GKH_TCR/TRB.yaml --database /home/mchernys/references/KITDB/TRB/ ../data/libraries/PRJNA300878/TRB/ ../data/igdiscover22/PRJNA300878/TRB/

# barcode length adjustment
igdiscover batch config --set barcode_length_3prime 0 --set barcode_length_5prime 10 ../data/igdiscover22/PRJNA300878/TRA/
igdiscover batch config --set barcode_length_3prime 0 --set barcode_length_5prime 10 ../data/igdiscover22/PRJNA300878/TRB/

# run
igdiscover batch run ../data/igdiscover22/PRJNA300878/TRA/
igdiscover batch run ../data/igdiscover22/PRJNA300878/TRB/
```

### Step 3: Run CHMMAIRRa
### NOTE: This script runs CHMMAIRRa on all real datasets in the paper, so comment out the datasets you don't want to run.
```
julia +1.10.5 scripts/run_CHMMAIRRa.jl
```

### Step 4: Collect CHMMAIRRa output
```
igdiscover batch collect --file final/CHMMAIRRa_out.tsv.gz ../data/igdiscover22/PRJNA300878/ ../data/igdiscover22/PRJNA300878/collected/collected_CHMMAIRRa_out.tsv
pigz -p 20 ../data/igdiscover22/PRJNA300878/collected/collected_CHMMAIRRa_out.tsv
```
