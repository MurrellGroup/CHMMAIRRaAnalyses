# Corcoran et al. 45 individuals TRABDG (doi.org/10.17044/scilifelab.14579091, doi.org/10.1016/j.immuni.2023.01.026)

### For data access, email datacenter@scilifelab.se with the DOI

### Step 1: Split reads into TRA/B/D/G files. All dependencies are included in the IgDiscover conda environment.

```
python3 ./scripts/tcr_locus_splitter.py
```

### Step 2: Copy the [default IgDiscover yaml](https://gitlab.com/gkhlab/igdiscover22/-/blob/main/src/igdiscover/igdiscover.yaml?ref_type=heads) from modify yaml file for the 4 TCR chains using the following commands.
### NOTE: If you are using IgDiscover v1.0.4, you can use the yaml files included in this repository as they are already modified.

```
# improves CDR3 detection in TCRs
# don't use D gene criteria for loci lacking D genes
# pseudo CDR3 detection doesn't work on TCRs
# can't use expansion indicators on TCRs
# can't use expansion indicators on TCRs

# TRA
igdiscover config --set sequence_type TCR \
--set barcode_length_3prime 21 \
--set aux human_gl.aux \
--set ignore_d true \
--set cdr3_location "[-230, -60]" \
--set preprocessing_filter.v_coverage  90 \
--set germline_filter.len_maxfreq_CDR3 1.0 \
--set germline_filter.seq_maxfreq_CDR3 1.0 \
--set pre_germline_filter.unique_js 4 \
--set germline_filter.unique_js 8 \
--set low_gene '["TRAV23/DV6", "TRAV7", "TRAV18", "TRGV1", "TRAV8-6", "TRGV5P", "TRBV4-3"]' \
--file ./IgDiscover_yamls/Corcoran_2023/TCRA.yaml

# TRB
igdiscover config --set sequence_type TCR \
--set barcode_length_3prime 21 \
--set aux human_gl.aux \
--set ignore_d true \
--set cdr3_location "[-230, -60]" \
--set preprocessing_filter.v_coverage  90 \
--set germline_filter.len_maxfreq_CDR3 1.0 \
--set germline_filter.seq_maxfreq_CDR3 1.0 \
--set germline_filter.unique_js 4 \
--set low_gene '["TRAV23/DV6", "TRAV7", "TRAV18", "TRGV1", "TRAV8-6", "TRGV5P", "TRBV4-3"]' \
--file ./IgDiscover_yamls/Corcoran_2023/TCRB.yaml

# TRD
igdiscover config --set sequence_type TCR \
--set barcode_length_3prime 21 \
--set aux human_gl.aux \
--set ignore_d true \
--set cdr3_location "[-130, -80]" \
--set preprocessing_filter.v_coverage  90 \
--set germline_filter.len_maxfreq_CDR3 1.0 \
--set germline_filter.seq_maxfreq_CDR3 1.0 \
--set pre_germline_filter.unique_js 1 \
--set germline_filter.unique_js 1 \
--set low_gene '["TRAV23/DV6", "TRAV7", "TRAV18", "TRGV1", "TRAV8-6", "TRGV5P", "TRBV4-3"]' \
--file ./IgDiscover_yamls/Corcoran_2023/TCRD.yaml

# TRG
igdiscover config --set sequence_type TCR \
--set barcode_length_3prime 21 \
--set aux human_gl.aux \
--set ignore_d true \
--set cdr3_location "[-130, -80]" \
--set preprocessing_filter.v_coverage  90 \
--set germline_filter.len_maxfreq_CDR3 1.0 \
--set germline_filter.seq_maxfreq_CDR3 1.0 \
--set pre_germline_filter.unique_js 1 \
--set germline_filter.unique_js 1 \
--set low_gene '["TRAV23/DV6", "TRAV7", "TRAV18", "TRGV1", "TRAV8-6", "TRGV5P", "TRBV4-3"]' \
--file ./IgDiscover_yamls/Corcoran_2023/TCRG.yaml
```

### NOTE: We added IMGT pseudogene TRGVs to our database following our observation that pseudogene TRGV7*01 was marked as a high frequency chimera
### NOTE: The KITDB database used here is v0.0.1 of the database at https://gkhlab.gitlab.io/tcr/sequences/
### Step 3: init and run IgDiscover
```
nohup igdiscover batch igdiscover --num-pooled 2 --igdiscover-yaml ./IgDiscover_yamls/Corcoran_2023/TRA.yaml --database /home/mchernys/references/KITDB/TRA/ ../data/libraries/GKH_TCR/TRA/ ../data/igdiscover22/GKH_TCR/TRA/ &
nohup igdiscover batch igdiscover --num-pooled 4 --igdiscover-yaml ./IgDiscover_yamls/Corcoran_2023/TRB.yaml --database /home/mchernys/references/KITDB/TRB/ ../data/libraries/GKH_TCR/TRB/ ../data/igdiscover22/GKH_TCR/TRB/ &
nohup igdiscover batch igdiscover --num-pooled 4 --igdiscover-yaml ./IgDiscover_yamls/Corcoran_2023/TRD.yaml --database /home/mchernys/references/KITDB/TRD/ ../data/libraries/GKH_TCR/TRD/ ../data/igdiscover22/GKH_TCR/TRD/ &
nohup igdiscover batch igdiscover --num-pooled 4 --igdiscover-yaml ./IgDiscover_yamls/Corcoran_2023/TRG.yaml --database /home/mchernys/references/KITDB/TRG/ ../data/libraries/GKH_TCR/TRG/ ../data/igdiscover22/GKH_TCR/TRG/ &
nohup igdiscover batch igdiscover --num-pooled 4 --igdiscover-yaml ./IgDiscover_yamls/Corcoran_2023/TRD.yaml --database /home/mchernys/references/KITDB/TRA/ ../data/libraries/GKH_TCR/TRAD/ ../data/igdiscover22/GKH_TCR/TRAD/ &
```


### Step 4: Run CHMMAIRRa
### NOTE: This script runs CHMMAIRRa on all real datasets in the paper, so comment out the datasets you don't want to run.
```
julia +1.10.5 scripts/run_CHMMAIRRa.jl
```

### Step 5: Collect CHMMAIRRa output
```
igdiscover batch collect --file final/CHMMAIRRa_out.tsv.gz ../data/igdiscover22/GKH_TCR/ ../data/igdiscover22/GKH_TCR/collected/collected_CHMMAIRRa_out.tsv
pigz -p 20 ../data/igdiscover22/GKH_TCR/collected/collected_CHMMAIRRa_out.tsv
```

