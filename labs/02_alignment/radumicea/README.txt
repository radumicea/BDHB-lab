1) Create a codespace for the main branch of a fork of bozdogalex/BDHB-lab.
2) Download or reuse downloaded FASTA files with sequences.
3) Extract homologous TP53 locus from NG_017013.2, NC_000017.11 and NC_060941.1
(the 3 accessions in the downloaded FASTA file my_tp53.fa). GPT helped find the
relevant genomic locations from the NCBI TP53 gene page (GeneID 7157):
python task_0.py --in_fasta /workspaces/BDHB-lab/data/work/radumicea/lab01/my_tp53.fa --out_fasta /workspaces/BDHB-lab/data/work/radumicea/lab01/my_tp53_extracted.fa
4) Run Task 1 on the extracted regions: python task_1.py --fasta /workspaces/BDHB-lab/data/work/radumicea/lab01/my_tp53_extracted.fa
5) Run Task 2 on the extracted regions (by default on the first 2 sequences, but you can specify which pairs using the args --i1 and -- i2):
python task_2.py --fasta /workspaces/BDHB-lab/data/work/radumicea/lab01/my_tp53_extracted.fa
6) For Task 3, upload /workspaces/BDHB-lab/data/work/radumicea/lab01/my_tp53_extracted.fa to https://www.ebi.ac.uk/jdispatcher/msa/clustalo.