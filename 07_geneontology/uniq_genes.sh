# how many unique genes 

cat top_degs_sseqid.txt | awk -F "." '{print $1}' | sort | uniq | wc -l
# 1619 unique genes, same as UniProt ID Mapper and DAVID 

cat top_degs_sseqid.txt | awk -F "." '{print $1}'| sort | uniq > top_degs_unique.txt 



### get baseline gene set ### i.e. all annotated in the whelk transcriptome 
### PANTHER's max is 100,000

cat swissprot.txt | cut -f 2 | cut -d '|' -f2  | sort | uniq | wc -l 

