#!/bin/bash


cd /data/mcgaugheyd/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq

bash ~/git/NGS_db/build_metadata.sh Huf | cut -f4 -d',' | grep -v Trek | awk -v OFS="" '{print "scp trek.nhgri.nih.gov:", $1, " ."}' > ~/git/ipsc_rpe_atac/scp_call.sh
bash ~/git/ipsc_rpe_atac/scp_call.sh

bash ~/git/NGS_db/build_metadata.sh Huf > ~/git/ipsc_rpe_atac/metadata.csv
