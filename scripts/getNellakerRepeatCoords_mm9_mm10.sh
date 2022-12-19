#!/usr/bin/env bash

set -u

## This script gets the TEV files from Nellaker et al Genome Biol 2012
## and converts into a table that contains mm9 and mm10 coordinates for each TEV.
## The script needs the program liftOver from UCSC
## and the mm9ToMm10.over.chain.gz downloaded from UCSC and saved in ~/Projects/mouse-piRNA-variation/data/forLiftOver/

cd ~/Projects/mouse-piRNA-variation/ # Project dir
mkdir -p data/NellakerGenomeBiol2012/ # Folder for the data from Nellaker
cd data/NellakerGenomeBiol2012/

wget https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2012-13-6-r45/MediaObjects/13059_2012_2884_MOESM13_ESM.TGZ # Get the sup file from Nellaker
tar xzvf 13059_2012_2884_MOESM13_ESM.TGZ  # Uncompress

mkdir -p parsed # Make a directory for the new files

cd bed_files # Folder with files from nellaker

# Add TEV type to the end of each TEV record
cat ETn.shared.tab | awk -F"\t" '{print $0"\tETn"}' >  ../parsed/ETn.shared.withTEVtype.tab
cat IAP-I.shared.tab | awk -F"\t" '{print $0"\tIAP-I"}' >  ../parsed/IAP-I.shared.withTEVtype.tab
cat IS2.shared.tab | awk -F"\t" '{print $0"\tIS2"}' >  ../parsed/IS2.shared.withTEVtype.tab
cat LINE.shared.tab | awk -F"\t" '{print $0"\tLINE"}' >  ../parsed/LINE.shared.withTEVtype.tab
cat LINE_frag.shared.tab | awk -F"\t" '{print $0"\tLINE_frag"}' >  ../parsed/LINE_frag.shared.withTEVtype.tab
cat MaLR.shared.tab | awk -F"\t" '{print $0"\tMaLR"}' >  ../parsed/MaLR.shared.withTEVtype.tab
cat MuLV.shared.tab | awk -F"\t" '{print $0"\tMuLV"}' >  ../parsed/MuLV.shared.withTEVtype.tab
cat RLTR1B.shared.tab | awk -F"\t" '{print $0"\tRLTR1B"}' >  ../parsed/RLTR1B.shared.withTEVtype.tab
cat RLTR10.shared.tab | awk -F"\t" '{print $0"\tRLTR10"}' >  ../parsed/RLTR10.shared.withTEVtype.tab
cat RLTR45.shared.tab | awk -F"\t" '{print $0"\tRLTR45"}' >  ../parsed/RLTR45.shared.withTEVtype.tab
cat VL30.shared.tab | awk -F"\t" '{print $0"\tVL30"}' >  ../parsed/VL30.shared.withTEVtype.tab
cat SINE.shared.tab | awk -F"\t" '{print $0"\tSINE"}' >  ../parsed/SINE.shared.withTEVtype.tab                    

cd ../parsed 

cat *.shared.withTEVtype.tab > allTEVs.shared.withTEVtype.tab # put all TEVs in one file
cat allTEVs.shared.withTEVtype.tab | sort -k1,1 -k2,2n -k3,3n |  grep -v start | awk -F"\t" '{print $0"\tTEV"(n+1);n++}' > allTEVs.shared.withTEVtype.withTEVID.tab # Add an arbitrary TEVID (needed to join mm9 to mm10 coords)
cat allTEVs.shared.withTEVtype.tab | sort -k1,1 -k2,2n -k3,3n |  grep -v start | awk -F"\t" '{print $1"\t"($2-1)"\t"$3"\tTEV"(n+1)"\t.\t "$4;n++}' > allTEVs.shared.TEVID.bed # Generate a bed file for liftOver
liftOver allTEVs.shared.TEVID.bed ~/Projects/piC-variation/data/public/forLiftOver/mm9ToMm10.over.chain.gz allTEVs.shared.TEVID.mm10.bed allTEVs.shared.TEVID.Unmapped # Convert mm9 to mm10
join --nocheck-order -1 24 -2 4 -a 1 -a 2 -o auto -e NA allTEVs.shared.withTEVtype.withTEVID.tab allTEVs.shared.TEVID.mm10.bed  > allTEVs.shared.TEVID.mm9.mm10.tab # Add mm10 coords to Nellaker file with all the rest of the info on each TEV
head -1 allTEVs.shared.withTEVtype.tab | cut -f  5-22 | tr '\t' ' ' | awk '{print "TEVID chr_mm9 start_mm9 end_mm9 strand_mm9 "$0" TEtype chr_mm10 start_mm10 end_mm10 score strand_mm10"}' > header.txt # Get the header (it contains the strain names in the right order)
cat header.txt allTEVs.shared.TEVID.mm9.mm10.tab > tmp.txt  # Add header to final file
mv tmp.txt ../allTEVs.shared.TEVID.mm9.mm10.tab # Rename the final file

# Cleanup
cd ..
rm bed_files/*
rmdir bed_files
rm parsed/*
rmdir parsed
rm 13059_2012_2884_MOESM13_ESM.TGZ
gzip allTEVs.shared.TEVID.mm9.mm10.tab
