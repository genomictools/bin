
# HMM
# Create a new copy of hhall.hmm, replacing B1_mean: 100.000000 with 0.000000
# wget -o input/hhall.hmm https://raw.githubusercontent.com/WGLab/PennCNV/master/lib/hhall.hmm
# sed 's/100.000000/0.000000/g' input/hhall.hmm > input/hhall_0.hmm

# ANNOTATION
# wget -o input/cytoBand.txt.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz
# wget -o input/knownGene.txt.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/knownGene.txt.gz
# wget -o input/refGene.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
# wget -o input/refLink.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hgFixed/database/refLink.txt.gz

# GC
# hg38.gc5Base.short.txt

# common_all.vcf.gz