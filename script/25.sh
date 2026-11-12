sample=$1
chr=$2
mkdir -p data/$sample
if [ ! -e data/$sample/$sample.chr$chr.vcf.gz ]; then
  bcftools view --threads 10 -s $sample -Oz -o data/$sample/$sample.chr$chr.vcf.gz op/chr$chr.op.vcf.gz && tabix data/$sample/$sample.chr$chr.vcf.gz
fi

bcftools view -s ^$sample --threads 10 -Oz -o data/$sample/loo.chr$chr.vcf.gz ../panel2/chr$chr.filt.vcf.gz && tabix data/$sample/loo.chr$chr.vcf.gz

java -jar beagle/beagle.27Feb25.75f.jar ref=data/$sample/loo.chr$chr.vcf.gz gt=data/$sample/$sample.chr$chr.vcf.gz map=gmap/$chr.t2t.gmap.beagle.txt out=data/$sample/$sample.chr$chr.op.b5 nthreads=10 window=100

tabix -f data/$sample/$sample.chr$chr.op.b5.vcf.gz

rm data/$sample/loo.chr$chr.vcf.gz

