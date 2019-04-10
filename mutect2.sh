#!/bin/sh

##GATK bundle download
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.gz -O /path/to/GATK/bundle/ucsc.hg19.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz -O /path/to/GATK/bundle/dbsnp_138.hg19.vcf.gz
wget https://ndownloader.figshare.com/files/7354246 -O /path/to/GATK/bundle/TP53.sorted.bed
wget https://ndownloader.figshare.com/files/7354213 -O /path/to/GATK/bundle/CosmicCodingMuts.chr.sort.head.vcf

##ANNOVAR database files download
export PATH=$PATH:/path/to/annovar

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
annotate_variation.pl -buildver hg19 -downdb genomicSuperDups humandb/ 
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp138 humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar cosmic70 humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/ 
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20160302 humandb/



#define the reference file path
GATK=/path/to/GATK/GenomeAnalysisTK.jar
REF=/path/to/GATK/bundle/ucsc.hg19.fasta
DBSNP=/path/to/GATK/bundle/dbsnp_138.hg19.vcf
COSMIC=/path/to/GATK/bundle/CosmicCodingMuts.chr.sort.head.vcf
BED=/path/to/GATK/bundle/TP53.sorted.bed

#create a Panel of Normals (PoN) vcf file from the 4 low-grade tumour samples
for pathandfile in /path/to/EGA/normal/*.clean.recal.TP53.bam ; do
	basewithpath=${pathandfile%.clean.recal.TP53.*}
	basenopath=$(basename $basewithpath)

java -jar $GATK \
	-T MuTect2 \
	-R $REF \
	-I:tumor $(echo $basewithpath).clean.recal.TP53.bam \
	--dbsnp $DBSNP \
	--cosmic $COSMIC \
	--artifact_detection_mode \
	-L $BED \
	-o $(echo $basewithpath).clean.recal.TP53.normal.vcf
done


java -jar $GATK \
	-T CombineVariants \
	-R $REF \
	-V /path/to/EGA/normal/GB544-10_S18.clean.recal.TP53.normal.vcf -V /path/to/EGA/normal/GB624-11_S81.clean.recal.TP53.normal.vcf -V /path/to/EGA/normal/GB730-12_S41.clean.recal.TP53.normal.vcf -V /path/to/EGA/normal/GB909-13_S90.clean.recal.TP53.normal.vcf \
	-minN 2 \
    --setKey "null" \
    --filteredAreUncalled \
    --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
    -L $BED \
    -o /path/to/EGA/normal/MuTect2_PON.vcf
    
    
    

#call somatic variants
for pathandfile in /path/to/EGA/tumor/*.clean.recal.TP53.bam ; do
	basewithpath=${pathandfile%.clean.recal.TP53.*}
	basenopath=$(basename $basewithpath)
	
java -jar $GATK \
	-T MuTect2 \
	-R $REF \
	--dbsnp $DBSNP \
	--cosmic $COSMIC \
	-dt NONE \
	--input_file:tumor $(echo $basewithpath).clean.recal.TP53.bam \
	--intervals $BED \
	-PON /path/to/EGA/normal/MuTect2_PON.vcf \
	-o $(echo $basewithpath).clean.recal.TP53.vcf

done




#ANNOVAR annotation
mkdir /path/to/EGA/tumor/ANNOVAR
cp /path/to/EGA/tumor/*.vcf /path/to/EGA/tumor/ANNOVAR

##convert file format
for pathandfile in /path/to/EGA/tumor/ANNOVAR/*.vcf ; do
	filename=${pathandfile%.*} 
	convert2annovar.pl --format vcf4 --includeinfo --withzyg $pathandfile > $(echo $filename).avinput
done

##add sample information
for pathandfile in /path/to/EGA/tumor/ANNOVAR/*.avinput ; do 
	filename=$(basename $pathandfile)
	mainfilename=${filename%.clean.recal.TP53.avinput}
	pathandmain=${pathandfile%.*}
	awk -v var1="$mainfilename" '{OFS = "\t" ; print $0, var1}' $pathandfile > $(echo $pathandmain).avinputs
done

##merge file
cat /path/to/EGA/tumor/ANNOVAR/*.avinputs > /path/to/EGA/tumor/ANNOVAR/TP53.Tonly.avinput
	
##annotate
table_annovar.pl /path/to/EGA/tumor/ANNOVAR/TP53.Tonly.avinput /path/to/annovar/humandb/ -buildver hg19 -out /path/to/EGA/tumor/ANNOVAR/TP53.Tonly -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,snp138,dbnsfp30a,cosmic70,exac03,clinvar_20160302 -operation g,r,r,f,f,f,f,f,f,f --otherinfo
