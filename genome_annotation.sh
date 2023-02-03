### Repeats annotation ###
# Repeatmodeler
cmd="cd /bio/pmgonza1/mnark/output/quickmerge/canuNplatnus; module load RepeatModeler/1.0.11;
BuildDatabase -name genome -engine ncbi qm2.def.fasta;
RepeatModeler -srand 2017 -engine ncbi -pa 8 -database genome"
echo "$cmd" | qsub -N RM_canuPiFSC -l nodes=1:ppn=8,walltime=999:00:00,vmem=40gb -q default


# Repeatmasker
cmd="cd /bio/pmgonza1/mnark/output/quickmerge/canuNplatnus; source activate py2.7; 
RepeatMasker -lib genome-families.fa -gff -xsmall -norna -pa 8 qm2.def.fasta"
echo $cmd | qsub -N rmsk_canuPiFSC -l nodes=1:ppn=8,walltime=30:00:00,vmem=40gb -q default
mv qm2.def.fasta.out qm2.def.rmsk.fasta

### Flybase protein identification ###
wget ftp://ftp.flybase.net/genomes/dwil/dwil_r1.05_FB2016_05/gff/dwil-all-r1.05.gff.gz
wget ftp://ftp.flybase.net/genomes/dwil/dwil_r1.05_FB2016_05/fasta/dwil-all-chromosome-r1.05.fasta.gz

gunzip dwil-all-r1.05.gff.gz
gunzip dwil-all-chromosome-r1.05.fasta.gz
outname=dwil-all-r1.05
agat_sp_keep_longest_isoform.pl --gff dwil-all-r1.05.gff -o ${outname}.longest_isoform.gtf
agat_sp_extract_sequences.pl -p --clean_final_stop --clean_internal_stop -g ${outname}.longest_isoform.gtf \
  -f dwil-all-chromosome-r1.05.fasta.gz -o ${outname}.longest_isoform.agat.faa


# Exonerate
exonerate -q ${outname}.longest_isoform.agat.faa -t qm2.def.rmsk.fasta \
  --model protein2genome --querytype protein --targettype dna --showvulgar no \
  --softmaskquery yes --softmasktarget yes --minintron 20 --maxintron 3000 \
  --showalignment no --showtargetgff yes --showcigar no --geneseed 250 --score 250 \
  --verbose 0 --gff3 yes > exonerate.gff3


### Genome guided trancriptome assembly and merging with exonerated proteins ###
# Trimmommatic
for fastq in /LUSTRE/usuario/pgonzale/atlas/output/raw/RNA-seq/*1P.fq.gz; do
java -Xmx15G -Xms15G -jar software/Trimmomatic-0.32/trimmomatic-0.32.jar \
           PE -threads 8 \
          -basein $fastq \
          -baseout output/trimm/Trimmommatic/${fastq%1P.fq.gz}.fq.gz \
          ILLUMINACLIP: software/Trimmomatic-0.32/adapters/TruSeq2-SE.fa:2:30:10 \
          SLIDINGWINDOW:4:15 MINLEN:20
done

# HISAT2
hisat2-build qm2.def.rmsk.fasta qm2.def.rmsk.hisat2_index
for fastq in /LUSTRE/usuario/pgonzale/atlas/output/trimm/Trimmommatic/*1P.fq.gz; do
prefix=$(basename -s 1P.fq.gz $fastq)
hisat2 --new-summary -x qm2.def.rmsk.hisat2_index \
              -1 $fastq \
              -2 ${fastq%1P.fq.gz}2P.fq.gz \
              -p 8 \
              --summary-file ${prefix}.txt \
              -k 100 \
              | samtools view -bS -F 4 -F 8 -F 256 - > HISAT/mapping/${prefix}.bam
done

# Stringtie
for bam in HISAT/mapping/*.out.bam; do
  prefix=$(basename ${bam%.bam})
  stringtie $bam -o StringTie/individual/${prefix}.str.gtf
done
stringtie --merge StringTie/individual/*.gtf -o StringTie/merged/merged.gtf

### PASA ### 
# filter long transcripts fusing multiple genes 
# only those with ortholog to the same gene and spanning neighbor genes will stay
gffread exonerate.gff3 -C -T -o- | grep "mRNA" -v - > noMRNA.gtf
perl /storage/software/PASApipeline-2.0.2/misc_utilities/gtf_to_gff3_format.pl noMRNA.gtf ../../../data/genome/dwl_pilot.fa > noMRNA.gff3
perl /storage/software/PASApipeline-2.1.0/misc_utilities/pasa_gff3_validator.pl noMRNA.gff3
# Prepare stringTie gtf
perl /storage/software/PASApipeline-2.1.0/misc_utilities/cufflinks_gtf_to_alignment_gff3.pl merged.gtf > stpasa.gff3
perl /storage/software/PASApipeline-2.1.0/misc_utilities/pasa_gff3_validator.pl stpasa.gff3
perl /storage/software/PASApipeline-2.1.0/misc_utilities/cufflinks_gtf_genome_to_cdna_fasta.pl merged.gtf ../../../data/genome/dwl_pilot.fa > stpasa.fa
perl /storage/software/PASApipeline-2.1.0/scripts/Load_Current_Gene_Annotations.dbi -c \
  /storage/software/PASApipeline-2.1.0/pasa_conf/conf.txt -C -R \
  -g qm2.def.rmsk.fasta -P stpasa.gff3
perl /storage/software/PASApipeline-2.1.0/scripts/Launch_PASA_pipeline.pl -c \
  /storage/software/PASApipeline-2.1.0/pasa_conf/annotCompare.config -C -R \
  -g qm2.def.rmsk.fasta --cufflinks_gtf stpasa.gff3



### lncRNAs - FeelNc ###
# https://github.com/tderrien/FEELnc
FEELnc_filter.pl -i /storage/Dwillistoni/output/pilot_tissues/StringTie/merged/merged.gtf -a stpasa.gff3 > candidate_lncRNA.gtf
FEELnc_codpot.pl -i candidate_lncRNA.gtf -a stpasa.gff3 -g qm2.def.rmsk.fasta --mode=shuffle
FEELnc_classifier.pl -i feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -a stpasa.gff3 > candidate_lncRNA_classes.txt

  
#### misc RNAs - Rfam ####
rfam=data/dbs/Rfam.cm
fasta=qm2.def.rmsk.fasta
cmpress $rfam
cmsearch --rfam --cpu 2 --tblout genome.rfam $rfam $fasta
awk 'NR>2{ if($10 == "-"){start=$9+1; end=$8+1} else {start=$8; end=$9}; printf "%stinfernalt%st%dt%dt%7.2ft%st.tRfamID=%sn" ,$1,$3,start,end,$15,$10,$4}' genome.rfam | gzip -c > rfam.gff3.gz

#### tRNAs - tRNAscan ####
tRNAscan-SE --thread 2 -E -m trnascan.stats -o trnascan.tab $fasta
awk '$3 ~ /^[0-9]+$/ { if($3 > $4){strand = "-"; start=$4+1; end=$3+1} else {strand = "+"; start=$3; end=$4}; printf "%sttRNAscan-SEttRNAt%dt%dt%7.2ft%st.taa=%s;anti=%sn" ,$1,start,end,$9,strand,$5,$6}' tRNAs.tab | gzip -c > tRNAs.gtf.gz

#### rRNAs - RNAmmer ####
rnammer -S euk -m tsu,lsu,ssu -gff rRNAs.gff < $fasta
