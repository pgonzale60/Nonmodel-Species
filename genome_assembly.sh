# trim the reads from adapter and low quality bases
java -Xmx15G -Xms15G -jar software/Trimmomatic-0.32/trimmomatic-0.32.jar \
           PE -threads 8 \
          -basein illumina_dna_sample_R1.fq.gz \
          -baseout Dwil_illmn_trimm.fq.gz \
          ILLUMINACLIP: software/Trimmomatic-0.32/adapters/TruSeq2-SE.fa:2:30:10 \
          SLIDINGWINDOW:4:15 MINLEN:20
### Platanus assembly ###
module load mchakrab/platanus/1.2.1
platanus assemble -o platanus_Dwil -f \
  Dwil_illmn_trimm* -t 12 -m 300 \
  >& /pub/pmgonza1/Dwil/platanusAssemRun1.log
  
### Explore contamination with blobtools ###

cmd="/software/diamond-0.9.17/diamond blastx --max-target-seqs 25 --verbose \
  -c 1 --threads 20 -q platanus_Dwil_contig.fa \
  -d dbs/diamond/uniprot_ref_proteomes -a \
  platanus_Dwil_contig_uniref;
source activate blobtools_env;
/software/diamond-0.9.17/diamond view --daa platanus_Dwil_contig_uniref.daa > platanus_Dwil_contig_uniref.diamond
blobtools taxify --map_col_sseqid 0 --map_col_taxid 2 -f platanus_Dwil_contig_uniref.diamond \
  -m dbs/diamond/uniprot_ref_proteomes.taxids -o blobtools/"
echo $cmd | qsub -N diamond -l nodes=1:ppn=20,walltime=99:00:00,vmem=50gb -q default

blobtools create -i platanus_Dwil_contig.fa -b contigs.sorted.bam \
  -t platanus_Dwil_contig_uniref.diamond.taxified.out -t \
  platanus_Dwil_contig_uniref.blastout -x bestsumorder 
blobtools view -i blobDB.json --hits --rank all -x bestsumorder
blobtools blobplot -i blobDB.json -x bestsumorder

### Canu assembly ###
# First correct the PacBio reads
cmd="module load canu;
canu -correct -d canu -p Dwil minReadLength=1000 minOverlapLength=500 rawErrorRate=0.3 \
  correctedErrorRate=0.045 genomeSize=250m useGrid=false -pacbio-raw data/pacbio_raw/*.fastq.gz"

echo "$cmd" | qsub -N correctReads -l nodes=1:ppn=30,walltime=9999:00:00,vmem=200gb -q ensam

# Then assemble these corrected reads
cmd="module load canu;
canu -assemble -p Dwil -d output/canu/trim/run1 genomeSize=250m -pacbio-corrected \
  canu/correct/Dwil.correctedReads.fasta.gz minOverlapLength=500 rawErrorRate=0.3 \
  correctedErrorRate=0.045 genomeSize=250m useGrid=false"

echo "$cmd" | qsub -N canuAssem -l nodes=1:ppn=30,walltime=9999:00:00,vmem=230gb -q ensam


### Quickmerge ###
# from https://github.com/mahulchak/quickmerge
module load quickmerge/2017.08.3 MUMmer/4.0.0b

cd canuNplatnus

merge_wrapper.py --clean_only platanus_Dwil_contig.fa canuDwil.fa
mv self_oneline.fa canuDwil.qmFormat.fa
mv hybrid_oneline.fa platanus_Dwil_contig.qmFormat.fa
nucmer -t 64 -l 100 -p qm2 platanus_Dwil_contig.qmFormat.fa canuDwil.qmFormat.fa
delta-filter -i 95 -r -q qm2.delta > qm2.rq.delta
quickmerge -d qm2.rq.delta -q platanus_Dwil_contig.qmFormat.fa -r Dwil_contig.qmFormat.fa -hco 5.0 -c 1.5 -l n -ml m
mv merged.fasta qm2.def.fasta

### Pilon polishing ###

genome=canuNplatnus/qm2.def.fasta
reads=trimmomatic/interleaved/Dwil_illmn_trimm.fq

iter=1
echo -e "#$ -N bwa$iter\n#$ -q bio,free*,pub*\n#$ -pe one-node-mpi 64\n#$ -V \
  \n\ncd canuNplatnus/\nbwa index $genome\nbwa mem -t 64 $genome -p $reads \
  | samtools view -@ 64 -bS - | samtools sort -@ 64 -o pilon$iter.bam \
  -\nsamtools index -@ 64 pilon$iter.bam" > bwa$iter.sh
qsub bwa$iter.sh

echo -e "#$ -N pilon$iter\n#$ -q bio,free*,pub*\n#$ -pe one-node-mpi 64\n#$ -V \
  \n\ncd canuNplatnus/\npilon --genome $genome --frags pilon$iter.bam \
  --output pilon$iter --threads 64 --changes" > pilon$iter.sh
qsub -hold_jid $previous_job_id pilon$iter.sh
