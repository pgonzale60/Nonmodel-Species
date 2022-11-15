### Canu assembly ###
cd /LUSTRE/usuario/pgonzale/DwilAssem
cmd="cd /LUSTRE/usuario/pgonzale/DwilAssem; module load canu;
canu -correct -d output/canu/correct/run1 -p Dwil minReadLength=1000 minOverlapLength=500 rawErrorRate=0.3 correctedErrorRate=0.045 genomeSize=250m useGrid=false -pacbio-raw data/raw/*.fastq.gz"

echo "$cmd" | qsub -N correctReads -l nodes=1:ppn=30,walltime=9999:00:00,vmem=200gb -q ensam

cmd="cd /LUSTRE/usuario/pgonzale/DwilAssem; module load canu;
canu -assemble -p Dwil -d output/canu/trim/run1 genomeSize=250m -pacbio-corrected output/canu/correct/run1/Dwil.correctedReads.fasta.gz minOverlapLength=500 rawErrorRate=0.3 correctedErrorRate=0.045 genomeSize=250m useGrid=false"

echo "$cmd" | qsub -N canuAssem -l nodes=1:ppn=30,walltime=9999:00:00,vmem=230gb -q ensam

### Platanus assembly ###
module load mchakrab/platanus/1.2.1
platanus assemble -o /pub/pmgonza1/Dwil/output/platanus/rundef/Dwil -f /pub/pmgonza1/Dwil/output/trimmomatic/unzip/Dwil_illmn_trimm_* -t 12 -m 300 >& /pub/pmgonza1/Dwil/platanusAssemRun1.log
platanus scaffold -o /pub/pmgonza1/Dwil/output/platanus/rundef/Dwil -c /pub/pmgonza1/Dwil/output/platanus/rundef/Dwil_contig.fa -b /pub/pmgonza1/Dwil/output/platanus/rundef/Dwil_contigBubble.fa -IP1 /pub/pmgonza1/Dwil/output/trimmomatic/unzip/Dwil_illmn_trimm_1P.fq /pub/pmgonza1/Dwil/output/trimmomatic/unzip/Dwil_illmn_trimm_2P.fq  -t 12 2> /pub/pmgonza1/Dwil/platanusScaffRun1.log
platanus gap_close /pub/pmgonza1/Dwil/output/platanus/rundef/Dwil -c /pub/pmgonza1/Dwil/output/platanus/rundef/Dwil_scaffold.fa -IP1 /pub/pmgonza1/Dwil/output/trimmomatic/unzip/Dwil_illmn_trimm_1P.fq /pub/pmgonza1/Dwil/output/trimmomatic/unzip/Dwil_illmn_trimm_2P.fq  -t 12 2> /pub/pmgonza1/Dwil/platanusGapRun1.log

### Explore contamination with blobtools ###

cmd="/LUSTRE/usuario/pgonzale/software/diamond-0.9.17/diamond blastx --max-target-seqs 25 --verbose \
-c 1 --threads 20 -q /LUSTRE/usuario/pgonzale/DwilAssem/output/platanus/Dwil_contig.fa -d /LUSTRE/usuario/pgonzale/dbs/diamond/uniprot_ref_proteomes -a /LUSTRE/usuario/pgonzale/DwilAssem/output/diamond/Dwil_contig_uniref;
source activate blobtools_env;
/LUSTRE/usuario/pgonzale/software/diamond-0.9.17/diamond view --daa /LUSTRE/usuario/pgonzale/DwilAssem/output/diamond/Dwil_contig_uniref.daa > /LUSTRE/usuario/pgonzale/DwilAssem/output/diamond/Dwil_contig_uniref.diamond
blobtools taxify --map_col_sseqid 0 --map_col_taxid 2 -f /LUSTRE/usuario/pgonzale/DwilAssem/output/diamond/Dwil_contig_uniref.diamond -m /LUSTRE/usuario/pgonzale/dbs/diamond/uniprot_ref_proteomes.taxids -o /LUSTRE/usuario/pgonzale/DwilAssem/output/blobtools/"
echo $cmd | qsub -N diamond -l nodes=1:ppn=20,walltime=99:00:00,vmem=50gb -q default

blobtools create -i /pub/pmgonza1/Dwil/output/platanus/rundef/Dwil_contig.fa -b /pub/pmgonza1/Dwil/output/bwa/platanusDef/contigs.sorted.bam -t /pub/pmgonza1/Dwil/output/blobtools/platanusDef/Dwil_contig_uniref.diamond.taxified.out -t /pub/pmgonza1/Dwil/output/blast/platanusDef/Dwil_contig.blastout -x bestsumorder 
blobtools view -i blobDB.json --hits --rank all -x bestsumorder
blobtools blobplot -i blobDB.json -x bestsumorder

### Quickmerge ###
module load quickmerge/2017.08.3 MUMmer/4.0.0b

cd /bio/pmgonza1/mnark/output/quickmerge/canuNplatnus

merge_wrapper.py --clean_only Dwil_contig.fa canuDwil.fa
mv self_oneline.fa canuDwil.qmFormat.fa
mv hybrid_oneline.fa Dwil_contig.qmFormat.fa
nucmer -t 64 -l 100 -p qm2 Dwil_contig.qmFormat.fa canuDwil.qmFormat.fa
delta-filter -i 95 -r -q qm2.delta > qm2.rq.delta
quickmerge -d qm2.rq.delta -q canuDwil.qmFormat.fa -r Dwil_contig.qmFormat.fa -hco 5.0 -c 1.5 -l n -ml m
mv merged.fasta qm2.def.fasta

### Pilon polishing ###

genome=/bio/pmgonza1/mnark/output/quickmerge/canuNplatnus/qm2.def.fasta
reads=/pub/pmgonza1/Dwil/output/trimmomatic/interleaved/Dwil_illmn_trimm.fq

iter=1
echo -e "#$ -N bwa$iter\n#$ -q bio,free*,pub*\n#$ -pe one-node-mpi 64\n#$ -V\n\ncd /bio/pmgonza1/mnark/output/quickmerge/canuNplatnus/\nbwa index $genome\nbwa mem -t 64 $genome -p $reads | samtools view -@ 64 -bS - | samtools sort -@ 64 -o pilon$iter.bam -\nsamtools index -@ 64 pilon$iter.bam" > bwa$iter.sh
qsub bwa$iter.sh

echo -e "#$ -N pilon$iter\n#$ -q bio,free*,pub*\n#$ -pe one-node-mpi 64\n#$ -V\n\ncd /bio/pmgonza1/mnark/output/quickmerge/canuNplatnus/\npilon --genome $genome --frags pilon$iter.bam --output pilon$iter --threads 64 --changes" > pilon$iter.sh
qsub -hold_jid 5362017 pilon$iter.sh
