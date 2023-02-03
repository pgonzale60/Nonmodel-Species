The genome_assembly.sh script contains the commands to generate a new genome assembly using PacBio long sequencing reads for scaffolding and Illumina short reads for polishing. This pipeline is largely based on that devised by Chakraborty and colleagues in their publication Contiguous and accurate de novo assembly of metazoan genomes with modest long read coverage (doi: 10.1093/nar/gkw654). The underlying structure is organized as listed below. The programs were non-default parameter values were used have these specificied inside parenthesis. 



1. With Trimmommatic v0.32 (`ILLUMINACLIP: adaptersE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:20`), trim Illumina reads to remove adapter sequence, low confidence bases and exclude reads which are too short after trimming.
2. The trimmed reads are used by Platanus to generate a genome assembly.
3. The single pass Continuous Long Reads (CLR) PacBio reads were used by Canu v1.4 (`minOverlapLength=500 rawErrorRate=0.3 correctedErrorRate=0.045 genomeSize=250m`) to generate an assembly with lower nucleotide fidelity but higher contiguity.
4. Optional. The tool [blobtools](https://blobtools.readme.io/docs) can be used to detect and remove sequences assembled from contaminant sources (e.g. bacteria) in the above assemblies.
5. The assembled contigs are then combined with Quickmerge (`delta-filter -i 95`), using the Platanus assembly as ground sequence while the canu assembly provides the linkage information where the assemblies are redundant or the nucleotide sequence where no sequence from the platanus assembly could be assigned. The latest version of quickmerge and the detailed description of its use can be found in the following [link](https://github.com/mahulchak/quickmerge).
6. Polishing steps of the resulting assembly using the trimmed Illumina short reads with Pilon v1.2.
