RetroSeq
====================
What is RetroSeq?
-------
RetroSeq is a tool for discovery and genotyping of transposable element variants (TEVs) (also known as mobile element insertions) from next-gen sequencing reads aligned to a reference genome in BAM or CRAM format. The goal is to call TEVs that are not present in the reference genome but present in the sample that has been sequenced. It should be noted that RetroSeq can be used to locate any class of viral insertion in any species where whole-genome sequencing data with a suitable reference genome is available.

If you want to find if a transposable element that is present in the reference genome and not present in your sample (a deletion in structural variation terminology), then you should use a structural variation deletion caller such as [Breakdancer](http://gmt.genome.wustl.edu/breakdancer/current/).

RetroSeq is a two phase process, the first being the read pair discovery phase where discorandant mate pairs are detected and assigned to a TE class (Alu, SINE, LINE, etc.) by using either the annotated TE elements in the reference and/or aligned with Exonerate to the supplied library of viral sequences.

Dependencies
-------------
RetroSeq makes use of the [bedtools package](http://code.google.com/p/bedtools/), [samtools](https://github.com/samtools/samtools), [Exonerate](http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate), and various common unix tools such as sort and grep. RetroSeq will check if these are installed and available in the default path.

NOTE: RetroSeq has been primarily tested on BAM files produced by [MAQ](http://maq.sourceforge.net/) and [BWA](http://bio-bwa.sourceforge.net/). There is no guarantee it will work on BAM files from other aligners.

Using RetroSeq
================
1 Discovery Phase
------------------
The goal here is to pass through the BAM (or CRAM) and identify discordant read pairs that might support a TE insertion. You can either supply a tab delimited file specifying a set of TE types (e.g. Alu, LINE etc.) and the corresponding BED file of locations where these are in the reference genome (-refTEs parameter). Alternatively, you can provide a tab delimited file specifying a set of viral/TE types and the corresponding fasta file with a set of consensus sequences for these (-eref parameter).

Usage: retroseq.pl -discover -bam <string> -eref <string> -output <string> [-srmode] [-q <int>] [-id <int>] [-len <int> -noclean]
    
    -bam        BAM or CRAM file of paired reads mapped to reference genome
    -output     Output file to store candidate supporting reads (required for calling step)
    [-eref      Tab file with list of transposon types and the corresponding fasta file of reference sequences (e.g. SINE   /home/me/refs/SINE.fasta). Required when the -align option is used.]
    [-refTEs    Tab file with TE type and BED file of reference elements. These will be used to quickly assign discordant reads the TE types and avoid alignment. Using this will speed up discovery dramatically.]
    [-noclean   Do not remove intermediate output files. Default is to cleanup.]
    [-q         Minimum mapping quality for a read mate that anchors the insertion call. Default is 30. Parameter is optional.]
    [-id        Minimum percent ID for a match of a read to the transposon references. Default is 90.]
    [-len       Minimum length of a hit to the transposon references. Default is 36bp.]
    [-rgs       Comma separated list of readgroups to operate on. Default is all.]
    [-exd       Fofn of BED files of regions where discordant mates falling into will be excluded e.g. simple repeats, centromeres, telomeres]
    [-align     Do the computational expensive exonerate PE discordant mate alignment step]

The discovery stage will produce an output file (specified by the -output parameter) with lists of read pair names per TE type. Multiple discovery stage outputs can be inputted to the calling phase (i.e. to allow parallelisation of the discovery phase).

2 Calling
----------
The calling phase takes one or more outputs from the discovery phase, clusters the reads, and carries out various checks on the breakpoints to decide if a TEV is present. You can provide a list of locations to ignore per TE type - this would typically be the list of locations of the reference elements of that type (-filter option).

Usage: retroseq.pl -call -bam <string> -input <string> -ref <string> -output <string> [-srinput <SR candidates file> -filter <BED file> -cleanup -reads <int> -depth <int> -hets]
    
    -bam            BAM file OR cram OR BAM fofn
    -input          Either a single output file from the PE discover stage OR a prefix of a set of files from discovery to be combined for calling OR a fofn of discovery stage output files
    -ref            Fasta of reference genome
    -output         Output file name (VCF)
    [-hets          Call heterozygous insertions. Default is homozygous.]
    [-filter        Tab file with TE type and BED file of reference elements. These will be filtered out from the calling.]
    [-region        Call a particular chromosome only (chr) OR region (chr:start-end) only]
    [-depth         Max average depth of a region to be considered for calling. Default is 200.]
    [-reads         It is the minimum number of reads required to make a call. Default is 5. Parameter is optional.]
    [-q             Minimum mapping quality for a read mate that anchors the insertion call. Default is 30. Parameter is optional.]
    [-ignoreRGs     Single read group name OR a file of readgroups that should be ignored. Default is none.]
    [-noclean       Do not remove intermediate output files. Default is to cleanup.]

Output
-------
The final TE calls from RetroSeq are in [VCF format](http://vcftools.sourceforge.net/). The calls are annotated with information on number of supporting reads (GQ tag). The FL tag ranges from 1-8 and gives information on the breakpoint with 8 being the most confident calls and lower values indicating calls that don’t meet the breakpoint criteria for reasons such as lack of 5’ or 3’ reads.

More information on the VCF format and tags is on the [tutorial page](https://github.com/tk2/RetroSeq/wiki/RetroSeq-Tutorial)

Human CEU Trio Example
-----------------------
There is a page on the [wiki](https://github.com/tk2/RetroSeq/wiki/1000-Genome-CEU-Trio-Analysis) describing how the CEU human trio was analysed in our [paper](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/bts697?ijkey=TV8gQVdMmWsWM7n&keytype=ref)

Citing RetroSeq
---------------
If you find RetroSeq useful, please cite our [paper](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/bts697?ijkey=TV8gQVdMmWsWM7n&keytype=ref).

Questions or help
-----------------
For questions or support contact Thomas Keane (tk2 --at-- sanger.ac.uk) 
