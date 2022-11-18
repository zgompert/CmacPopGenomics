# CmacPopGenomics
Genomic analysis of host adaptation in *Callosobruchus maculatus*. Various data sets included here. This is mostly for preliminary results for a NSF proposal.

# PoolSeq data 

See [cmac_pool_DNA_plates](https://drive.google.com/drive/folders/15dMkPsn46LDWEHZt-oxzHbwVGTYr3PV3?usp=share_link) for details on samples and sample concentrations. The table below provides a summary of samples.

| Sample | Sample size|
|--------|------------|
| BF1F1C | 20 |
| BZ1F1C | 20 |
| CA1F1C | 20 |
|BF2F20L | 20 |
|BF4F20L | 20 |
|BF5F20L | 20 |
|BZ4F20L | 20 |
|BZ5F20L | 20 |
|CA4F20L | 20 |
|CA5F20L | 20 |

The pool-seq data are currently in /uufs/chpc.utah.edu/common/home/gompert-group2/data/cmac_poolseq/. The data were generated and cleaned up by BGI. Reports for each round of sequencing are attached: [BGI_F22FTSUSAT0771_ANImdrcR_report_en.pdf](https://github.com/zgompert/CmacPopGenomics/files/9954230/BGI_F22FTSUSAT0771_ANImdrcR_report_en.pdf).

# PoolSeq Alignment and Variant Calling

I am aligning the DNA sequence data to the *C. maculatus* genome on [NCBI](https://www.ncbi.nlm.nih.gov/genome/?term=txid64391[Organism:noexp]), which comes from [Sayadi et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31740847/). I am using bwa-mem2 for this, which is basically just a sped up version of bwa mem that also works directly with gzipped files https://github.com/bwa-mem2/bwa-mem2. I am using bwa-mem2 version 2.0pre2. My overall approach roughly follows the suggested workflow for `samtools` [as described here](https://www.htslib.org/workflow/wgs-call.html).

First, I (re)indexed the reference genome.

```bash
## index genome with bwa-mem2
/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index /uufs/chpc.utah.edu/common/home/gompert-group1/data/cmac_qtl_AR/goran_genome/reference/CAACVG01.fasta
```

Then, I set up the alignment. The submission scrip is (from /uufs/chpc.utah.edu/common/home/gompert-group2/data/cmac_poolseq/Scripts):

```bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=bwa-mem2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)


cd /uufs/chpc.utah.edu/common/home/gompert-group2/data/cmac_poolseq/Alignment

perl ../Scripts/BwaMemFork.pl ../F22FTSUSAT0771_ANImdrcR/soapnuke/clean/*/*1.fq.gz 
```

Which runs the following:

```perl
#!/usr/bin/perl
#
# alignment with bwa mem 
#


use Parallel::ForkManager;
my $max = 12;
my $pm = Parallel::ForkManager->new($max);
my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group1/data/cmac_qtl_AR/goran_genome/reference/CAACVG01.fasta";

FILES:
foreach $fq1 (@ARGV){
	$pm->start and next FILES; ## fork
	$fq2 = $fq1;
	$fq2 =~ s/_1\.fq\.gz/_2.fq.gz/ or die "failed substitution for $fq1\n";
        $fq1 =~ m/clean\/([A-Za-z0-9]+)/ or die "failed to match id $fq1\n";
	$ind = $1;
	$fq1 =~ m/([A-Za-z_\-0-9]+)_1\.fq\.gz$/ or die "failed match for file $fq1\n";
	$file = $1;
        system "/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 1 -k 19 -r 1.5 -R \'\@RG\\tID:Cmac-"."$ind\\tLB:Cmac-"."$ind\\tSM:Cmac-"."$ind"."\' $genome $fq1 $fq2 | samtools sort -@ 2 -O BAM -o $ind"."_$file.bam - && samtools index -@ 2 $ind"."_$file.bam\n";

	$pm->finish;
}

$pm->wait_all_children;
```
I am pipping the results on to samtools (version 1.16) to compress, sort and index the alignments.

Next, I am name sorting the alignments and fixing read pairing information with `samtools`. I am doing this to get things ready for indel realignment and removing duplicates (my next steps). For this, I ran:

```bash
#!/bin/sh
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=wolf-kp
#SBATCH --partition=wolf-kp
#SBATCH --job-name=samt
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)

cd /scratch/general/nfs1/cmac

perl /uufs/chpc.utah.edu/common/home/gompert-group2/data/cmac_poolseq/Scripts/FixMatesFork.pl *bam
```

Which runs the following:

```perl
#!/usr/bin/perl
#
# fix mates, sort and index with samtools
#


use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $bam (@ARGV){
        $pm->start and next FILES; ## fork
        $bam =~ m/^([A-Za-z0-9]+)/ or die "failed to match $bam\n";
        $base = $1;
        system "samtools collate -o co_$base.bam $bam /scratch/general/nfs1/cmac/t$bam\n";
        system "samtools fixmate -m co_$base.bam fix_$base.bam\n";
        system "samtools sort -o sort_$base.bam fix_$base.bam\n";
        $pm->finish;
}

$pm->wait_all_children;
```

Note that this is being done in a scratch directory. I will move back key files for long-term storage.

The next step was to perform indel realignment. Note, this is not available in `GATK` version 4 as this is done as part of variant calling. However, my plan is to call variants with bcftools (much faster) and thus I do want to perform a separate indel realignment. I am using `GATK` version 3.8 for this. 

As a precoursos to indel realignment I created a dict reference file for the genome.

```bash
## dict file with Picard tools 2.22.0
java -jar /uufs/chpc.utah.edu/sys/installdir/picard/2.22.0/picard.jar CreateSequenceDictionary R=/uufs/chpc.utah.edu/common/home/gompert-group1/data/cmac_qtl_AR/goran_genome/reference/CAACVG01.fasta
```
Next, I used `samtools` to index the sorted bam files and `GATK` for indel realignment.

```bash
#!/bin/sh
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=gatkre
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)
module load gatk
## version 3.8

cd /scratch/general/nfs1/cmac

perl /uufs/chpc.utah.edu/common/home/gompert-group2/data/cmac_poolseq/Scripts/RealnFork.pl sort_*bam
```

```perl
#!/usr/bin/perl
#
# indel realignment with GATK version v3.8-0-ge9d806836 
# using this because GATK4 does not have indel realignment 
#


use Parallel::ForkManager;
my $max = 12;
my $pm = Parallel::ForkManager->new($max);

my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group1/data/cmac_qtl_AR/goran_genome/reference/CAACVG01.fasta";

FILES:
foreach $bam (@ARGV){
	$pm->start and next FILES; ## fork
	$bam =~ m/sort_([A-Za-z0-9]+)/ or die "failed to match $bam\n";
	$base = $1;
	#system "samtools index -b $bam\n";
	print "gatk -T RealignerTargetCreator  -I $bam -o $base.intervals -R $genome &\n";
	#print "gatk -T IndelRealigner -R $genome -targetIntervals $base.intervals -I $bam -o realn_$base.bam\n";
	
	$pm->finish;
}

$pm->wait_all_children;
```

