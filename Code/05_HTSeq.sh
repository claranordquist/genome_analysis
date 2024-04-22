# Need to take away the fasta stuff at the end of the gff files for each bin

# Syntax: htseq-count [options] <alignment_files> <gff_file>
# Options that need to be specified are:
# -f bam That the input files are in bam format, default is sam
# -r pos That the bam files are sorted according to coordinates (from samtools sort), default is by name. It's because we have paired-end reads
# -t CDS That the feature type (3rd column in the gff file) is CDS, default is exon
# -o <output file> The name of the (sam-) file to which the program will write the output
