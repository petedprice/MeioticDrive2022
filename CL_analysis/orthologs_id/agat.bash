source  activate agat

agat_sp_keep_longest_isoform.pl -gff $1/$2.gtf -o $1/${2}_longest.gtf

agat_sp_extract_sequences.pl -gff $1/${2}_longest.gtf --fasta $1/$2.fasta -t cds -o $1/${2}_longest.fasta

