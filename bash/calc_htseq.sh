#! /bin/bash -p

bam_dir=$1
gff_file=$2

outputextension=".count"

shopt -s nullglob
for bam_path in "$bam_dir"/*.bam; do
    bam_file=${bam_path##*/}
    bam_name=${bam_file%.bam}
    echo "Calculating number of reads for $bam_file"
    samtools index "$bam_file"
    htseq-count -f bam -m union -t gene -i ID -a 0 "$bam_file" "$gff_file" > "$bam_name$outputextension"
    echo "saved as "$bam_name$outputextension""
done
echo "All done! :)"