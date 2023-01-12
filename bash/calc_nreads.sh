#! /bin/bash -p

bam_dir=$1

shopt -s nullglob
for bam_path in "$bam_dir"/*.bam; do
    bam_file=${bam_path##*/}
    bam_name=${bam_file%.bam}
    echo "Calculating number of reads for $bam_file"
    samtools view -c -F 260 "$bam_file"
done
echo "All done! :)"