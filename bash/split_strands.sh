#! /bin/bash -p

bam_dir=$1

plusextension="_forward.bam"
minusextension="_reverse.bam"

shopt -s nullglob
for bam_path in "$bam_dir"/*.bam; do
    bam_file=${bam_path##*/}
    bam_name=${bam_file%.bam}
    echo "Splitting forward strand for $bam_file"
    samtools view -F 16 -o "$bam_name$plusextension" "$bam_file"
    echo "saved as $bam_name$plusextension"
    echo "Splitting reverse strand for $bam_file"
    samtools view -f 16 -o "$bam_name$minusextension" "$bam_file"
    echo "saved as $bam_name$minusextension"

done
echo "All done!"
