fasta_formatter -i data/pATLAS_download_3071981314485723071.fas -o data/pATLAS_download_3071981314485723071_formatted.fas 
paste - - < data/pATLAS_download_3071981314485723071_formatted.fas | shuf | head -n 10  | tr '\t' '\n' > data/pATLAS_random_10.fasta
mkdir -p data/split_plasmids
rm data/split_plasmids/*
awk '/^>/{split($1,a,">")}{print >> "data/split_plasmids/"a[2]".fa"}' < data/pATLAS_random_10.fasta
