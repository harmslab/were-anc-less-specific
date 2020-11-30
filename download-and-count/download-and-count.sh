#!/bin/bash 

file_with_sra=${1}
if [ ! "${file_with_sra}" ]; then
    echo "specify a file with sra ids"
    exit
fi

prefetch_binary=`which prefetch`
if [ ! "${prefetch_binary}" ]; then
    echo "could not find prefetch"
    echo "please install: https://github.com/ncbi/sra-tools/"
    exit
fi

hops_count_binary=`which hops_count`
if [ ! "${hops_count_binary}" ]; then
    echo "could not find hops_enrich"
    echo "please install: https://github.com/harmslab/hops_enrich"
    exit
fi

while read line; do
    sra=`echo "${line}" | awk '{print $1}'`
    name=`echo "${line}" | awk '{print $2}'`

    echo ${sra} ${name} 
    
    prefetch ${sra}
    fasterq-dump ${sra}
    mv ${sra}.fastq ${name}.fastq
    gzip ${name}.fastq

done < ${file_with_sra}

# Count peptides in each fastq file, applying quality control
for read_file in *.fastq.gz; do 
    count_file=`echo $read_file | sed 's/.fastq.gz/.counts/'`

    echo Creating ${count_file}
    hops_count ${read_file} -o ${count_file}
done

# Create pooled samples
./pool_replicates.py hA5_conv_1.counts hA5_conv_2.counts hA5_conv_pooled.counts
./pool_replicates.py hA5_comp_1.counts hA5_comp_2.counts hA5_comp_pooled.counts
./pool_replicates.py hA6_conv_1.counts hA6_conv_2.counts hA6_conv_pooled.counts
./pool_replicates.py hA6_comp_1.counts hA6_comp_2.counts hA6_comp_pooled.counts
./pool_replicates.py aA5A6_conv_1.counts aA5A6_conv_2.counts aA5A6_conv_pooled.counts
./pool_replicates.py aA5A6_comp_1.counts aA5A6_comp_2.counts aA5A6_comp_pooled.counts
./pool_replicates.py alt_conv_1.counts alt_conv_2.counts alt_conv_pooled.counts
./pool_replicates.py alt_comp_1.counts alt_comp_2.counts alt_comp_pooled.counts

# Create pools of every peptide seen for each protein to generate clusters
./pool_replicates.py hA5_conv_pooled.counts hA5_comp_pooled.counts hA5_all_pooled.counts
./pool_replicates.py hA6_conv_pooled.counts hA6_comp_pooled.counts hA6_all_pooled.counts
./pool_replicates.py aA5A6_conv_pooled.counts aA5A6_comp_pooled.counts aA5A6_all_pooled.counts
./pool_replicates.py alt_conv_pooled.counts alt_comp_pooled.counts alt_all_pooled.counts
