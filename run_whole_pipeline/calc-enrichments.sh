#!/bin/bash

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

# Cluster sequences seen by each protein
for all_seq_file in *_all_pooled.counts; do
    out_file=`echo ${all_seq_file} | sed 's/.counts/.cluster/'`
    echo "Clustering ${out_file}"
    awk '{print $1}' ${all_seq_file} > tmp
    hops_cluster tmp -s 2 -e 1 -o ${out_file}
    rm -f tmp
done

# Calculate enrichments
echo "Calculating hA5 enrichments"
hops_enrich hA5_conv_1.counts hA5_comp_1.counts -f hA5_all_pooled.cluster -m 6 -b hA5_1
hops_enrich hA5_conv_2.counts hA5_comp_2.counts -f hA5_all_pooled.cluster -m 6 -b hA5_2
hops_enrich hA5_conv_pooled.counts hA5_comp_pooled.counts -f hA5_all_pooled.cluster -m 6 -b hA5_pooled

echo "Calculating hA6 enrichments"
hops_enrich hA6_conv_1.counts hA6_comp_1.counts -f hA6_all_pooled.cluster -m 6 -b hA6_1
hops_enrich hA6_conv_2.counts hA6_comp_2.counts -f hA6_all_pooled.cluster -m 6 -b hA6_2
hops_enrich hA6_conv_pooled.counts hA6_comp_pooled.counts -f hA6_all_pooled.cluster -m 6 -b hA6_pooled

echo "Calculating aA5A6 enrichments"
hops_enrich aA5A6_conv_1.counts aA5A6_comp_1.counts -f aA5A6_all_pooled.cluster -m 6 -b aA5A6_1
hops_enrich aA5A6_conv_2.counts aA5A6_comp_2.counts -f aA5A6_all_pooled.cluster -m 6 -b aA5A6_2
hops_enrich aA5A6_conv_pooled.counts aA5A6_comp_pooled.counts -f aA5A6_all_pooled.cluster -m 6 -b aA5A6_pooled

echo "Calculating alt enrichments"
hops_enrich alt_conv_1.counts alt_comp_1.counts -f alt_all_pooled.cluster -m 6 -b alt_1
hops_enrich alt_conv_2.counts alt_comp_2.counts -f alt_all_pooled.cluster -m 6 -b alt_2
hops_enrich alt_conv_pooled.counts alt_comp_pooled.counts -f alt_all_pooled.cluster -m 6 -b alt_pooled
