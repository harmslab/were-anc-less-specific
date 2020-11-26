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

