#! /usr/bin/env bash
# depth.sh

# path to cram files
cramPath="/Users/jonathan_stevens/ABO/1000G_data/cram/"
cramSuffix=".extract_for_bloodantigens-full_gene_master.cram"

# path to output
depthOutPath="/Users/jonathan_stevens/ABO/1000G_data/depth/"
outSuffix="_depth.txt"

# create output directory
mkdir -p $depthOutPath

# loop through all cram files in $cramPath
function depth {
    for cram in $cramPath*$cramSuffix
    do
        # Remove path from filename
        pathRemoved="${cram/$cramPath/}"
        # Remove suffix to get file name
        sampleName="${pathRemoved/$cramSuffix/}"
        samtools depth -a -r chr9:133255176-133385146 \
        $cramPath$sampleName$cramSuffix \
        -o $depthOutPath$sampleName$outSuffix
    done
}

depth 1>depth.log 2>depth.err &


