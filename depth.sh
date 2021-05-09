#! /usr/bin/env bash
# depth.sh

# Process all cram files in a directory
# Determine ABO gene sequence depth with samtools depth command
# samtools depth -a -r chr9:133255176-133385146 [*.cram] -o [*.txt]

# path to cram files
cramPath="/Users/jonathan_stevens/ABO/1000G_data/cram/"
# cram file extension
cramSuffix=".extract_for_bloodantigens-full_gene_master.cram"

# path to output directory
depthOutPath="/Users/jonathan_stevens/ABO/1000G_data/depth/"
# output file extension
outSuffix="_ABO_pos_read_depth_of_coverage.txt"

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
        echo $sampleName$outSuffix
    done
}

depth 1>depth.out 2>depth.err &


