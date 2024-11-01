#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i <fasta_file> -b <busco_file> -a <asm_stats_file> -f <blob_file> -o <output_name>"
    echo
    echo "Options:"
    echo "  -i <fasta_file>      Path to the fasta file"
    echo "  -b <busco_file>      Path to the BUSCO file"
    echo "  -a <asm_stats_file>  Path to the assembly stats file"
    echo "  -f <blob_file>       Path to the blob file"
    echo "  -o <output_name>     Desired name for the output file"
    echo "  -h                   Display this help message"
    exit 1
}

# Check if seqkit is installed
if ! command -v seqkit &> /dev/null; then
    echo "Error: seqkit is not installed. Please install seqkit to proceed."
    exit 1
fi

# Check if no arguments were passed and display usage
if [ $# -eq 0 ]; then
    usage
fi

# Read command-line arguments
while getopts "i:b:a:f:o:h" flag; do
    case "${flag}" in
        i) fasta=${OPTARG};;
        b) busco_file=${OPTARG};;
        a) asm_stats_filename=${OPTARG};;
        f) blob_file=${OPTARG};;
        o) output_name=${OPTARG};;
        h) usage;;
        *) usage;;
    esac
done

# Check if all required arguments are provided
if [ -z "${fasta}" ] || [ -z "${busco_file}" ] || [ -z "${asm_stats_filename}" ] || [ -z "${blob_file}" ] || [ -z "${output_name}" ]; then
    echo "Error: Missing required arguments."
    usage
fi

############## Logging braces open #####################
{
###### Filtering assembly based on contigs size

contigs_siz_fil() {
    local fileA=$1
    local fileB=$2
    
    # Extracting contigs from fileA
    fil_contigs=$(sed '1,/>>>>>>> Coverage per contig/d' "$fileA" | awk -F'\t' 'NR>1 && $3 < 1000 {print $2}')
    
    # contigs to compare (example contigs)
    contigs=($fil_contigs)
    
    # Using awk to filter contigs that do not match the  column of fileB
    unmatched_contigs=$(awk -F'\t' -v contigs="${contigs[*]}" '
    BEGIN {
        #Split the contigs into an array
        split(contigs, elementArray, " ")
        #Create an associative array to store contigs
        for (i in elementArray) {
            comparecontigs[elementArray[i]] = 1
        }
    }
        # Skip lines starting with #
         /^#/ {
         next
         }
    {
        #Check if the element is in the files column
        if ($3 in comparecontigs) {
            delete comparecontigs[$3]
        }
    }
    END {
        #Print contigs that were not found in the file
        for (e in comparecontigs) {
            print e
        }
    }' "$fileB")
    
    echo "$unmatched_contigs"    
}

######### Filtering assembly based on contaminated contigs

contigs_cont_fil() {
    local fileA=$1
    #Extract contaminated contigs
    cont_contigs=$(awk -F "," 'NR > 1 {
    # Remove double quotes from fields
    gsub(/"/, "", $6);
    gsub(/"/, "", $7);
    
    # Check if 6th column is neither "Arthropoda" nor "no-hit"
    if ($6 != "Arthropoda" && $6 != "no-hit") {
        print $7;
        }
    }' "$fileA") 
    
    echo "$cont_contigs"

}


echo -e "File $asm_stats_filename is loaded\n"
echo -e "File $busco_file is loaded\n"
echo -e "File $fasta is loaded\n"
echo -e "File $blob_file is loaded\n"

# Call the size filter function
siz_fil_contigs=$(contigs_siz_fil "$asm_stats_filename" "$busco_file")
#echo "$siz_fil_contigs" #> size_fil_contig.txt

# Call the contaminantion filter function
cont_fil_contigs=$(contigs_cont_fil "$blob_file")
#echo "$cont_fil_contigs"

echo -e "$siz_fil_contigs\n$cont_fil_contigs" > contigs_rmvd.txt

#SeqKit to filter scaffold
seqkit grep -v -f contigs_rmvd.txt $fasta > $output_name

echo -e "Filtered fasta is stored in $output_name\n"

ori_contig_num=$(grep -c "^>" $fasta)
fil_contig_num=$(grep -c "^>" $output_name)

echo -e "Number of contigs\nBefore filteration:$ori_contig_num\nAfter filteration :$fil_contig_num\n"

}  2>&1 | tee -a file.log

####################### Logging braces closed #####################