#!/bin/sh

upint=5000
downint=5000
query=EFB0058.fa
genomespath=""
samplename=""
run_all=true
step=""

while getopts ":u:d:q:g:n:s:h" opt; do
	case $opt in
		u)
			upint="$OPTARG"
			;;
		d)
			downint="$OPTARG"
			;;
		q)
			query="$OPTARG"
			;;
		g)
			genomespath="$OPTARG"
			;;
		h)
			echo "usage:  ./RADS.sh -u [upstream integer (optional, default = 5000)] -d [downstream integer(optional, default = 5000)] -q [fasta query file (optional, default = EFB0058.fa)] -g [path to directory containing directories with genome information. (Required. Often ncbi_dataset/data/)] -n [sample name to appended to all file and directory names. (Optional)] -s [step to run individually. (Optional)]"
			echo ""
			echo "[REQUIRED OPTIONS]"
			echo "-------------------------------------------------------------"
			echo "-g|--genomespath	path to directory containing directories with genomes information. Each genome should have its own directory containing .fna files"
			echo ""
			echo "[OPTIONAL ARGUMENTS]"
			echo "------------------------------------------------------------"
			echo "-h|--help		print help message and exit"
            		echo "-u|--upint	integer for number of nucleotides upstream of ORF to extract. Default=5000 if not specified"
			echo "-d|--downint	integer for number of nucleotides downstream of ORF to extract. Default=5000 if not specified"
			echo "-q|--query	amino acid fasta (.faa/.fa) file for ORF(s) to use as query. Default=EFB0058.fa if not specified"
			echo "-n|--name		sample name to be added to all output files and directories"
			echo "-s|--step		run an individual step of the pipeline. Options:"
            echo "         		1|mvgenomes --> move genomes into genomes_${samplename}"
            echo "         		2|translate --> translate genomic sequences into amino acid sequences"
			echo "				3|makedbs --> make diamond databases from amino acid sequences"
			echo "				4|blast --> blast for query in diamond databases"
			echo "				5|extractcontigs --> extract contigs of specified size (default 10kb) around query"
			echo "				6|contigorfprocessing --> search for ORFs and run InterProScan on contigs around query"
			echo "				7|cotranscription --> search for putatively cotranscribed genes around query"
            		exit 1
			;;
		n)
			if [ -n "$OPTARG" ]; then  # Check if an argument was provided to -n
                                samplename="$OPTARG"
                        fi
                        ;;
        	s)
            		step="$OPTARG"
            		run_all=false
            		shift 2
            		;;
		\?)
			echo "option -$OPTARG requires an argument. Usage: bash ./RADS.sh -u [upstream integer (optional, default = 5000)] -d [downstream integer(optional, default = 5000)] -q [fasta query file (optional, default = EFB0058.fa)] -g [path to directory containing directories with genome information. Often ncbi_dataset/data/)" >&2
			exit 1
			;;
		:)
			echo "option -$OPTARG requires an argument. Usage: bash ./RADS.sh -u [upstream integer (optional, default = 5000)] -d [downstream integer(optional, default = 5000)] -q [fasta query file (optional, default = EFB0058.fa)] -g [path to directory containing directories with genome information. Often ncbi_dataset/data/)" >&2
			exit 1
			;;
	esac
done


cp -r ${genomespath} genomes_${samplename}/

translate( ){
    mkdir genomes_translated_${samplename}/
    for i in $(ls genomes_${samplename}/*.fna)
    do
        prodigal \
        -p meta \
        -i ${i} \
        -o genomes_translated_${samplename}/${i}prodigal.txt \
        -a genomes_translated_${samplename}/${i}translated.faa
    done
    echo "genomes translated"
}

makedbs( ){
    mkdir diamonddbs_${samplename}/
    for i in $(ls genomes_translated_${samplename}/*.faa)
    do
        diamond makedb \
        --in ${i} \
        --db diamonddbs_${samplename}/${i}.db
    done
}

blast( ){
    mkdir blast_results_30_${samplename}/
    for i in $(ls diamonddbs_${samplename})
    do
        if [ ! -f blast_results_30_${samplename}/EFB0058_blast_${i}.txt ]; then
	        diamond blastp \
			-d ${i} \
			--query $query \
			--out blast_results_30_$samplename/EFB0058_blast_${i}.txt \
			--outfmt 6 qseqid sseqid length nident \
			--max-target-seqs 0 \
			--id 30
        fi
    done

    for i in $(ls blast_results_30_${samplename}/)
    do
        cat blast_results_30_${samplename}/${i} >> master${samplename}.txt
    done
}

extractcontigs ( ){
    mkdir EFB0058_ORFS_$samplename/
    mkdir EFB0058_flanks_$samplename/
    mkdir bedfiles_$samplename/
    mkdir contigs_$samplename/
    mkdir EFB0058_formatted_coordinates_$samplename/
    mkdir EFB0058_final_coordinates_$samplename/


    for i in $(ls genomes_${samplename}/)
    do
        if [ -s blast_results_30_${samplename}/EFB0058_blast_${i}translated.faa.db.dmnd.txt ]; then
            echo ${i}
            
			cut -f2 blast_results_30_${samplename}/EFB0058_blast_${i}translated.faa.db.dmnd.txt > EFB0058_ORFS_${samplename}/${i}_EFB0058_ORFs.txt
            echo "ORF IDs extracted"
            
			seqkit grep -f EFB0058_ORFS_${samplename}/${i}_EFB0058_ORFs.txt genomes_translated_${samplename}/${i}translated.faa | seqkit seq -n > EFB0058_ORFS_${samplename}/${i}_EFB0058_hits_coordinates.txt
            echo "Coordinates extracted"
            
			awk 'BEGIN{OFS="\t"} {F="#"} {up=$3; down=$5; if (up>down) print down, up; else print up, down}' EFB0058_ORFS_${samplename}/${i}_EFB0058_hits_coordinates.txt > EFB0058_formatted_coordinates_${samplename}/${i}_EFB0058_formatted_coordinates.txt
            awk -v upint=$upint -v downint=$downint 'BEGIN{OFS="\t"} {up=$1-upint; down=$2+downint; if (up>0) print up, down; else print 0, down}' EFB0058_formatted_coordinates_${samplename}/${i}_EFB0058_formatted_coordinates.txt > EFB0058_final_coordinates_${samplename}/${i}_EFB0058_final_coordinates.txt
            echo "Coordinates formatted"
            
			cut -f 2 blast_results_30_${samplename}/EFB0058_blast_${i}translated.faa.db.dmnd.txt | cut -f 1 -d "_" > EFB0058_flanks_${samplename}/${i}_EFB0058_contigs.txt
            echo "Contig names extracted"
            
			paste EFB0058_flanks_${samplename}/${i}_EFB0058_contigs.txt EFB0058_final_coordinates_${samplename}/${i}_EFB0058_final_coordinates.txt > bedfiles_${samplename}/${i}.bed
            echo "Bed file created"
            
			seqkit subseq --bed bedfiles_${samplename}/${i}.bed genomes_${samplename}/${i} > contigs_${samplename}/${i}_EFB0058_contigs.fna
            echo "Contigs created for ${i} with flanks $upint up and $downint down"
        else
            echo "File ${i} is empty, skipping..."
        fi
    done
}

contigorfprocessing( ){
    for i in $(contigs_${samplename}/*.fna)
    do
        cat contigs_${samplename}/${i} >> allcontigsconcatenated_${samplename}.fna
    done

    seqkit seq -m 50 allcontigsconcatenated_${samplename}.fna > allcontigsconcatenated_filtered_${samplename}.fna
    prodigal \
    -p meta \
    -i allcontigsconcatenated_filtered_${samplename}.fna \
    -o allcontigsconcatenated_${samplename}.txt \
    -a allcontigsconcatenated_${samplename}.faa
    sed 's/*//g' allcontigsconcatenated_${samplename}.faa > interposcaninput_${samplename}.faa
    bash my_interproscan/interproscan-*/interproscan.sh -i interposcaninput_${samplename}.faa
}


if [ -z "$genomespath" ]; then
	echo "Error: Option -g (genomes path) is required." >&2
	exit 1
fi

if [ -z "$upint" ]; then
	echo "Error: Option -u (nucleotides upstream) is required as an integer." >&2
	exit 1
fi

if [ -z "$downint" ]; then
	echo "Error: option -d (nucleotides downstream) is required as an integer." >&2
	exit 1
fi

if $run_all; then
    echo "Thank you for using RADS! RADS will now run with the following parameters:"
    echo "upstream set to: $upint"
    echo "downstream set to: $downint"
    echo "query set to: $query"
    echo "genomes path set to: $genomespath"
    echo "sample name set to: $samplename"
    echo "running entire pipeline..."
    translate
    makedbs
    extractcontigs
    contigorfprocessing
else
    case "$step" in
      2|translate)
		    echo "Thank you for using RADS! RADS will now translate input genomes..."
        translate ;;
	    3|makedbs)
		    echo "Thank you for using RADS! RADS will now make Diamond databases for input genomes..."
		    makedbs ;;
        4|blast)
		    echo "Thank you for using RADS! RADS will now BLAST for your query sequence..."
		    echo "query set to: $query"
		    blast;;
        5|extractcontigs)
		    echo "Thank you for using RADS! RADS will now extract contigs surrounding your query sequence..."
		    echo "upstream set to: $upint"
		    echo "downstream set to: $downint"
		    extractcontigs;;
        6|contigorfprocessing)
		    echo "Thank you for using RADS! RADS will now process your contigs..."
		    contigorfprocessing;;
        *) echo "Error: invalid step argument";;
    esac
fi
