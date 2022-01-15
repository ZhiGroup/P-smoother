#!/bin/bash

if [ $# -eq 0 ]; then
	echo "Program: P-smoother"
	echo ""
	echo "Contact: Degui Zhi [degui.zhi@uth.tmc.edu] or Shaojie Zhang [shzhang@cs.ucf.edu]"
	echo ""
	echo "Usage: ./P-smoother.sh [options] parameters"
	echo ""
	echo "Required Parameters:"
	echo -e "\t--inputVCF <file>\tVCF file"
	echo -e "\t--map <file>\t\tGenetic mapping file"
	echo ""
	echo "Optional Parameters:"
	echo -e "\t--writeTo <filename>\tOutput filename and location (parameter can be full file path or just filename) [Default = VCF filename]"
	echo -e "\t--length <integer>\tMismatch correcting block length (in units of sites) [Default = 20]"
	echo -e "\t--width <integer>\tNumber of haplotypes in block [Default = 20]"
	echo -e "\t--gap <integer>\t\tGap size (in units of site) [Default = 1]"
	echo -e "\t--rho <float>\t\tMinimum allele frequency threshold for mismatch correction [Default = 0.05]"
	echo -e "\t--checkpoint <integer>\tConsole output every n sites [Default = 100000]"
	exit 1
fi

OPTIONS=c:i:o:l:w:g:m:r:
LONGOPTS=checkpoint:,inputVCF:,writeTo:,length:,width:,gap:,map:,rho:

PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
eval set -- "$PARSED" 

checkpoint=100000 inputVCF="" writeTo="" length=20 width=20 gap=1 map="" rho=0.05
while true; do
	case "$1" in
		-c|--checkpoint)
			checkpoint="$2"
			shift 2
			;;
		-i|--inputVCF)
			inputVCF="$2"
			shift 2
			;;
		-o|--writeTo)
			writeTo="$2"
			shift 2
			;;
		-l|--length)
			length="$2"
			shift 2
			;;
		-w|--width)
			width="$2"
			shift 2
			;;
		-g|--gap)
			gap="$2"
			shift 2
			;;
		-m|--map)
			map="$2"
			shift 2
			;;
		-r|--rho)
			rho="$2"
			shift 2
			;;
		--)
			shift
			break
			;;
	esac
done

if [ "$inputVCF" == "" ]; then
	echo "The VCF input file must be specified with the required option --inputVCF <file>"
	exit 1
fi

if [ "$map" == "" ]; then
	echo "The genetic mapping file must be specified with the required option --map <file>"
	exit 1
fi

basename=$(basename $inputVCF)
filename="${basename%.*}"
if [ "$writeTo" = "" ]; then
	writeTo="$filename"
fi

echo "Running rPBWT..."
./rPBWT "$inputVCF" "$writeTo" "$checkpoint"
echo "Running PBWT..."
./PBWT "$inputVCF" "$writeTo" "$map" "$checkpoint" "$length" "$width" "$gap" "$rho"
rm "${writeTo}.rpbwt" "${writeTo}.meta"
echo "P-smoother Finished."

