# PS-cluster (P-smoother cluster)

## Introduction
PS-cluster is an efficient method to extract multi-way IBD clusters while tolerating recent mutations and genotyping errors. PS-cluster integrates the process of smoothing a panel with P-smoother and then extracting clusters with positional Burrows-Wheeler transform. The input data for PS-cluster is phased genotype data (in VCF format).

## Dependencies
- C++ (at least GCC 5)  
- GNU Make  

## Installation
To install the program clone the repository to a local folder using:

`git clone https://github.com/ZhiGroup/P-smoother.git`

Enter the repository folder and compile the program:

`cd P-smoother/PS-cluster`  
`make`

## Usage Instructions
Type  
`./PS-cluster.sh`  
by itself to show the help page.  

|         Option         |                         Parameter                        |                                                    Description                                                    |
|:----------------------:|:--------------------------------------------------------:|:-----------------------------------------------------------------------------------------------------------------:|
| `-i` or `--inputVCF`    | Full or relative file path to input VCF file             | Sets the input VCF file on which to run PS-cluster.                                                                   |
| `-m` or `--map`    	 | Full or relative file path to Genetic Mapping file       | Sets the Genetic Mapping file.                                                                                    |
| `-o` or `--writeTo`            | Full or relative file path and filename for output files | Sets the location and filename for biPBWT output files. The default option is the VCF filename.                   |
| `-l` or `--length`     | Mismatch correction block length (in units of sites)                    | Sets the minimum length requirement for mismatch correction blocks. The default value is 20 sites.                          |
| `--length_min`     | Minimum target block length (in units of centimorgans)                    | Sets the minimum length requirement for reported blocks. The default value is 1 cM.                          |
| `-w` or `--width`      | Mismatch correction block width                                              | Sets the minimum number of haplotypes required for a mismatch correction block. The default value is 20 haplotypes.   |
| `--width_min`      | Minimum target block width                                              | Sets the minimum number of haplotypes required for a block to be reported. The default value is 10 haplotypes.   |
| `-g` or `--gap`        | Gap size (in units of sites)                             | Sets the size of the gap in the block. The default value is 1 site.                                               |
| `-c` or `--checkpoint` | Checkpoint interval `n`                                  | biPBWT will print a checkpoint message to console every `n` sites. The default value is 100,000.                  |

An example:  
`./PS-cluster.sh --inputVCF "example.vcf" --map "example.map" --writeTo "output" --length 20 --length_min 1 --width 20 --width_min 10 --gap 1 --checkpoint 100000`  

## Genetic Mapping Format
The format of the genetic mapping file must be 2 tab-separated fields per line: "site number" "genetic mapping".

## Results
When finished executing, PS-cluster will generate a file with the extension ".blocks".

The file ".blocks" represents each block on its own line with five space-separated fields `<starting site of block> <ending site of block> <starting physical location of block> <ending physical location of block> <starting genetic location of block> <ending genetic location of block> <width of block>` followed by space seperated ID's of all the haplotypes in the block. IDs are suffixed with either "-0" or "-1" indicating the first and second haplotype of the individual ID, respectively.

Note that PS-cluster generates 3 intermediate files with the extensions ".rpbwt", ".sites", and ".meta" during execution that are deleted upon completion. The ".rpbwt" file stores the reverse positional prefix array and reverse divergence array in binary format (each value takes 4 bytes) and takes up around four times the disk space of the VCF file. The ".sites" file stores the physical position of each site on its own line. The ".meta" files stores two space-separated values M and N representing the number of haplotypes and the number of sites in the VCF file.

## Citation
If you found our work useful in your research, please consider citing the following paper:

