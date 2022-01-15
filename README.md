# P-smoother 

## Introduction
P-smoother is an efficient method to correct recent mutations and genotyping errors in biobank-scale haplotype panels. Utilizing P-smoother to "smooth" a panel allows for downstream haplotype matching tasks to be error tolerant and more accurate. The input data for P-smoother is phased genotype data (in VCF format).

## Dependencies
- C++ (at least GCC 5)  
- GNU Make  

## Installation
To install the program clone the repository to a local folder using:

`git clone https://github.com/ZhiGroup/P-smoother.git`

Enter the repository folder and compile the program:

`cd P-smoother`  
`make`

## Usage Instructions
Type  
`./P-smoother.sh`  
by itself to show the help page.  

|         Option         |                         Parameter                        |                                                    Description                                                    |
|:----------------------:|:--------------------------------------------------------:|:-----------------------------------------------------------------------------------------------------------------:|
| `-i` or `--inputVCF`    | Full or relative file path to input VCF file             | Sets the input VCF file on which to run biPBWT.                                                                   |
| `-m` or `--map`    	 | Full or relative file path to Genetic Mapping file       | Sets the Genetic Mapping file.                                                                                    |
| `-o` or `--writeTo`            | Full or relative file path and filename for output files | Sets the location and filename for biPBWT output files. The default option is the VCF filename.                   |
| `-l` or `--length`     | Mismatch correction block length (in units of sites)                    | Sets the minimum length requirement for mismatch correction blocks. The default value is 20 sites.                          |
| `-w` or `--width`      | Mismatch correction block width                                              | Sets the minimum number of haplotypes required for a mismatch correction block. The default value is 20 haplotypes.   |
| `-g` or `--gap`        | Gap size (in units of sites)                             | Sets the size of the gap in the block. The default value is 1 site.                                               |
| `-c` or `--checkpoint` | Checkpoint interval `n`                                  | biPBWT will print a checkpoint message to console every `n` sites. The default value is 100,000.                  |

An example:  
`./P-smoother.sh --inputVCF "example.vcf" --map "example.map" --writeTo "output" --length 20 --width 20 --gap 1 --checkpoint 100000`  

## Genetic Mapping Format
The format of the genetic mapping file must be 2 tab-separated fields per line: "site number" "genetic mapping".

## Results
When finished executing, P-smoother will generate a VCF file with the extension ".smooth.vcf".

Note that P-smoother generates 2 intermediate files with the extensions ".rpbwt" and ".meta" during execution that are deleted upon completion. The ".rpbwt" file stores the reverse positional prefix array and reverse divergence array in binary format (each value takes 4 bytes) and takes up around four times the disk space of the VCF file. The ".meta" files stores two space-separated values M and N representing the number of haplotypes and the number of sites in the VCF file.

## PS-cluster (P-smoother cluster)
For an efficent multi-way IBD cluster detection algorithm with an integrated P-smoother, see the PS-cluster folder.

## Citation
If you found our work useful in your research, please consider citing the following paper:

