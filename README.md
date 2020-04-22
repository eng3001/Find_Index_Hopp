# Find_Index_Hopp

***Author:*** Wyatt Eng

***Goal:*** This script is designed to identify index swapped records from the 'Unidentified' output files post demultiplexing.

### Required Input
- -U1 : The path to a fastq file containing the undetermined reads from R1.
- -U2 : The path to a fastq file containing the undetermined reads from R2.
- -I : File of known indexes. With each index separated by a new line.

Example Run:
```python3 Sort_Undetermined_V3.py -U1 Undetermined_S0_L008_R1_001.fastq.gz -U2 Undetermined_S0_L008_R2_001.fastq.gz -I seq_indexes.txt```