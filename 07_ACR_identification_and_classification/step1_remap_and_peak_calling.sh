Follow the remapping and peak-calling workflow in 02_mapping and 03_high_quality_nuclei.

Note:
B73 is no longer used as the universal reference genome.
Map each parental inbred to its own reference genome (B73, Ki3 or Oh43) and each hybrid to the corresponding pseudo-hybrid genome (B73-Ki3, B73-Oh43 or Ki3-Oh43).
The BAM files will then be split by cell type based on barcode and peak calling will be performed separately for each cell type.
