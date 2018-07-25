import pysam
import pandas as pd
import sys
import csv

## get the filename from the command line input
if(len(sys.argv) != 2):
    print("Requires 1 argument: path to the vcf file created from write_variant_file.py")
    sys.exit()
else:
    vcf_file = sys.argv[1]


## read in the variant file
with pysam.VariantFile(vcf_file) as variant_file:
    variants = [v for v in variant_file.fetch()]

## put the pertinent information from the vcf file into a dataframe
variants_arr = [
    ("%s:%d" % (rec.chrom, rec.pos),
     rec.samples[i].name,
     rec.samples[i]['AO'][0],
     rec.samples[i]['RO'],
     rec.samples[i]['DP']
    ) 

    for rec in variants
    for i in range(0, len(rec.samples))
];
variants_df = pd.DataFrame.from_records(variants_arr, columns=['region', 'cell_barcode', 'alt_count', 'ref_count', 'depth']);
variants_df = variants_df.fillna(0)


## predict the tumor cells with a simple cutoff algorithm
alt_sum = variants_df.groupby('cell_barcode')['alt_count'].agg('sum')
depth_sum = variants_df.groupby('cell_barcode')['depth'].agg('sum')
positive_ratio = alt_sum/depth_sum

info_ret_tumor = positive_ratio[positive_ratio >= .1]
info_ret_normal = positive_ratio[positive_ratio <= .1]
tumor_barcodes_pred = info_ret_tumor.index    
normal_barcodes_pred = info_ret_normal.index


## write the tumor barcodes to a csv file called tumor_barcodes.csv
with open('tumor_barcodes.csv', 'w') as tumor_file:
    for barcode in tumor_barcodes_pred:
        tumor_file.write(barcode + "\n")

## write the normal barcodes to a csv file called normal_barcodes.csv
with open('normal_barcodes.csv', 'w') as normal_file:
    for barcode in normal_barcodes_pred:
        normal_file.write(barcode + "\n")
