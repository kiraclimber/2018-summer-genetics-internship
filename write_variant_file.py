## import necessary libraries
import pysam
import sys
import numpy as np
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import itertools

## merges the two dictionaries x and y. if there are duplicate keys, the value from dictionary y is kept
def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z


## filters out cells with a low UMI count.
## bam_file = path to the bam file
## cutoff_val = only cells with a UMI count larger than cutoff_val are kept (optional, should specify either this or
##   cutoff_percentile)
## cutoff_percentile = only cells with a UMI count in the cutoff_percentile are kept (optional, should specify either
##   this or cutoff_val)
## returns the list of barcodes above the cutoff value or percentile, as well as a dictionary of the UMI count for  
##   every barcode
def filter_by_umi(bam_file, cutoff_val, cutoff_percentile):
    barcode_umi = {}
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        for read in samfile.fetch():
            if(read.has_tag('CB') and read.has_tag('UB')):
                barcode = read.get_tag('CB')
                molecule = read.get_tag('UB')
                if barcode in barcode_umi:
                    if(molecule in barcode_umi[barcode]):
                        barcode_umi[barcode][molecule] += 1
                    else:
                        barcode_umi[barcode][molecule] = 1
                else:
                    barcode_umi[barcode] = {molecule:1}
    umi_counts = {}
    for barcode in barcode_umi:
        umi_counts[barcode] = len(barcode_umi[barcode])
    barcodes = []
    if(cutoff_val == 0):
        cutoff_val = np.percentile(umi_counts.values(), cutoff_percentile) #get cutoff_val from the cutoff_percentile
    for barcode in umi_counts:
        if(umi_counts[barcode] > cutoff_val):                
            barcodes.append(barcode)
    
    return barcodes, umi_counts



## This function takes a path to a bam file, and a list of variants, 
## and returns a dataframe that represents the pileup results at those 
## variant positions, with cell_barcode information preserved.
def pileup(bam_fn, variants):
    with pysam.AlignmentFile(bam_fn, "rb") as samfile:   
        pileup_recs = [
            ("%s:%d" % (rec.chrom, rec.pos),
            #("%s:%d" % (samfile.get_reference_name(pc.reference_id), pc.pos), 
             r.alignment.get_tag('CB'),
             r.alignment.query_name,
             r.alignment.mapping_quality,
             r.alignment.query_sequence[r.query_position:r.query_position + max(len(rec.ref), len(rec.alts[0])) + 1],
             rec.ref,
             rec.alts[0] #len(rec.alts) = 1 for every variant
            ) 
        
            for rec in variants 
            for pc in samfile.pileup(rec.chrom, rec.pos - 1, rec.pos, truncate=True) 
            for r in pc.pileups if not r.is_del and not r.is_refskip and r.alignment.has_tag('CB')
        
            if len(rec.alts) == 1 and len(rec.alts[0]) == 1
        ];
    return(pd.DataFrame.from_records(pileup_recs, columns=["region", "cell_barcode", "query_name", "mapping_quality", "query_base", "ref_base", "alt_base"]));



#####----------------------------end of methods, start of script-----------------------#####
if(len(sys.argv) != 4):
    print("Requires 3 arguments: vcf file path, bam file path, and output vcf file name")
    sys.exit()
else:
    vcf_file = sys.argv[1]
    bam_file = sys.argv[2]
    file_name = sys.argv[3]


## read in the original vcf file
with pysam.VariantFile(vcf_file) as variant_file:
    variants = [v for v in variant_file.fetch()]
    
filtered_variants = [v for v in variants if v.samples[0]['AO'][0] < 2]


## get the pileup from the bam file at the variant sites
pileup_at_variants = pileup(bam_file, filtered_variants);


## get all the information from the bam file that we need to make the new vcf file
var_stat = pileup_at_variants.query_base.str.slice(0, 1) == pileup_at_variants.alt_base
genotyped_pileup = pd.concat([pileup_at_variants['region'], pileup_at_variants['cell_barcode'], var_stat], axis=1, keys=['region', 'cell_barcode', 'var_stat'])
grouped_stats = genotyped_pileup.groupby(['cell_barcode', 'region'])
kept_cell_pileup = genotyped_pileup.groupby('cell_barcode').filter(lambda g: len(g['cell_barcode'])>3)

## calculate the stats we need to put into the vcf file
counts = kept_cell_pileup.groupby(kept_cell_pileup.columns.tolist(),as_index=False).size().to_frame()
counts_d = counts.reset_index()
counts_d.columns = ['region', 'cell_barcode', 'var_stat', 'count']

barcode_depth = counts_d.groupby(['region', 'cell_barcode'])['count'].agg('sum').to_frame()
barcode_depth_d = barcode_depth.reset_index()
barcode_depth_d.columns = ['region', 'cell_barcode', 'depth']

total_depth = counts_d.groupby('region')['count'].agg('sum').to_frame()
total_depth_d = total_depth.reset_index()
total_depth_d.columns = ['region', 'total_depth']

alternate_depth = counts_d[counts_d['var_stat'] == True].groupby('region')['count'].agg('sum').to_frame()
alternate_depth_d = alternate_depth.reset_index()
alternate_depth_d.columns = ['region', 'alternate_depth']

reference_depth = counts_d[counts_d['var_stat'] == False].groupby('region')['count'].agg('sum').to_frame()
reference_depth_d = reference_depth.reset_index()
reference_depth_d.columns = ['region', 'reference_depth']

barcode_depth_count = barcode_depth_d.merge(counts_d, on=['region', 'cell_barcode'], how='outer')
barcode_depth_total_depth_count = barcode_depth_count.merge(total_depth_d, on='region', how='outer')
with_alternate = barcode_depth_total_depth_count.merge(alternate_depth_d, on='region', how='outer')
with_ref_and_alt = with_alternate.merge(reference_depth_d, on='region', how='outer')


## filter out cells with low UMI count
barcodes, umi_counts = filter_by_umi("/uufs/chpc.utah.edu/common/home/ucgdlustre/work/u0991678/bowtell-singlecell/cellranger/14629X4/outs/possorted_genome_bam.bam", 15000, 0)
#barcodes.append(0)  ## append 0 to the list of possible barcodes because we fill in NA barcode values with 0
with_ref_and_alt = with_ref_and_alt.fillna(0)
with_ref_and_alt = with_ref_and_alt[with_ref_and_alt['cell_barcode'].isin(barcodes)]


## create variant data frame (with important information from the original vcf file) so we can merge it with the other data frame (that has bam file info)  
pileup_recs = [
    itertools.chain(("%s:%d" % (rec.chrom, rec.pos), 
     rec.id,
     rec.ref,
     rec.alts[0],
     rec.qual
    ),(rec.info[key][0] if type(rec.info[key]) == tuple else rec.info[key] for key in rec.info))

    for rec in variants
];
variants_df = pd.DataFrame.from_records(pileup_recs, columns=["region", "id", "ref_base", "alt_base", "quality"] +
 [key for key in variants[0].info]);

## merge variants_df and the dataframe we got from the bam file
variants_info = variants_df.merge(with_ref_and_alt, on='region', how='outer', sort=False)



#######------------------------------- write the vcf file ----------------------------------########
## write the vcf file
vcf_in = pysam.VariantFile(vcf_file, 'r', drop_samples=True)
vcf_in.subset_samples([])
vcf_out = pysam.VariantFile(file_name, 'w', header=vcf_in.header, drop_samples=True)


vcf_info = variants_info
vcf_info = vcf_info.fillna(0)

cell_barcodes = vcf_info['cell_barcode'].unique().tolist()
del cell_barcodes[cell_barcodes.index(0)]


## add each barcode as a sample
for barcode in cell_barcodes:
    vcf_out.header.add_sample(barcode)

## sort and group the vcf info dataframe
groups = vcf_info.groupby('region').groups
for key in groups:
    samples_me=[]
    for i in range(0, len(cell_barcodes)):
        samples_me.append({'GT':None, 'DP':None, 'AO':None, 'RO':None})
    for index in groups[key]:
        if(vcf_info.loc[index]['cell_barcode'] != 0):
            position = cell_barcodes.index(vcf_info.loc[index]['cell_barcode'])
            samples_me[position]['DP'] = int(vcf_info.loc[index]['depth'])
            if(vcf_info.loc[index]['var_stat'] == True):
                samples_me[position]['AO'] = int(vcf_info.loc[index]['count'])
                samples_me[position]['GT'] = vcf_info.loc[index]['alt_base']
            else:
                samples_me[position]['RO'] = int(vcf_info.loc[index]['count'])
                if(samples_me[position]['GT'] == None):
                    samples_me[position]['GT'] = vcf_info.loc[index]['ref_base']
    row = vcf_info.loc[groups[key][0]]

    ## add in the old info keys and their values
    info_old = {}
    for key in filtered_variants[0].info:
        try: # it doesn't like string values for some reason so we just... ignore those ones
            info_old[key] = int(row[key])
        except ValueError:
            continue
    # need to merge in this order so we change the values of DP, AO, and RO to reflect what we want
    info_all = merge_two_dicts(info_old, {'DP':int(row['total_depth']),
                                   'AO':int(row['alternate_depth']),
                                   'RO':int(row['reference_depth'])})
    entry = vcf_out.header.new_record(contig = row['region'].split(':')[0],
                            start = int(row['region'].split(':')[1])-1, 
                            stop = int(row['region'].split(':')[1]),
                            alleles = [row['ref_base'], row['alt_base']],
                            qual = row['quality'],
                            info = info_all,
                            samples=samples_me)
    vcf_out.write(entry)
    
vcf_out.close()

