
import argparse
import pysam
import collections
import sys

#python TAPE_10X_read2fromBAM.py \
#    --input_bam /outs/possorted_genome_bam.bam \
#    --output_file CellByTape_10X_test.csv \
#    --UMI_threshold 2


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Script to generate a file with counts of all TargetBC-TAPE barcodes found within each cell in a 10X experiment.')
    parser.add_argument('--input_bam', '-i', help='Position sorted BAM (or list of bams) from 10X pipestance.')
    parser.add_argument('--output_file', '-o', help='Tab delimited file with cell, mutation barcode, read count, umi count. All observed barcodes correctable to a whitelist are reported.')
    parser.add_argument('--UMI_threshold', type=float, default=2, help='Threshold for calling a UMI real.')

    args = parser.parse_args()

    input_bam=args.input_bam
    UMI_threshold=args.UMI_threshold
    output_file_name=args.output_file

    
    ############################################################
    # tallying up UMIs observed in cells for each mBC
    ############################################################
    
    mutation_barcodes = {}
    
    read_number=0
    read_mapped=0
    read_no_feature=0
    read_no_corrected_cBC_umi=0
    read_before_umiFilter=0
    read_final=0

    for read in pysam.Samfile(input_bam):
        read_number += 1

        if not read.is_unmapped:
            read_mapped += 1
            continue

        # get various part of the read in the bam file
        seq = read.seq.upper()
        tags = dict(read.tags)
        cell = tags.get('CB', None)
        umi = tags.get('UB', None)
        feature = tags.get('fb', None)

        # skip read if no corrected UMI or cBC from cell ranger
        if not cell or not umi:
            read_no_corrected_cBC_umi += 1
            continue

        if not feature:
            read_no_feature += 1

        # leave barcode error correction out for now, can implement once we have accepted lists.
        corrected_barcode=seq[29:41]+','+seq[58:64]+','+seq[78:84]+','+seq[98:104]+','+seq[118:124]+','+seq[138:144]+','+seq[158:164]

        # list of all barcodes already listed w/ cBC=cell
        barcodes_in_cell = mutation_barcodes.get(cell, dict())

        # add the corresponding UMI to cBC
        if corrected_barcode not in barcodes_in_cell:
            barcodes_in_cell[corrected_barcode] = []
        barcodes_in_cell[corrected_barcode].append(umi)

        # update full dictionary with updated list of mBC/umis
        mutation_barcodes[cell] = barcodes_in_cell
    
    print("read number:",read_number)
    print("mapped reads (not expected):",read_mapped)
    print("reads without corrected cBC+UMI:",read_no_corrected_cBC_umi)
    print("reads without feature flag:",read_no_feature)
    print("non mapped reads w/ corrected cBC+UMI:",read_number-read_mapped-read_no_corrected_cBC_umi-read_no_feature)
    #print("reads without exact match in GFP:", read_no_search_seq_match, read_no_search_seq_match/(read_number-read_mapped-read_no_corrected_cBC_umi))    
        

    ############################################################
    # generate output file 
    ############################################################
        
    output_file=open(output_file_name, 'w')
    original_stdout = sys.stdout
    output_file.write(','.join(['cBC','mBC','Site1','Site2','Site3','Site4','Site5','Site6','n_reads','n_UMI','\n']))
    #,'list_reads_per_UMIs','list_reads_per_UMIs_filtered'

    for cell in mutation_barcodes:

        # get the list of all UMI (each appearance corresponding to a read), irrespective of mBC, in a given cell (for a cBC)
        all_umis = []
        for mBC in mutation_barcodes[cell]:
            all_umis.extend(mutation_barcodes[cell][mBC])

        # tally up read counts for each UMI (again, irrespective of associated mBC)
        umi_counts_all_bc = collections.Counter(all_umis)    

        # loop through barcodes
        for mBC in mutation_barcodes[cell]:
            # counts of reads and UMIs for a given mBC, generate the list of read counts per UMI
            umi_counts_bc=collections.Counter(mutation_barcodes[cell][mBC])
            n_UMIs=len(umi_counts_bc.keys())
            n_reads=sum(umi_counts_bc.values())
            #list_umi_counts=list(umi_counts_bc.items())
            #str_list_umi_counts=[umi+"_"+str(count) for umi,count in list_umi_counts]

            # output to file: only print if non 0 filtered umis>0. 
            if n_UMIs >= UMI_threshold: 
                read_final += 1
                output_file.write(','.join([cell,mBC,str(n_reads),str(n_UMIs)]))
                output_file.write('\n')
                sys.stdout = output_file
                sys.stdout = original_stdout

    print("Final read number saved: ",read_final)


        
        
        
        
