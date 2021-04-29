############################################################
# allele_counts.py
############################################################

import sys
import os
import pysam
import gzip

min_mapq = 30

#####################################################################
# Parse input arguments
#####################################################################

usage_string = "Usage: "+sys.argv[0]+" <sample_bam> <output_file> <somatic_vcf> [<germline_vcf>]"


try:
    bam_input = os.path.abspath(sys.argv[1])
    out_file =  os.path.abspath(sys.argv[2])
    vcf_somatic = os.path.abspath(sys.argv[3])
 
    if len(sys.argv) > 4:
        vcf_germline = os.path.abspath(sys.argv[4])
    else:
        vcf_germline = None
except:
    raise SystemExit(usage_string)



if not os.access(bam_input, os.R_OK):
    raise SystemExit("Error: BAM input file '"+bam_input+\
                     "' does not exist or unreadable.\n\n"+usage_string)

if not os.access(vcf_somatic, os.R_OK):
    raise SystemExit("Error: Somatic VCF file '"+vcf_somatic+\
                     "' does not exist or unreadable.\n\n"+usage_string)

if not os.access(vcf_germline, os.R_OK) and vcf_germline is not None:
    raise SystemExit("Error: Germline VCF file '"+vcf_germline+\
                     "' does not exist or unreadable.\n\n"+usage_string)

opath = os.path.dirname(out_file)
if not os.path.exists(opath):
    opath = os.path.dirname(out_file)
    os.mkdir(opath) 
    print("=> Created output path '"+opath+"'\n")
    

#####################################################################
# Functions
#####################################################################


def skip_vcf_header(f):
    line = ''
    while not line.startswith('#CHROM'):
        line = f.readline()
    return(f)


def get_base_counts(sf, chr, pos, ref, alt):
    ref_count = 0
    alt_count = 0
    depth = 0

    fetched_qnames = list() # don't count same base from read pairs twice.

    for read in sf.fetch(reference=chr, start=pos-1, end=pos-1+len(alt)):

        if read.is_secondary or read.is_supplementary or read.is_duplicate or \
           read.mapping_quality < min_mapq or not read.is_proper_pair or \
           read.qname in fetched_qnames:
            continue

        rpos = read.reference_start + read.query_alignment_start
        idx_start = pos - 1 - rpos
        if idx_start >= 0 and idx_start+len(alt) <= read.query_alignment_length:
            qseq = read.query_alignment_sequence[idx_start:idx_start+len(alt)]
            if qseq == ref:
                ref_count += 1
            if qseq == alt:
                alt_count += 1
            depth += 1
            fetched_qnames.append(read.qname)

    return([ref_count, alt_count, depth])


def parse_vcf_line_germline(line, samfile):
    fields = line[:-1].split('\t')
    reference_data = dict(zip(fields[-2].split(":"), fields[-1].split(":")))
    reference_gt = reference_data['GT']
    return(parse_vcf_line(line, samfile, GT=reference_gt))


def parse_vcf_line(line, samfile, GT=None):

    # split the line into the variant annotations:
    (chr,pos,annotation,ref,alts) = line[:-1].split('\t')[0:5] # chr, position, annotation, ref, alts
    pos = int(pos)

    if ',' in alts:
        alts = alts.split(',')
    else:
        alts = [alts]


    # parse gt argument:
    if GT is not None:
        gt_idx = [int(e) for e in GT.split('/')]

        for i in gt_idx:
            if i > len(alts) or i < 0:
                raise ValueError("Invalid GT field.")

        if len(gt_idx) > 2:
            raise ValueError("Invalid GT field.")

        if gt_idx[0] == gt_idx[1]:
            return(list())               # Ignore homozygous sites.

        alleles = [ref] + alts
        ref = alleles[gt_idx[0]]
        alts = [alleles[gt_idx[1]]]


    # get counts for all valid alleles:
    return_value = list()
    for alt in alts:

        if len(ref) != len(alt): # InDels are skiped
            continue

        counts = get_base_counts(samfile, chr, pos, ref, alt) # ref, alt, depth
        return_value.append([chr,pos,ref,alt,annotation]+counts)

    return(return_value)


def concat_elements(elements, delim='\t'):
    line = ''
    for element in elements:
        line += str(element) + delim
    line = line[:-1] + '\n'
    return(line)



#####################################################################
# Main
#####################################################################


with pysam.Samfile(bam_input, 'rb') as samfile:
    with gzip.open(out_file,'wt') as output:

        header_elements = ["sample","chr","pos","ref","alt","annotation","ref_count","alt_count","depth","source"]
        output.write(concat_elements(header_elements))

        sample = os.path.basename(bam_input).replace(".bam", "") # get sample name from input file name.

        with skip_vcf_header(gzip.open(vcf_somatic, 'rt')) as somatic_calls:
            for site in somatic_calls:
                for return_line in parse_vcf_line(site, samfile):
                    out_line = concat_elements([sample]+return_line+["S"])
                    output.write(out_line)
        
        if vcf_germline is not None:
            with skip_vcf_header(gzip.open(vcf_germline, 'rt')) as germline_calls:
                for line in germline_calls:
                    for return_line in parse_vcf_line_germline(line, samfile):
                        out_line = concat_elements([sample]+return_line+["GL"])
                        output.write(out_line)

