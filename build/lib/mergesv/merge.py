import csv
import argparse
from argparse import RawTextHelpFormatter
from operator import itemgetter

__version__='0.0.1'
__usage__="""
                                                  
o     o                            .oPYo. o     o 
8b   d8                            8      8     8 
8`b d'8 .oPYo. oPYo. .oPYo. .oPYo. `Yooo. 8     8 
8 `o' 8 8oooo8 8  `' 8    8 8oooo8     `8 `b   d' 
8     8 8.     8     8    8 8.          8  `b d'  
8     8 `Yooo' 8     `YooP8 `Yooo' `YooP'   `8'   
<------------>            8     <----------->     
  <------------>       ooP'  <----------------->  
      <--------->      <-------->        <------->



Version:   {}
Author:    Dan Averbuj

About:     Merge SVs from a BED file based on reciprical overlap and distance to SV starts and ends
Usage:     mergesv <-i BED> [options]

Required Arguments:
	-i    FILE    BED file

Options:
	-o    FILE    output file [bed-file.mergesv.bed]
	-r    FLOAT   reciprocal overlap [0.8]
	-s    INT     distance (bp) from start and end values to be considered a merge [1000]
	-c            Include number of merged SVs in output

""".format(__version__)

#####################################
# FUNCTIONS   FUNCTIONS   FUNCTIONS #   
#####################################

def get_args():
	parser = argparse.ArgumentParser(usage=__usage__, add_help=False)
	parser.add_argument('-i'  ,type=str,   required=True,   default=None, help="Input BED file (chr, start, end, svtype)")
	parser.add_argument('-o'  ,type=str,   required=False,  default=None, help="Output file w/ directory")
	parser.add_argument('-r'  ,type=float, required=False,  default=0.8,  help="Reciprocal overlap of SVs to be merged [0.8]")
	parser.add_argument('-s'  ,type=int,   required=False,  default=1000, help="Limit of length between SV starts and ends [1000]")
	parser.add_argument('-c'  ,action="store_true", 					  help="Includes number of merged SVs in output")

	args = parser.parse_args()

	return args

def is_overlap(ref, sv):
	ref_start, ref_end = ref[1], ref[2]
	sv_start, sv_end = sv[1], sv[2]
	if sv_start > ref_end:
		return False
	else:
		return True

def is_recip_overlap(ref, sv):
	# checks for reciprocal overlap and within 1kb
	ref_start, ref_end = ref[1], ref[2]
	ref_size = int(ref_end) - int(ref_start)
	sv_start, sv_end = sv[1], sv[2]
	sv_size = int(sv_end) - int(sv_start)

	# sv overlap
	if ref_end >= sv_end: sv_overlap = 1.0
	else: sv_overlap = (int(ref_end)-int(sv_start))/sv_size

	# ref overlap
	if sv_end >= ref_end:
		ref_overlap = (int(ref_end)-int(sv_start))/ref_size
	else: ref_overlap = (int(sv_end)-int(sv_start))/ref_size

	# within 1kb start and end
	start_diff = abs(int(sv_start) - int(ref_start))
	end_diff = abs(int(sv_start) - int(ref_start))

	# check if x% reciprocal overlap
	if sv_overlap >= ro and ref_overlap >= ro and start_diff < edge_size and end_diff < edge_size:
		return True
	else:
		return False

def ci_positions(val):
	# returns ci positions globally, not relative to the SVs start and end
	start, end, cipos, ciend = val[1], val[2], val[4], val[5]
	cipos_s, cipos_e = val[4].split(',')
	ciend_s, ciend_e = val[5].split(',')
	cpsp, cpep = (start + int(cipos_s)), (start + int(cipos_e))
	cesp, ceep = (end + int(ciend_s)), (end + int(ciend_e))
	return cpsp, cpep, cesp, ceep

def merge_and_delete(merged_list, svlist, values):
	## determine consensus merge values and add to list

	# used for first iteration of merge where no cipos, ciend, or n_merge exist
	if len(values[0]) == 4:
		n_merge = len(values)
		if n_merge == 1:
			cipos, ciend = '.', '.'
			chrom, start, end, svtype = values[0][0], values[0][1], values[0][2], values[0][3]
		else:
			chrom, svtype = values[0][0], values[0][3]
			starts = []
			ends = []
			for val in values:
				starts.append(val[1])
				ends.append(val[2])
			starts.sort()
			ends.sort()
			if len(starts) % 2 == 1 and len(ends) % 2 == 1:
				start = starts[int((len(starts)-1)/2)]
				end = ends[int((len(ends)-1)/2)]
			elif len(starts) % 2 == 0 and len(ends) % 2 == 0:
				start = starts[int(((len(starts)-1)/2)-.5)]
				end = ends[int(((len(ends)-1)/2)-.5)]
			else: print('MISTAKE - line 147')

			cipos_start = int(int(starts[0])-int(start))
			cipos_end = int(int(starts[-1])-int(start))
			ciend_start = int(int(ends[0])-int(end))
			ciend_end = int(int(ends[-1])-int(end))
			cipos = '{},{}'.format(cipos_start, cipos_end)
			ciend = '{},{}'.format(ciend_start, ciend_end)

		if args.c:
			merged_list.append([chrom, start, end, svtype, cipos, ciend, n_merge])
		else:
			merged_list.append([chrom, start, end, svtype, cipos, ciend])

	# used after the first merge loop, when args.c (-c) was not flagged
	elif len(values[0]) == 6:
		if len(values) == 1:
			merged_list.append( values[0] )
		else:
			chrom, svtype = values[0][0], values[0][3]
			starts = []
			ends = []
			cipos_start, cipos_end = [], []
			ciend_start, ciend_end = [], []
			for val in values:
				starts.append(val[1])
				ends.append(val[2])
				if val[4] == '.' and val[5] == '.':
					continue
				cps, cpe, ces, cee = ci_positions(val)
				cipos_start.append(cps)
				cipos_end.append(cpe)
				ciend_start.append(ces)
				ciend_end.append(cee)				
			starts.sort()
			ends.sort()
			if len(starts) % 2 == 1 and len(ends) % 2 == 1:
				start = starts[int((len(starts)-1)/2)]
				end = ends[int((len(ends)-1)/2)]
			elif len(starts) % 2 == 0 and len(ends) % 2 == 0:
				start = starts[int(((len(starts)-1)/2)-.5)]
				end = ends[int(((len(ends)-1)/2)-.5)]
			else: print('MISTAKE - line 161')
			# gets the smallest _start positions and the largest _end positions and determine cistart
			# and ciend based on the chosen start and end values
			ps = min(cipos_start) - start
			pe = max(cipos_end) - start
			es = min(ciend_start) - end
			ee = max(ciend_end) - end
			cipos = '{},{}'.format(ps, pe)
			ciend = '{},{}'.format(es, ee)

			merged_list.append( [chrom, start, end, svtype, cipos, ciend] )

			
	# used after the first merge loop, when args.c (-c) was flagged
	elif len(values[0]) == 7:
		if len(values) == 1:
			merged_list.append( values[0] )
		else:
			chrom, svtype = values[0][0], values[0][3]
			# determine the new n_merge, the start and end value of the most previously-merged SV,
			# and compile the global cipos and ciend positions
			n_merge = 0
			start_end = (values[0][1], values[0][2], values[0][6])
			cipos_start, cipos_end = [], []
			ciend_start, ciend_end = [], []
			for val in values:
				n_merge += val[6]
				if val[6] > start_end[2]:
					start_end = (val[1], val[2], val[6])
				if val[4] == '.' and val[5] == '.':
					continue
				cps, cpe, ces, cee = ci_positions(val)
				cipos_start.append(cps)
				cipos_end.append(cpe)
				ciend_start.append(ces)
				ciend_end.append(cee)
			start = start_end[0]
			end = start_end[1]
			# gets the smallest _start positions and the largest _end positions and determine cistart
			# and ciend based on the chosen start and end values
			ps = min(cipos_start) - start
			pe = max(cipos_end) - start
			es = min(ciend_start) - end
			ee = max(ciend_end) - end
			cipos = '{},{}'.format(ps, pe)
			ciend = '{},{}'.format(es, ee)

			merged_list.append( [chrom, start, end, svtype, cipos, ciend, n_merge] )

	else:
		print("Length of values to merge are wrong")
		
	# delete values from svlist
	for v in values:
		while v in svlist:
			svlist.remove(v)
	return

def merge_list(sv_list):
	## Loop through sv list to actively merge SVs
	merged_list = []
	while len(sv_list) > 0:
		ilen = len(sv_list)
		i = 0
		ref_allocated = False
		to_merge = []

		## Loop through SV list and compile SVs to be merged. At the end, write the merged SV from those to file and delete from SV list
		for sv in sv_list:
			i += 1
			# initialize reference SV (first in list)
			if not ref_allocated:
				ref = sv
				to_merge.append(ref)
				ref_allocated = True
				# if looking at last SV in the list
				if i == ilen: # i == ilen when last SV in the list
					merge_and_delete(merged_list, sv_list, to_merge)
					break
				continue
			## compare following SVs until no overlap to ref is found
			# if there is no overlap to ref, write merged SV to file and delete the merging SVs from the list of all SVs
			elif not is_overlap(ref, sv):
				merge_and_delete(merged_list, sv_list, to_merge)
				break
			# if there is 80% reciprocal overlap to the ref SV and start and end positions are within 1kb, add to list of SVs to be merged
			elif is_recip_overlap(ref, sv):
				to_merge.append(sv)
				if i == ilen:
					merge_and_delete(merged_list, sv_list, to_merge)
					break
			elif i == ilen:
				merge_and_delete(merged_list, sv_list, to_merge)
				break
			continue
	return merged_list

def write_to_file(owriter, merged_list):
	for sv in merged_list:
		owriter.writerow( sv )


####################################
# MAIN   MAIN   MAIN   MAIN   MAIN #   
####################################

def main()

	args = get_args()

	bed = args.i
	ro = args.r
	edge_size = args.s
	if not args.o:
		if not bed.endswith('.bed'):
			sys.stderr.write("Warning: input file does not end in '.bed'. Adding default suffix ('.bed') to end of file name")
			outfile = bed + '.mergesv.bed'
		else:
			outfile = bed[:bed.rfind('.bed')] + '.mergesv.bed'
	else:
		outfile = args.o

	bed_fh = open(bed)
	out_fh = open(outfile, 'w', newline='')
	out_writer = csv.writer(out_fh, delimiter='\t')

	## Create dictionary separating by chromosome and svtype
	split_dict = {}
	for sv in csv.reader(bed_fh, dialect='excel-tab'):
		chrom, svtype = sv[0], sv[3]
		if not chrom in split_dict:
			split_dict[chrom] = {'DEL':[], 'DUP':[], 'INV':[]}
		split_dict[chrom][svtype].append([sv[0], int(sv[1]), int(sv[2]), sv[3]])

	## Loop through each chromosome and SV type list, merging all SVs to one outfile
	for chrom in split_dict:
		for svtype in split_dict[chrom]:
			print(chrom, svtype)
			# Sort SVs by start position
			split_dict[chrom][svtype].sort(key=itemgetter(1))
			sv_list = split_dict[chrom][svtype]

			# recursive loop to merge SVs, ensuring to merge previously-merged SVs
			while True:
				premerge_length = len(tuple(sv_list))
				merged_list = merge_list(sv_list)
				# if the SVs did not merge any further, break out of loop and write output to file
				if len(merged_list) == premerge_length:
					write_to_file(out_writer, merged_list)
					break
				else:
					print("Size of old list: {}\nSize of merged list: {}".format(premerge_length, len(merged_list)))
					sv_list = merged_list

	out_fh.close()