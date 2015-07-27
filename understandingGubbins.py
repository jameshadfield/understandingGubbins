#!/usr/bin/env python

import pdb
import argparse
import re,sys,os,string
from pprint import pprint
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_j/jh22/git/understandingGubbins']))
from gubbins_parser import Rec,MultipleRec,RecAln,RecGenic,EMBLparser

def get_user_options():
	parser = argparse.ArgumentParser(description="making sense of gubbins' rec file")
	parser.add_argument('-r', required=True, action='store', dest="tabfile", help="Output from gubbins. In the old versions this was *_rec.tab, in the newer versions its *.embl.", default="", metavar="FILE")
	parser.add_argument("-e", action="store", dest="emblfile", help="emblfile", default="", metavar="FILE")
	parser.add_argument("--taxa", action="store", dest="taxafile", help="analyse only these taxa", default="", metavar="FILE")
	parser.add_argument("--invert", action="store_true", dest="inverttaxafile", help="exclude the taxa (--taxa) rather than include", default=False)


	table = parser.add_argument_group('create tables for that nature paper of yours')
	table.add_argument("--table", action="store", dest="tablefile", help="outfile for table creation", default="", metavar="FILE")

	stats = parser.add_argument_group('s t a t i s t i c s')
	stats.add_argument("--blocksize", action="store_true", dest="blocksize", help="report length of recombination fragments", default=False)
	# stats.add_argument("--perc_recomb", action="store_true", dest="perc_recomb", help="report percentage of the genome involved in recombination events", default=False)

	hotspot = parser.add_argument_group('create hotspot plots for iCANDY')
	hotspot.add_argument("--hotspot", action="store", dest="hotspotfile", help="outfile for hotspot plot file (needs genome length)", default="", metavar="FILE")
	hotspot.add_argument("--missingplot", action="store", dest="missingdataplotfile", help="outfile for missing data plot file (needs genome length and num taxa)", default="", metavar="FILE")
	hotspot.add_argument("--inversemissing", action="store_true", dest="inversemissing", help="missing data inversed, i.e. shows proportion of data", default=False)
	hotspot.add_argument("-l", action="store", dest="genome_length", help="genome_length", default="", metavar="INT")
	hotspot.add_argument("-n", action="store", dest="num_taxa", help="number of taxa", default="", metavar="INT")

	merge = parser.add_argument_group('merge the tabfile blocks')
	merge.add_argument("-m", action="store", dest="merge", help="Merge the tabfiles in one way or another (this file is the merged output)", default="", metavar="FILE")
	merge.add_argument("--mergemethod", action="store", dest="merge_method", help="block merge method to use ['collapsed' (default) or 'heat' or 'overlap' or 'compare']. Compare only works if you provide multiple rec files (-r 1.tab,2.tab,...)", default="collapsed", metavar="STRING")
	merge.add_argument("--brewer", action="store", dest="brewer", help="if using heatmap, name your ColourBrewer colour scheme [default: Spectral]", default="Spectral", metavar="STR")
	merge.add_argument("--colours", action="store", dest="colours", help="if using heatmap this is the maximum number of colours to be used. Each colour is one level with levels >= colours -> final coulour. If using 'collapsed' then this is the colour code (Artemis) for them all. For 'overlap' this should be 2 numbers comma-seperated.", default="", metavar="INT/STR")

	beast = parser.add_argument_group('remove recombinations from an alignment for BEAST')
	beast.add_argument("-a", "--aln", action="store", dest="beastalnin", help="multi-fasta alignment file", default="", metavar="FILE")
	beast.add_argument("-b", "--beast", action="store", dest="beastprefix", help="prefix for BEAST alignment outputs (recombination regions removed)", default="", metavar="FILE")
	beast.add_argument("-d", "--dates", action="store", dest="dates", help="Dates csv file name (Name,date)", default="", metavar="FILE")
	beast.add_argument("-t", "--tree", action="store", dest="treefile", help="Gubbins treefile to create starting tree for BEAUTi (optional)", default="", metavar="FILE")


	genes = parser.add_argument_group('score genes based on recombination events')
	genes.add_argument("-g", action="store", dest="genes", help="file to save recombination events per genes (normalised by default)", default="", metavar="FILE")
	genes.add_argument("--genepos", action="store", dest="genepos", help="file to save gene positions", default="", metavar="FILE")

	return parser.parse_args()



# def basic_EMBL_parsing(emblfile):
# 	print "Parsing EMBL(-like) file"
# 	lines = open(emblfile).read().splitlines()
# 	lines = [x for x in lines if x.startswith("FT")] ## filtering
# 	geneDB = []
# 	for block in generate_blocks_from_lines(lines,['gene','CDS']):
# 		h = {'name':None,'x1':None,'x2':None,'strand':'+','locus_tag':None}
# 		m = re.search(r"(\d+)\.\.(\d+)",block[0])
# 		h['x1']=int(m.group(1))
# 		h['x2']=int(m.group(2))
# 		if h['x1'] > h['x2']:
# 			h['x1'],h['x2'] = h['x2'],h['x1']
# 		if 'complement' in block[0]:
# 			h['strand'] = '-'
# 		## rest of the lines...
# 		for line in block: ## won't match line 1 anyways
# 			if '/gene' in line:
# 				h['name'] = line.split('"')[1]
# 			if '/locus_tag' in line:
# 				h['locus_tag'] = line.split('"')[1]
# 		geneDB.append(h)
# 	return geneDB


# def generate_blocks_from_lines(text,matches): ## generator
# 	"""the string "misc_feature" seperates the block. returns (yields) list of lines."""
# 	buff = []
# 	for line in text:
# 		feature = line.split()[1] ## e.g. CDS,gene, but also /product=... /gene=...
# 		if feature in matches and buff!=[]:
# 			if buff:
# 				yield buff
# 				buff = [line]
# 		else:
# 			buff.append(line)
# 	if buff:
# 		yield buff

# def draw_table_of_blocks_and_taxa_and_genes(options):
# 	print "Constructing table of recombination regions"
# 	if options.emblfile:
# 		geneDB = basic_EMBL_parsing(options.emblfile)
# 		rec.append_genes_into_regions(geneDB)
# 	## open fh for writing:
# 	if options.tablefile is not "stdout":
# 		fh=open(options.tablefile,'w')
# 	else:
# 		fh=sys.stdout
# 	fh.write("REGION\tCO-ORDS\tTAXA\tSNPs\tGENES\n")
# 	for idx,region in enumerate(rec.regions):
# 		try:
# 			genes = region['genes']
# 		except KeyError:
# 			genes = ["N/A"]
# 		fh.write("{}\t{}..{}\t{}\t{}\t{}\n".format(idx+1,region['x1'],region['x2'],",".join(region['taxa']),region['SNPs'],",".join(genes)))
# 	if fh is not sys.stdout:
# 	    fh.close()



if __name__ == "__main__":
	options = get_user_options()
	rec = Rec(options.tabfile.split(',')[0])
	if options.taxafile:
		# b4 = (len(rec.regions),len(rec.regions_by_taxa_unmerged))
		rec.subset_taxa(options.taxafile,options.inverttaxafile)
		print "Subsetting on {}.\n\treduced the number or recombination regions {}->{} and the number of taxa (with recombinations) from {}->{}".format(options.taxafile,len(rec.regions),len(rec.regions),len(rec.regions_by_taxa_unmerged),len(rec.regions_by_taxa_unmerged))
	if options.tablefile:
		print "no longer working :("
		# draw_table_of_blocks_and_taxa_and_genes(options)
	if options.hotspotfile:
		rec.hotspot_plot(options.hotspotfile,options.genome_length)
	if options.missingdataplotfile:
		rec.missing_data_plot(options.missingdataplotfile,options.genome_length,options.num_taxa,options.inversemissing)
	if options.emblfile:
		embl = EMBLparser(options.emblfile)

	# ------ stats and stuff -------- #
	if options.blocksize:
		rec.report_block_stats(rec.regions)


	# ------- mergeing stuff -------- #
	if options.merge:
		if options.merge_method in ['collapsed','heat','overlap']:
			rec.merge_regions(options.merge_method,colourscheme=options.brewer,colours=options.colours)
			rec.write_regions_merged(options.merge)
		elif options.merge_method == 'compare':
			tabfiles = options.tabfile.split(',')
			assert(len(tabfiles)>1)
			recs = [ Rec(tabfile) for tabfile in tabfiles ]
			multiplerec = MultipleRec(*recs) ## MultipleRec.__init__() willcollapse the individual regions
			multiplerec.merge_tabfiles()
			multiplerec.write_regions_merged(options.merge)


	# ------- BEAST (alignment) stuff -------- #
	if options.beastalnin and options.beastprefix:# and options.dates:
		beast = RecAln(rec,options.beastalnin)
		beast.set_alignment_for_BEAST()
		beast.write_full_alignment(options.beastprefix)
		if options.dates:
			beast.prepare_BEAST_alignment(options.dates,options.beastprefix,options.treefile)
		else:
			beast.prepare_BEAST_alignment(None,options.beastprefix,options.treefile)


	# ------- genes most / more frequently observed as recombining -------- #
	if options.genes:
		## rec and embl objects exist already
		## if taxa subsetting is selected this, too, has already been done
		recgene = RecGenic(rec,embl)
		recgene.normalise_recombination_counts(by='numtaxa')
		recgene.write_recombination_counts(options.genes)
	if options.genepos:
		embl.write_gene_positions(options.genepos)























