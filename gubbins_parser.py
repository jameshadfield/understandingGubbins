import pdb
from pprint import pprint
import sys,re,os
from math import ceil
sys.path.append('/nfs/users/nfs_j/jh22/usr/lib/python2.7/site-packages/brewer2mpl')
import brewer2mpl
from Bio import SeqIO
from subprocess import check_call,CalledProcessError
import numpy as np
import scipy.stats


class RecBase(object):
	""" This should never be called directly.
	It contains methods that Rec and MultipleRec inherit
	"""
	def __init__(self):
		pass ## no init

	def _create_breakpoint_hash(self,DB):
		HSplus,HSminus={},{}
		for region in DB:
			try:
				HSplus[region['x1']] += 1
			except KeyError:
				HSplus[region['x1']] = 1
			try:
				HSminus[region['x2']] += 1
			except KeyError:
				HSminus[region['x2']] = 1
		return HSplus,HSminus

	def hotspot_plot(self,outfile,genome_length):
		""" IN: self.regions
		rec=Rec() calls this and rec.regions exists but not RecBase.regions. No Problem.
		"""
		## two hashs to speed up the hotspot plot creation
		HSplus,HSminus = self._create_breakpoint_hash(self.regions)
		count=0
		with open(outfile,'w') as fh:
			fh.write("#BASE\tR_COUNT\n")
			for pos in xrange(1,int(genome_length)):
				if pos in HSplus:
					count += HSplus[pos]
				if pos in HSminus:
					count -= HSminus[pos]
				fh.write("{}\t{}\n".format(pos,count))

	def write_regions_merged(self,filename): #,db=None):
		""" writes self.regions_merged  to a tabfile """
		# if not db:
		db=self.regions_merged
		with open(filename,'w') as fh:
			for taxa in db:
				for x in db[taxa]:
					fh.write(  self._print_tab_line("{}..{}".format(x['x1'],x['x2']),"misc_feature")    )
					### is the colour an RGB????
					if 'colour' in x:
						fh.write(  self._print_tab_line("/colour={}".format(x['colour']))               )
					elif 'RGB' in x:
						fh.write(  self._print_tab_line("/colour=\"{} {} {}\"".format( x['RGB'][0],x['RGB'][1],x['RGB'][2] ) ) )
					fh.write(  self._print_tab_line("/taxa=\"{}\"".format(taxa))                        )


	def _print_tab_line(self,contents,field=False):
		if field:
			pad = 16 - len(field)
			return "FT   "+field+" "*pad+contents+"\n"
		else:
			return "FT"+" "*19+contents+"\n"

	def report_block_stats(self, regions):

		print "Recombination blocks identified: {}".format(len(regions))

		lens = [ abs(x['x2'] - x['x1']) for x in regions  ]
		meanLen = int(np.mean(lens))
		ci = scipy.stats.norm.interval(0.95, loc=meanLen, scale=scipy.stats.sem(lens))
		#return (meanLen,int(ci[0]),int(ci[1]))
		print "Mean block size: {}bp".format(meanLen)
		print "95% CI: {}bp - {}bp".format(int(ci[0]),int(ci[1]))






























class Rec(RecBase):
	def __init__(self,tabfile):
		super(Rec, self).__init__() ## even though this does nothing
		self.regions,self.missing_data_regions=self._parse_tabfile(tabfile)
		## self.missing_data_regions will be None unless gubbins was run with the -G option
		self.regions_by_taxa_unmerged=self._associate_regions_with_taxa()
		self.tabfile=str(tabfile) ## use os.path.basename here e.t.c.
		self.taxa = set()
		for r in self.regions:
			for t in r['taxa']:
				self.taxa.add(t)


	def subset_taxa(self,taxalistfile,invert):
		taxalist = []
		with open(taxalistfile,'rU') as fh:
			for line in fh:
				if line.startswith('#'):
					continue
				taxalist.append(line.strip())
		print "Num in taxalist: {}".format(len(taxalist))
		oldregions = self.regions
		self.regions = []
		for region in oldregions:
			if invert:
				newtaxa = [x for x in region['taxa'] if x not in taxalist]
			else:
				newtaxa = [x for x in region['taxa'] if x in taxalist]
			if newtaxa:
				region['taxa'] = newtaxa
				self.regions.append(region)

		self.regions_by_taxa_unmerged=self._associate_regions_with_taxa()


	def _parse_tabfile(self,tabfile):
		""" called during class initialisation. returns """
		lines = open(tabfile).read().splitlines()
		regions = []
		MDregions = []
		for block in self._generate_blocks_from_lines(lines):
			## new version of gubbins has "missing data" as blocks. must filter these out!
			if filter(lambda x:'missing data' in x, block):
				MDregions.append(self._parse_missing_data_block(block))
			else: ## normal block with recombination data
				regions.append(self._parse_block(block))
		if MDregions:
			return regions,MDregions
		else:
			return regions,None

	def _generate_blocks_from_lines(self,text): ## generator
		"""the string "misc_feature" seperates the block. returns (yields) list of lines."""
		buff = []
		for line in text:
			if 'misc_feature' in line and buff!=[]:
				if buff:
					yield buff
					buff = [line]
			else:
				buff.append(line)
		if buff:
			yield buff

	def _parse_block(self,block):
		"""turns list of lines (i.e. block) into hash of salient features"""
		h = {'colour':None,'x1':None,'x2':None,'taxa':None,'SNPs':None,'SNPratio':None,'log':None,'pvalue':None,'nodeFrom':None,'nodeTo':None,'numTaxa':None}
		for line in block:
			if "misc_feature" in line:
				h['x1']=int(line.split()[2].split('..')[0])
				h['x2']=int(line.split('..')[1].split()[0])
				if h['x1'] > h['x2']:
					h['x1'],h['x2'] = h['x2'],h['x1']
			elif "colour" in line:
				try:
					h['colour'] = int(line.split('=')[1])
				except ValueError:
					h['colour'] = 1
			elif "node" in line:
				if '"' in line:
					h['nodeFrom'] = line.split('"')[1].split('->')[0]
					h['nodeTo'] = line.split('->')[1].split('"')[0]
				else:
					pass
			elif "SNP_count" in line:
				h['SNPs'] = int(line.split('=')[1].strip('"'))
			elif "neg_log_likelihood" in line:
				h['log'] = float(line.split('=')[1].strip('"'))
			elif "ecombination_to_background_SNP_ratio" in line:
				h['SNPratio'] = float(line.split('=')[1])
			elif "pvalue" in line:
				h['log'] = float(line.split('=')[1])
			elif "taxa" in line:
				taxa_bit = line.split('"')[1]
				commas = taxa_bit.count(',')
				if commas: ## simons version
					h['taxa'] = [x.strip() for x in taxa_bit.split(',')]
				else: ## run_gubbins
					h['taxa'] = [x.strip() for x in taxa_bit.split()]
				h['numTaxa'] = len(h['taxa'])
		## sanity checking
		if not (h['x1'] and h['x2'] and h['numTaxa']):
			print "WARNING: POORLY FORMED BLOCK PARSED..."
			pprint(block)
			pprint(h)
		# print "TAXA:"
		# pprint(h['taxa'])
		return h

	def _parse_missing_data_block(self,block):
		"""turns list of lines (i.e. block) into hash of salient features"""
		h = {'colour':None,'x1':None,'x2':None,'taxa':None}
		for line in block:
			if "misc_feature" in line:
				h['x1']=int(line.split()[2].split('..')[0])
				h['x2']=int(line.split('..')[1].split()[0]) + 1 	#### hack. simon to fix in gubbins
				if h['x1'] > h['x2']:
					h['x1'],h['x2'] = h['x2'],h['x1']
			elif "colour" in line:
				h['colour'] = int(line.split('=')[1])
			elif "taxa" in line:
				h['taxa'] = [x.strip().strip('"') for x in line.split('=')[1].split(',')]
		## sanity checking
		if not (h['x1'] and h['x2'] and h['taxa']):
			print "WARNING: POORLY FORMED MISSING DATA BLOCK PARSED..."
			pprint(block)
			pprint(h)
		return h


	def _associate_regions_with_taxa(self):
		""" regions (AoH) -> regions_by_taxa (H (keys: taxa) of A of H (keys: x1,x2,colour)) """
		regions_by_taxa = {}
		for region in self.regions:
			for taxa in region['taxa']:
				try:
					regions_by_taxa[taxa].append({"x1":region['x1'],"x2":region['x2'],"colour":region['colour']})
				except KeyError:
					regions_by_taxa[taxa]=[]
					regions_by_taxa[taxa].append({"x1":region['x1'],"x2":region['x2'],"colour":region['colour']})
		return regions_by_taxa



	# D E P R E C I A T E D
	def append_genes_into_regions(self,geneDB,key_name='name',backup_key_name='locus_tag'):
		"""geneDB must be a list of dicts with keys x1,x2,and name (create this yourself)"""
		genes = sorted(geneDB,key=lambda k:k['x2'])
		for idx,region in enumerate(self.regions):
			self.regions[idx]['genes']=[]
			x1,x2 = region['x1'],region['x2']
			for gene in genes:
				if ( gene['x1'] > region['x1'] and gene['x1'] < region['x2'] ) or ( gene['x2'] > region['x1'] and gene['x2'] < region['x2'] ):
					name="unknown gene"
					if gene[key_name]:
						name = gene[key_name]
					elif gene[backup_key_name]:
						name = gene[backup_key_name]
					self.regions[idx]['genes'].append(name)
				if gene['x2'] > region['x2']:
					break # speed
		return True # success


	def merge_regions(self,method='collapsed',singlecolour=2,colourscheme='Spectral',colours=False):
		"""
		This takes the regions_by_taxa and returns a regions_merged
		Each taxa is treated independently
		The algorithm is basically:
			(for each taxa) breakpoints stored in hashes
			new blocks created from breakpoints

		We have a lot of flexibility in how we create (read display) the new breakpoints
			collapsed: presence | absense of region
			heat: the colour is a function of how many overlaps there are
			overlap: two colours are used, 1 for singleton, 1 for overlap (of any depth)
			fragment: colours are maintained and and overlap results in a block being fragmented. very similar to heat but with subtle differences (e.g. adjoing, colours e.t.c.)
		The colour parsing is a bit difficult. Basically we discard it at the moment.

		Defaults: collapsed with colour red

		"""
		self.regions_merged_method = method ## for later examination, if needed
		self.regions_merged = {taxa:[] for taxa in self.regions_by_taxa_unmerged}
		print "merging tabfile {} (method={})".format(self.tabfile,method)
		for taxa in self.regions_merged:
			# print " -------- taxa: "+taxa+" ------------- "
			rOpen,rClose = self._create_breakpoint_hash(self.regions_by_taxa_unmerged[taxa])
			if method=='collapsed':
				if not colours:
					colours = 2 ## default = RED
				for region in self._generate_collapsed_regions(rOpen,rClose,max(rClose.keys())+1,int(colours)):
					self.regions_merged[taxa].append(region)
			if method=='heat':
				for region in self._generate_heat_regions(rOpen,rClose,max(rClose.keys())+1):
					self.regions_merged[taxa].append(region)
					# colorised seperately... see below!
			if method=='overlap':
				overlapCols = [ int(x) for x in colours.split(',')] if colours else [5,6]
				for region in self._generate_overlap_regions(rOpen,rClose,max(rClose.keys())+1,overlapCols):
					self.regions_merged[taxa].append(region)

		if method=='heat':
			if not colours:
				colours=11
			self._colourise_heatmap(colourscheme,int(colours))

	def _colourise_heatmap(self,colourscheme,maxColours):
		""" by this stage we already have a regions_merged structure with the depth. This method turns those depths into RGB calls using colorbrewer
		"""
		## what's the maximum depth?
		maxRec = max( [ max( [int(y['depth']) for y in self.regions_merged[taxa] ] ) for taxa in self.regions_merged ] )
		## what scheme are we employing?
		if colourscheme in ['BrBG','PRGn','PiYG','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral']:
			category="Diverging"
			maxCol = 11
			reverse = True
		elif colourscheme in ['Blues','BuGn','BuPu','GnBu','Greens','Greys','OrRd','Oranges','PuBu','PuBuGn','PuRd','Purples','RdPu','Reds','YlGn','YlGnBu','YlOrBr','YlOrRd']:
			category="Sequential"
			maxCol = 9
			reverse = False
		else:
			print "Invalid ColorBrewer scheme. Tabfile will not have any colours."
			return 1

		if maxColours < maxCol: maxCol = maxColours
		if maxRec < 3: maxRec = 3 ## brewer can't handle only 2

		# create depth -> RGB maps
		depth2RGB = {}
		if maxRec <= maxCol: ## easy! 1-1 map
			print "Colorising heatmap using the Colour Brewer scheme {} with {} colours for {} levels".format(colourscheme,maxRec,maxRec)
			brewerMap = brewer2mpl.get_map(colourscheme,category,maxRec,reverse=reverse).colors
			for x in xrange(1,maxRec+1):
				depth2RGB[x] = brewerMap[x-1]
		else:
			print "Colorising heatmap using the Colour Brewer scheme {} with {} colours for {} levels by binning all recombinations >= {} into the top colour".format(colourscheme,maxCol,maxRec,maxCol)
			## we will just use the max num of categories here!
			levels_per_colour = int(ceil(float(maxRec)/maxCol))
			brewerMap = brewer2mpl.get_map(colourscheme,category,maxCol,reverse=reverse).colors
			# depth,brewerLevel = 0,0
			# while depth <= maxRec:
			# 	for x in xrange(0,levels_per_colour):
			# 		depth+=1 ## starts at 1
			# 		depth2RGB[depth] =  brewerMap[brewerLevel]
			# 	brewerLevel+=1 ## starts at 0

			for depth in xrange(1,maxRec+1):
				if depth<maxCol:
					depth2RGB[depth] =  brewerMap[depth-1]
				else:
					depth2RGB[depth] =  brewerMap[maxCol-1]


		## apply it to the regions_merged
		for taxa in self.regions_merged:
			for idx,x in enumerate(self.regions_merged[taxa]):
				self.regions_merged[taxa][idx]['RGB'] = depth2RGB[ x['depth'] ]


	def _generate_collapsed_regions(self,rOpen,rClose,endpoint,colour):
		""" generator for the merge_regions function.
			yields blocks where depth>1, i.e. wherever there is any recombination block
			colour=int
		"""
		state=0 # NOT in a region to start with
		for pos in xrange(0,endpoint):
			state_on_entry = state
			if pos in rOpen:
				state+=rOpen[pos]
			if pos in rClose:
				state-=rClose[pos]
			if state_on_entry==0 and state!=0:
				region_start=pos
			elif state==0 and state_on_entry!=0:
				yield {'x1':region_start ,'x2':pos , 'colour':colour}

	def _generate_overlap_regions(self,rOpen,rClose,endpoint,colour):
		""" generator for the merge_regions function.
			yields blocks with based on depth, with the colours showing depth==1 or depth>1
			colour=list len 2
		"""
		state=0 # NOT in a region to start with
		baseCol = int(colour[0])
		overlapCol = int(colour[1])

		for pos in xrange(0,endpoint):
			state_on_entry = state
			if pos in rOpen:
				state+=rOpen[pos]
			if pos in rClose:
				state-=rClose[pos]
			if state_on_entry==0 and state!=0:
				region_start=pos
			elif state==0 and state_on_entry==1: # i.e. singleton block now finished
				yield {'x1':region_start ,'x2':pos , 'colour':baseCol}
			elif state<=1 and state_on_entry>1: # i.e. more than one overlap down to one or zero overlaps
				offset=0 if state==0 else 1 ## see _generate_heat_regions
				yield {'x1':region_start ,'x2':pos-offset , 'colour':overlapCol}
				region_start = pos

	def _generate_heat_regions(self,rOpen,rClose,endpoint):
		""" generator for the merge_regions function.
			yields blocks with the same depth (no colour yet, this comes later)
		"""
		depth=0 # NOT in a region to start with
		## the heatmap colours are just the depth, and are subsequently colourised
		for pos in xrange(0,endpoint):
			prevDepth = depth
			if pos in rOpen:
				depth+=rOpen[pos]
				# print "depth {} -> {} at position {}".format(prevDepth,depth,pos)
			if pos in rClose:
				depth-=rClose[pos]
				# print "depth {} -> {} at position {}".format(depth+1,depth,pos)
			if prevDepth != depth: ## state change
				if prevDepth==0: # i.e. entering a block
					region_start = pos
				else: # leaving a block
					# if the depth hasn't dropped to zero then this block finishes a base before this as the next block will take over at this point
					offset=0
					if depth>0:
						offset=1
					yield {'x1':region_start ,'x2':pos-offset , 'depth':prevDepth}
					region_start = pos


	def missing_data_plot(self,outfile,genome_length,num_taxa,inverse):
		## only called by Rec class, not Multiple Rec
		## two hashs to speed up the hotspot plot creation
		HSplus,HSminus = self._create_breakpoint_hash(self.missing_data_regions)
		count=0
		with open(outfile,'w') as fh:
			fh.write("#BASE\tMISSING\n")
			for pos in xrange(1,int(genome_length)):
				if pos in HSplus:
					count += HSplus[pos]
				if pos in HSminus:
					count -= HSminus[pos]
				frac = float(count)/int(num_taxa)
				if inverse:
					frac = 1-frac
				if frac>0:
					fh.write("{}\t{:.2f}\n".format(pos,frac))

# print "iCANDY.py -d line -w -l 2 -o "+options.tabfile+".recombHotspots.pdf "+options.tabfile+".plot "+options.tabfile+" <EMBL FILE>"
# print "iCANDY.py -d line -w -l 2 -o -q taxa -n "+options.tabfile+".pdf "+options.tabfile+" <EMBL FILE>"




	# def _remove_based_on_taxa(self,excludeset):
	# 	""" removes entries from self.rec and then, finally recalculeates self.regions_by_taxa_unmerged=self._associate_regions_with_taxa()
	# 	"""
	# 	for idx,block in enumerate(self.regions):
	# 		newtaxas = []
	# 		for taxa in block['taxa']:
	# 			if taxa not in excludeset:
	# 				newtaxas.append(taxa)
	# 		self.regions[idx]['taxa'] = newtaxas
	# 		self.regions[idx]['numTaxa'] = len(newtaxas)


	# 	self.regions[:] = [region for region in self.regions if region['numTaxa']!=0]
	# 	self.regions_by_taxa_unmerged=self._associate_regions_with_taxa()


	# def report_block_stats(self):

	# 	print "Recombination blocks identified: {}".format(len(self.regions))

	# 	lens = [ abs(x['x2'] - x['x1']) for x in self.regions  ]
	# 	meanLen = int(np.mean(lens))
	# 	ci = scipy.stats.norm.interval(0.95, loc=meanLen, scale=scipy.stats.sem(lens))
	# 	#return (meanLen,int(ci[0]),int(ci[1]))
	# 	print "Mean block size: {}bp".format(meanLen)
	# 	print "95% CI: {}bp - {}bp".format(int(ci[0]),int(ci[1]))

		# try:
		# 	if self.regions_merged_method!='collapsed':
		# 		self.merge_regions('collapsed')
		# except AttributeError:
		# 	self.merge_regions('collapsed')

		# print "Recombination blocks identified: {}".format(len(self.regions_merged))
		# pdb.set_trace()
		# bases_covered = {}
		# for taxa,regions in self.regions_merged.items():
		# 	bases_covered[taxa]= 0
		# 	for region in regions:
		# 		bases_covered[taxa] += int(region['x2']) - int(region['x1'])

		# meanLen = int(np.mean(bases_covered.values()))
		# ci = scipy.stats.norm.interval(0.95, loc=meanLen, scale=scipy.stats.sem(bases_covered.values()))
		# #return (meanLen,int(ci[0]),int(ci[1]))
		# print "Mean bases covered: {}bp".format(meanLen)
		# print "95% CI: {}bp - {}bp".format(int(ci[0]),int(ci[1]))





	def report_r_m_p_theta(self):
		pass








class MultipleRec(RecBase):
	""" this extends RecBase (just like Rec)
		IN: multiple Rec objects
	"""

	def __init__(self,*arg):
		super(MultipleRec, self).__init__() ## even though this does nothing
		# print "I was called with", len(arg), "arguments:", arg
		self.recs=arg
		### we have to ensure that the recombinations are collapsed (merged)
		for rec in self.recs:
			try:
				if rec.regions_merged_method!='collapsed':
					rec.merge_regions('collapsed')
			except AttributeError:
				rec.merge_regions('collapsed')

		## some basic sanity checking -- debugging purposes mainly
		if len(arg)==2:
			taxa_uniq_to_0, taxa_uniq_to_1, taxa_common = set(), set(), set()
			for taxa in self.recs[0].taxa:
				if taxa in self.recs[1].taxa:
					taxa_common.add(taxa)
				else:
					taxa_uniq_to_0.add(taxa)
			for taxa in self.recs[0].taxa:
				if taxa not in taxa_common and taxa not in taxa_uniq_to_0:
					taxa_uniq_to_1.add(taxa)
			print "{} taxa in common".format(len(taxa_common))
			print "{} only in {}: {}".format(len(taxa_uniq_to_0), self.recs[0].tabfile, taxa_uniq_to_0)
			print "{} only in {}: {}".format(len(taxa_uniq_to_1), self.recs[1].tabfile, taxa_uniq_to_1)
		else:
			print "comparing >2 files is difficult. you have been warned. "




	def merge_tabfiles(self,colours=False):
		"""
		IN: the Rec objects parsed at __init__
		remember gurat?
		produces regions (self.regions_merged of *this* object) with colours indicating which tabfile they came from // which combination of tabfiles contained them
		At maximum, 3 tabfiles can be parsed (7 colours). 4 needs 15 which is too hard.
		"""
		def sanity_check(): ## no self here!
			try:
				assert(len(self.recs)<=3)
			except AssertionError:
				print "Cannot merge more than 3 tabfiles."
				return
			# if self.colours:
			# 	try:
			# 		if len(self.recs)==2:
			# 			assert len(colours)>=3
			# 		elif len(self.recs)==3:
			# 			assert len(colours)>=7
			# 	except AssertionError:
			# 		print "Not enough colours given..."
			# 		return
			### check the taxa are the same in them all...
			# taxa = set(self.recs[0].regions_by_taxa_unmerged.keys())
			# try:
			# 	for obj in self.recs[1:]:
			# 		print "comparing {} to {}".format(taxa,set(obj.regions_by_taxa_unmerged.keys()) )
			# 		pdb.set_trace()
			# 		assert(taxa == set(obj.regions_by_taxa_unmerged.keys()))
			# except AssertionError:
			# 	print "Taxa in the tabfiles are different."
			# 	return

		def brew_colours(numTabs,tabs):
			# The colour format is:
			# tabfiles t1,t2:
			# 	[ t1 , t2 , t1 & t2]
			# 	  0     1       2
			# tabfiles t1,t2,t3:
			# 	[ t1, t2, t3, t1 & t2, t1 & t3, t2 & t3, t1 & t2 & t3 ]
			# 	   0   1   2     3       4         5          6
			if numTabs==2:
				s = "Colours are as follows:\n\t{} -> lime\n\t{} -> green\n\tboth -> blue"
				print s.format(tabs[0],tabs[1])
				return brewer2mpl.get_map("YlGnBu","Sequential",3,reverse=False).colors
			elif numTabs==3:
				s = "Colours are as follows:\n\t{0} -> brown\n\t{1} -> yellow\n\t{2} -> orange\n\t{0} & {1} -> purple\n\t{0} & {2} -> green\n\t{1} & {2} -> blue\n\tall 3 -> red"
				print s.format(tabs[0],tabs[1],tabs[2])
				# return brewer2mpl.get_map("Spectral","Diverging",7,reverse=True).colors
				return brewer2mpl.get_map("Set1","Qualitative",7,reverse=True).colors

		def state_to_colour(state,colours):
			"""state (array of bools) to colour, based on the table in brew_colours"""
			if state  ==[True,False,False]:
				return colours[0]
			elif state==[False,True,False]:
				return colours[1]
			elif state==[False,False,True]:
				return colours[2]
			elif state==[True,True,False]:
				if numTab==2: ## HACK. for the two tabfile situation TTF is the third colour ([2])
					return colours[2]
				return colours[3]
			elif state==[True,False,True]:
				return colours[4]
			elif state==[False,True,True]:
				return colours[5]
			elif state==[True,True,True]:
				return colours[6]

		def custom_overlapper(rHashes,colours,endpoint): ## generator
			state = [False,False,False] ## has what objects which are players in the current position
			## i.e. state = T,F,F means we're in a block from the frist tabfile, not second nor third
			for pos in xrange(0,endpoint):
				state_on_entry = state[:]
				# print "{}. before: {} now: {} sums: {} {}".format(pos,state_on_entry,state,sum(state_on_entry),sum(state))
				### are we entering any blocks at this position? (pair[0][pos] = 1|0 since collapsed)
				for idx,pair in enumerate(rHashes): ## scan the tabfiles (self.rec objects)
					if pos in pair[0]: ## in rOpen
						state[idx]=True
					if pos in pair[1]: ## in rClose
						state[idx]=False
				## if we've gone from [F,F,F] (i.e. nothingness) to something, then we open a block
				if sum(state_on_entry)==0 and sum(state)!=0:
					region_start=pos
				## otherwise a change of state means a change in colour
				elif state!=state_on_entry:
					offset=0 if sum(state)==0 else 1 ## see _generate_hear_regions
					yield {'x1':region_start ,'x2':pos-offset , 'RGB':state_to_colour(state_on_entry,colours)}
					region_start=pos



		print "Merging {} tabfiles...".format(len(self.recs))
		# sanity_check()
		tabfiles = [rec.tabfile for rec in self.recs]
		numTab = len(self.recs)
		colours=brew_colours(numTab,tabfiles) ## users CANNOT specify these at the moment
		self.regions_merged = {taxa:[] for taxa in self.recs[0].regions_by_taxa_unmerged} ## initialise

		for idx, taxa in enumerate(self.regions_merged):
			print "{}/{} {}{}".format(idx,len(self.regions_merged),taxa,"       \r"),
			# print str(taxa)+"               \r",
			sys.stdout.flush()
			## taxa by taxa, create hashes (plus, minus) for EACH rec file
			## i.e. rHashes[0] is (rOpen,rClose) (2 hashes) for the first rec object's regions_merged

			# rHashes = [ self._create_breakpoint_hash(rec.regions_merged[taxa]) for rec in self.recs ]
			# endpoint = max( [ max(pair[1].keys()) for pair in rHashes ] )   +1

			rHashes = [] # block open and close for each tabfile
			endpoint = 0
			for rec in self.recs:
				if taxa in rec.taxa:
					rHashes.append(self._create_breakpoint_hash(rec.regions_merged[taxa]))
					if max(rHashes[-1][1].keys())>endpoint:
						endpoint = max(rHashes[-1][1].keys()) + 1
				else:
					rHashes.append( ({}, {}) ) ## nothing for this taxa

			for region in custom_overlapper(rHashes,colours,endpoint):
				self.regions_merged[taxa].append(region)


		if len(self.recs)==2:
			regions_by_overlap_type = [{}, {}, {}]
			for taxa in self.regions_merged.keys():
				for block in self.regions_merged[taxa]:
					idx = colours.index(block['RGB'])
					try:
						regions_by_overlap_type[idx][taxa].append(block)
					except KeyError:
						regions_by_overlap_type[idx][taxa] = [block]

			for idx in xrange(0,3):
				# pdb.set_trace()
				if idx!=2:
					print "{}".format(tabfiles[idx])
				else:
					print "Both {} and {}".format(tabfiles[0],tabfiles[1])
				tmp = []
				for x in regions_by_overlap_type[idx].values(): tmp.extend(x)
				self.report_block_stats(tmp)


		print ""




class RecAln(object):
	"""
		IN: rec object (subsetted if needed), alignment filename, output filenmae
		FUNCTIONALITY: produce a alignment for beast with many many ?s
	"""

	def __init__(self,rec,alignfile):
		self.rec = rec
		self.alnfile = alignfile
		self.rec.merge_regions('collapsed')
		self.records = [] # this contains the sequences!
		self.taxalist = self.rec.regions_by_taxa_unmerged.keys() ## may have been subsetted
		with open(alignfile,'rU') as fh:
			for record in SeqIO.parse(fh, "fasta"):
				if record.id in self.taxalist:
					self.records.append(record)



	def set_alignment_for_BEAST(self,gapchar='?'):
		print "Removing recombination from the alignment"
		for idx,record in enumerate(self.records):
			self.records[idx].seq = self.records[idx].seq.tomutable()
			taxaname = record.id
			startinglen = len(record.seq)
			blocksprocessed = 0
			try:
				for block in self.rec.regions_merged[record.id]:
					insert = gapchar * (block['x2']+1 - block['x1'])
					self.records[idx].seq[block['x1']:block['x2']+1] = insert
					if len(self.records[idx].seq)!=startinglen:
						print "Len change: {} to {} from block {}..{} ({}bp)".format(startinglen,len(self.records[idx].seq),block['x1'],block['x2'],block['x2']+1 - block['x1'])
						# pdb.set_trace()
						sys.exit()
					blocksprocessed += 1
					# else:
					# 	print "block {}..{} ({}bp) parsed successfuly".format(block['x1'],block['x2'],block['x2']+1 - block['x1'])
			except KeyError:
				pass ## no recombination found for this taxa
			print "\t{} done ({} blocks)".format(taxaname,blocksprocessed)
			self.records[idx].seq = self.records[idx].seq.toseq()

	def write_full_alignment(self,outputalnprefix):
		self.outfile_full = outputalnprefix+".recombremoved.aln"
		print "Saving full alignment with recombinant regions removed to: "+str(self.outfile_full)
		with open(self.outfile_full,'w') as fh:
			SeqIO.write(self.records, fh, "fasta")
		# print "you may want to run \"snp_sites -o {} {}\"".format(outputalnprefix+".snp.aln",outputalnfile)

	def prepare_BEAST_alignment(self,dates,outputalnprefix,treefile):
		from ete2 import Tree
		fainfo = "~jh22/scripts/fainfo"
		fnull = open(os.devnull, "w")
		out_full = self.outfile_full
		out_simon = outputalnprefix+".BEAUTi.aln"
		out_tree = outputalnprefix+".BEAUTi.nex"
		out_taxamap = outputalnprefix+".BEAUTi.name.map"
		out_tmptree = outputalnprefix+".tmp.tre"

		## create the SNP alignment and the patterns block here
		print "preparing BEAST alignment via simons script ~sh16/scripts/BEAST/prepare_BEAST_alignment.py"
		if os.path.exists(out_simon):
			print "\tFile aleady exists. Skipping..."
		else:
			try:
				if dates:
					check_call(" ".join(["~sh16/scripts/BEAST/prepare_BEAST_alignment.py","-x","-n","-d",dates,"-i",out_full,"-o",out_simon]),shell=True)#,stdout=fnull, stderr=fnull)
				else:
					check_call(" ".join(["~sh16/scripts/BEAST/prepare_BEAST_alignment.py","-x","-n","-i",out_full,"-o",out_simon]),shell=True)#,stdout=fnull, stderr=fnull)
			except CalledProcessError:
				print "FAILED"
				return

		## now we (try and) take the tree subset it on the actual names in use here and create a nexus starting tree:
		print "Creating a starting tree from the gubbins tree:"
		try:
			check_call(" ".join([fainfo,"--names","-i",out_simon,"-o",out_taxamap]),shell=True,stdout=fnull)
			newtaxas=[]
			with open(out_taxamap,'rU') as fh:
				for line in fh:
					newtaxas.append(line.strip())
			r = re.compile(r"(.+)_(.+)")
			taxamap = {}
			for newtaxa in newtaxas:
				m = r.match(newtaxa)
				taxamap[m.group(1)] = newtaxa
			# print "taxamap: ",taxamap
			tree = Tree(treefile)
			tree.prune([x for x in tree.get_leaf_names() if x in taxamap.keys()])
			print "tree now has {} taxa".format(len(tree.get_leaf_names()))
			for leaf in tree.iter_leaves():
				leaf.name = "'{}'".format(taxamap[leaf.name])
			tree.resolve_polytomy(recursive=True) ## important for BEAST
			print "Writing tree file to "+out_tree
			tree.write(outfile=out_tmptree)
			with open(out_tmptree,'rU') as fh:
				tree_string = fh.readline().strip()
				# print "tree string: ",tree_string
			with open(out_tree,'w') as fh:
				fh.write('#NEXUS\nbegin taxa;\n\tdimensions ntax={};\n\ttaxlabels\n'.format(len(tree.get_leaf_names())))
				for newtaxa in newtaxas:
					fh.write("\t'{}'\n".format(newtaxa))
				fh.write(';\nend;\n\nbegin trees;\n\ttree gubbins = [&R] ')
				fh.write(tree_string)
				fh.write("\nend;\n")
			check_call(" ".join(["rm",out_taxamap,out_tmptree]),shell=True)
		except:
			print "\tfailed to construct the starting tree. Continuing..."


		print "Compare the differences via\nfainfo -i {}\nfainfo -i {}\nfainfo -i {}".format(self.alnfile,out_full,out_simon)
		print "Now you must create your xml_files (manually) in BEAUti using {}".format(out_simon)
		print "Afterwards, run something like this to introduce the number of constant sites and then run BEAST"
		print 'for x in *xml; do ~sh16/scripts/BEAST/replace_BEAST_blocks.py --no MLE -p '+out_simon+'.patterns -x ${x}; done'
		print 'for i in {1,2,3,4}; do bsub2g -q normal -o bsub.c${i}.txt beast1.8 -seed $RANDOM *_${i}.xml; done'



class RecGenic(object):
	"""
		IN: Rec object, gene struct (dict of dicts. key(outer):geneName, keys(inner):x1,x2,otherStuff)
			list of taxa to restrict output


	"""

	def __init__(self,recObj,emblObj):
		self.rec = recObj
		self.genes = emblObj.genes
		self.associate_recombination_events_with_genes()


	def _overlap(self,a,b):
		if a['x1'] > b['x2'] or a['x2'] < b['x1']:
			return False
		else:
			return True


	def associate_recombination_events_with_genes(self):
		print "Associating recombination events with genes"
		for taxa in self.rec.regions_by_taxa_unmerged:
			## super inefficient method but i don't care
			for region in self.rec.regions_by_taxa_unmerged[taxa]:
				for genename,geneinfo in self.genes.items():
					if self._overlap(geneinfo,region):
						try:
							self.genes[genename]['rcount'] += 1
						except KeyError:
							self.genes[genename]['rcount'] = 1

	def normalise_recombination_counts(self,by=None):
		if by=='numtaxa':
			numtaxa = len( self.rec.regions_by_taxa_unmerged.keys() )
			for genename in self.genes:
				try:
					self.genes[genename]['rcount'] = float(self.genes[genename]['rcount']) / numtaxa
				except KeyError:
					self.genes[genename]['rcount'] = 0

	def write_recombination_counts(self,filename,colheader="rcount"):
		def writeout(a,b):
			fh.write("{}\t{}\n".format(a,b))
		with open(filename,'w') as fh:
			writeout("gene",colheader)
			for genename,geneinfo in self.genes.items():
				try:
					rcount = geneinfo['rcount']
				except KeyError:
					rcount = 0
				writeout(genename,rcount)






class EMBLparser(object):
	""" generic functions cos everything else SUCKS """

	def __init__(self,filepath):
		print "Parsing {} as an EMBL file".format(filepath)
		# self.record = SeqIO.read(filepath,"embl")
		self.db = []
		lines = open(filepath).read().splitlines()
		for block in self._generate_blocks_from_lines(lines):
			# print "\n\nBLOCK RETURNED...\n",block,"\n"
			self.db.append(self._parse_block(block))
		self.genes = {}
		for h in self.db:
			if h['type'] == 'gene':
				if h['gene']:
					self.genes[h['gene']] = h
				elif h['locus_tag']:
					self.genes[h['locus_tag']] = h
				else:
					print "unparsable block!",h


	def _generate_blocks_from_lines(self,text): ## generator
		"""the string "FT   \w" seperates the block. returns (yields) list of lines."""
		buff = []
		for line in text:
			# print "line: {}".format(line)
			if not line.startswith('FT'):
				continue
			if not line.startswith('FT    '):
				# print "...NEW BLOCK DETECTED..."
				if buff:
					yield buff
				buff = [line]
			else:
				# print "...appending to stack..."
				buff.append(line)
		if buff:
			yield buff


	def _parse_block(self,block):
		"""turns list of lines (i.e. block) into hash of salient features"""
		def _x1x2parser(line):
			sec = line.split()[2].strip()
			r = re.compile(r"(\d+)\.\.(\d+)")
			m = r.findall(line)
			if len(m) == 1:
				return ( int(m[0][0]) , int(m[0][1]) )
			else:
				return ( int(m[0][0]) , int(m[-1][1]) )
		h = {'type':False, 'x1':False, 'x2':False, 'strand':False, 'locus_tag':False, 'gene':False}
		h['type'] = block[0].strip().split()[1]
		if 'complement' in block[0]:
			h['strand'] = '-'
		else:
			h['strand'] = '+'
		h['x1'],h['x2'] = _x1x2parser(block[0])
		for line in block:
			if '/locus_tag=' in line:
				h['locus_tag'] = line.split('"')[1]
			elif '/gene=' in line:
				h['gene'] = line.split('"')[1]
		# print "BLOCK PARSED AS:",h
		return h



	def write_gene_positions(self,outfile,selectedtype='gene'):
		genes = {}
		for h in self.db:
			if h['type'] == selectedtype:
				if h['gene']:
					genes[h['gene']] = h
				elif h['locus_tag']:
					genes[h['locus_tag']] = h
				else:
					print "unparsable block!",h
		with open(outfile,'w') as fh:
			fh.write("gene,x1,x2,strand\n")
			for key,value in genes.items():
				fh.write("{},{},{},{}\n".format(key,value['x1'],value['x2'],value['strand']))




