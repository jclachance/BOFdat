def metxref2dict(xref_file="chem_xref.tsv"):
	""" simple function to parse the metanetx file `chem_xref.tsv` and return 2 dictionaries, `met2mnx` and `mnx2met`
	
		met2mnx: {db_keyword -> {db_met_id -> mnx_id}}
		mnx2met: {mnx_id -> {db_keyword -> {[met_ids]}}
	
	"""
	met2mnx = {}
	mnx2met = {}
	with open(xref_file) as fh:
		for line in fh:
			if line.startswith("#"): continue
			met_id,mnx_id = line.split()[:2]
			try: db_id,met_id = met_id.split(":",1)
			except ValueError: continue #print "!!!",line; return met_id,mnx_id 
			if db_id not in met2mnx.keys(): met2mnx[db_id] = {}
			try: met2mnx[db_id][met_id]
			except: met2mnx[db_id][met_id] = mnx_id
			else: print "met2mnx: duplicate found in db %s for %s! %s and %s" %(db_id,met_id,met2mnx[db_id][met_id],mnx_id)
			mnx2met[mnx_id] = mnx2met.get(mnx_id,{})
			mnx2met[mnx_id][db_id] = mnx2met[mnx_id].get(db_id,[]) + [met_id]
			#except: mnx2met[mnx_id][db_id] = met_id
			#else: print "mnx2met: duplicate found in db %s for %s! %s and %s" %(db_id,mnx_id,mnx2met[mnx_id][db_id],met_id)
	return met2mnx,mnx2met
			
def locxref2dict(xref_file="comp_xref.tsv"):
	""" simple function to parse the metanetx file `comp_xref.tsv` and return 2 dictionaries, `loc2mnx` and `mnx2loc`
	
		loc2mnx: {db_keyword -> {db_loc_id -> mnx_id}}
		mnx2loc: {mnx_id -> {db_keyword -> {[loc_ids]}}
	
	"""
	loc2mnx = {}
	mnx2loc = {}
	with open(xref_file) as fh:
		for line in fh:
			if line.startswith("#"): continue
			loc_id,mnx_id = line.split()[:2]
			try: db_id,loc_id = loc_id.split(":",1)
			except ValueError: continue # print "!!!",line; return loc_id,mnx_id 
			if db_id not in loc2mnx.keys(): loc2mnx[db_id] = {}
			try: loc2mnx[db_id][loc_id]
			except: loc2mnx[db_id][loc_id] = mnx_id
			else: print "loc2mnx: duplicate found in db %s for %s! %s and %s" %(db_id,loc_id,loc2mnx[db_id][loc_id],mnx_id)
			mnx2loc[mnx_id] = mnx2loc.get(mnx_id,{})
			mnx2loc[mnx_id][db_id] = mnx2loc[mnx_id].get(db_id,[]) + [loc_id]
	return loc2mnx,mnx2loc
	
class Xref_handler(object):
	def __init__(self,chem_xref_file="chem_xref.tsv",comp_xref_file="comp_xref.tsv"):
		self.chem_xref_file = chem_xref_file
		self.comp_xref_file = comp_xref_file
		self.chem2mnx,self.mnx2chem = metxref2dict(self.chem_xref_file)
		self.comp2mnx,self.mnx2comp = locxref2dict(self.comp_xref_file)

	def crossDbQuery(self,db1,db2,my_id,type_="chem",return_list=False):
		""" given the dictionaries build using functions such as `metxref2dict`, the metabolite `met_id` is searched from db1 to db2.
			
			db1: database from which met_id should be taken (mnx, kegg, metacyc, lipidmaps, upa, seed, hmdb, biopath, reactome, umbbd, chebi, bigg)
			db2: database in which the equivalent of met_id should be searched
			met_id: the id of the metabolite from db1
			type_: the entity class which will be searched (i.e. chem or comp)
			return_list: if multiple ids are found in db2, only the first is returned. Set this option as true to return the whole list.
			
			Examples
			-------
			crossDbQuery("bigg","seed","atp") -> cpd00002
			crossDbQuery("bigg","seed","c",type_="comp") -> c0
		
			"""
		d1,d2 = map(lambda x: getattr(self,x), ["%s2mnx" %type_,"mnx2%s" %type_])
		
		try: mnx_id = d1[db1][my_id]
		except KeyError: raise KeyError("%s not found in %s! Check again and retry..." %(my_id,db1))
		try: new_id = d2[mnx_id][db2]
		except KeyError: raise KeyError("no equivalents of %s (mnx id: %s) found in %s" %(my_id,mnx_id,db2))
		if len(new_id) == 1:
			if not return_list: return new_id[0]
			else: return new_id
		else:
			print "WARNING: multiple elements found matching with %s in db %s!" %(new_id,db2)
			if return_list:
				print "All the elements will be returned"
				return new_id
			print "Only the last element will be returned!"
			return new_id[-1]		

	def translateBiggMetId(self,the_id,db,alternative_names=False,check_model=None):
		""" higher-level fxn to translate a full bigg id (eg datp_c) with the nomenclature of another db
		
			the_id: id of metabolite and compartment (eg datp_c)
			db: name of the db in which to find the equivalent of `the_id` (eg seed)
			alternative_names: by default, only on element is returned. if set true, a list of all the multiple names (either for metabolite or location,combined)
		 """
		import re
		met,loc = re.search("(.*)_(i?[a-z])$",the_id).groups()
		if alternative_names:
			from itertools import product
			new_mets = self.crossDbQuery("bigg",db,met,return_list=True)
			new_loc = self.crossDbQuery("bigg",db,loc,type_="comp",return_list=True)
			all_mets = ["%s_%s" %p for p in product(new_mets,new_loc)]
			if check_model != None: all_mets = filter(lambda x: x in check_model.metabolites,all_mets)
			return all_mets
		else:
			new_met = self.crossDbQuery("bigg",db,met)
			new_loc = self.crossDbQuery("bigg",db,loc,type_="comp")
			return ["%s_%s" %(new_met,new_loc)]
		
	def biggBiomass2mnxs(self,db,out_prefix=".",check_model=None,brute_force=False):
		""" convenient wrapper to translate all the bigg ids of biomass constituents used in BOFdat 
		
			db: name of db in which to find the equivalent of bigg metabolites
			out_prefix: directory in which to write the conversion tables
			brute_force: by default, when aliases are found in `db` for given metabolites/locations, only the first is returned.
						With `brute_force=True`, all aliases will be searched in a `check_model` and will be written only if present. Require a valid `check_model`  
			check_model: a model in which to search the translated compounds. Required if `brute_force=True`
		
		TODOs
		-------
		improve the model import (only read_sbml_model from cobra.io)
		
		"""
		import os
		from cobra.io import read_sbml_model
		# bigg names fetched directly from BOFdat/core/ files (eg grep -Po '(?<=metabolites\.)([a-zL_]*)(?=_c)' projects/BOFdat/BOFdat/core/protein.py)
		dntps = ["datp_c","dctp_c","dgtp_c","dttp_c","ppi_c"]
		ntps = ["atp_c","ctp_c","gtp_c","ttp_c","ppi_c"]
		aas = ["ala__L_c","cys__L_c","asp__L_c","glu__L_c","phe__L_c",
				"gly_c","his__L_c","ile__L_c","lys__L_c","leu__L_c",
				"met__L_c","asn__L_c","pro__L_c","gln__L_c","arg__L_c",
				"ser__L_c","thr__L_c","val__L_c","trp__L_c","tyr__L_c"]
				
		to_check = brute_force
		if to_check:
			if check_model == None: raise Exception("with brute_force=True a check model should be supplied")
			else:
				if type(check_model) == str: check_model = read_sbml_model(check_model)
				elif type(check_model) == cobra.core.model.Model: pass
				else: raise Exception("`check_model` should be either a valid .sbml file or a Model object")
				
		for l in zip([dntps,ntps,aas],["dntps","ntps","aas"]):
			met_list,suffix = l
			out_name = os.path.join(out_prefix,"%s.conversion_table.txt" %suffix)
			with open(out_name,'w') as wh:
				wh.write("bigg_id\t%s_id\n" %db)
				for met in met_list:
					new_mets = self.translateBiggMetId(met,db,alternative_names=to_check,check_model=check_model)
					for m in new_mets: wh.write("%s\t%s\n" %(met,m))

if __name__ == "__main__":

	usage = """
- USAGE
	python metaNetX_xref_lib.py db out_prefix [check_model] [xref_prefix]

- DESCRIPTION
	this script is a collection of functions to translate ids of compounds through the MetaNetX db.
	In the present form, it produces a set of *conversion_table.txt files to be used with BOFdat step1 fxns that use hardcoded bigg ids. 
	For this reason the id translation is from bigg to a gived db of interest. 
	Since multiple synonymous ids can be present, the user can provide a model to keep the synonyms present in the model.
	The xref files from MetaNetX should be available (can be downloaded from https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv and https://www.metanetx.org/cgi-bin/mnxget/mnxref/comp_xref.tsv)

	db: database used for the translation
	out_prefix: path in which output files (*conversion_table.txt) are written
	check_model: (OPTIONAL, default None) a .sbml file. If given, all possible synonyms of bigg metabolites in `db`, which are also present in the model, are reported. 
	xref_prefix: (OPTIONAL, default ./) the path in which the MetaNetX files are located 

- EXAMPLE INPUT:
	python metaNetX_xref_lib.py seed ./ ./toy_model.sbml 

- EXAMPLE OUTPUT (head *conversion_table.txt):

	==> aas.conversion_table.txt <==
	bigg_id	seed_id
	ala__L_c	cpd00035_c0
	cys__L_c	cpd00084_c0
	asp__L_c	cpd00041_c0
	glu__L_c	cpd00023_c0
	phe__L_c	cpd00066_c0
	gly_c	cpd00033_c0
	his__L_c	cpd00119_c0
	ile__L_c	cpd00322_c0
	lys__L_c	cpd00039_c0

	==> dntps.conversion_table.txt <==
	bigg_id	seed_id
	datp_c	cpd00115_c0
	dctp_c	cpd00356_c0
	dgtp_c	cpd00241_c0
	ppi_c	cpd00012_c0

	==> ntps.conversion_table.txt <==
	bigg_id	seed_id
	atp_c	cpd00002_c0
	ctp_c	cpd00052_c0
	gtp_c	cpd00038_c0
	ttp_c	cpd00357_c0
	ppi_c	cpd00012_c0

"""

	import sys,os
	args = sys.argv[1:]
	check_model,xref_prefix = None,'.'
	db,out_prefix = args[:2]
	if len(args) == 3:
		if os.path.isdir(args[2]): xref_prefix = args[2]
		else: check_model = args[2]
	elif len(args) == 4: check_model,xref_prefix = args[-2:]
	else: print usage; print "wrong usage! RETRY!!"; exit()
	metanetx_xrefs = map(lambda x: os.path.join(xref_prefix,x),["chem_xref.tsv","comp_xref.tsv"])
	xref = Xref_handler(*metanetx_xrefs)
	is_brute_force = check_model != None
	xref_handler.biggBiomass2mnxs(db,check_model=check_model,brute_force=is_brute_force)




