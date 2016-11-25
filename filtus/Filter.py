import FiltusUtils
import FiltusWidgets
import VariantFileReader

import itertools
import collections
import time
from operator import itemgetter

class Filter(object):
    def __init__(self,
                model = None,
                exclude_var_txt = None,
                restrict_var_txt = None,
                exclude_genes_txt = None,
                restrict_genes_txt = None,
                regions_txt = None,
                exclude_variants = None,
                restrict_to_variants = None,
                exclude_genes = None,
                restrict_to_genes = None,
                columnfilters = None,
                benignPairs = None,
                controls = None,
                regions = None,
                closePairLimit = 0,
                filterFile = None,
                filtus=None):
        
        if filterFile: # raises error if something goes wrong
            restrict_genes_file, exclude_genes_file, exclude_var_file, regions_file, columnfilters_file = self.readFromFile(filterFile)
            if restrict_genes_file: restrict_genes_txt = '; '.join(x for x in [restrict_genes_txt, restrict_genes_file] if x)
            if exclude_genes_file: exclude_genes_txt = '; '.join(x for x in [exclude_genes_txt, exclude_genes_file] if x)
            if exclude_var_file: exclude_var_txt = '; '.join(x for x in [exclude_var_txt, exclude_var_file] if x)
            if columnfilters_file: columnfilters = columnfilters_file + (columnfilters if columnfilters else [])
            
        errormessage = 'Error in filter initialization.\n\n'
        if exclude_genes_txt is not None:
            if exclude_genes is not None:
                raise ValueError(errormessage + "At least one of 'exclude_genes_txt' and 'exclude_genes' must be None")
            else:
                exclude_genes = FiltusWidgets.GeneFile.read(exclude_genes_txt)

        if restrict_genes_txt is not None:
            if restrict_to_genes is not None:
                raise ValueError(errormessage + "At least one of 'restrict_genes_txt' and 'restrict_to_genes' must be None")
            else:
                restrict_to_genes = FiltusWidgets.GeneFile.read(restrict_genes_txt)

        if exclude_var_txt is not None:
            if exclude_variants is not None:
                raise ValueError(errormessage + "At least one of 'exclude_var_txt' and 'exclude_variants' must be None")
            elif filtus is not None:
                if exclude_var_txt in filtus.storage:
                    exclude_variants = filtus.storage[exclude_var_txt]
                else:
                    exclude_variants = FiltusWidgets.VariantFile.read2(exclude_var_txt, filtus)
            else:
                dbReader = VariantFileReader.VariantFileReader()
                exclude_variants = dbReader.readNonVCF(exclude_var_txt, skiplines=4, sep="\t", chromCol="CHROM", posCol="POS", geneCol=None, gtCol=None).getUniqueVariants()
                
        if restrict_var_txt is not None:
            if restrict_to_variants is not None:
                raise ValueError(errormessage + "At least one of 'restrict_var_txt' and 'restrict_to_variants' must be None")
            elif filtus is not None:
                restrict_to_variants = FiltusWidgets.VariantFile.read2(restrict_var_txt, filtus)
            else:
                dbReader = VariantFileReader.VariantFileReader()
                exclude_variants = dbReader.readNonVCF(restrict_var_txt, skiplines=4, sep="\t", chromCol="CHROM", posCol="POS", geneCol=None, gtCol=None).getUniqueVariants()
                
        if regions_txt is not None:
            if regions is not None:
                raise ValueError(errormessage + "At least one of 'regions_txt' and 'regions' must be None")
            else:
                regions = FiltusWidgets.RegionFile.read(regions_txt)
        
        self.restrict_genes_txt, self.exclude_genes_txt, self.restrict_var_txt, self.exclude_var_txt, self.regions_txt = \
                restrict_genes_txt, exclude_genes_txt, restrict_var_txt, exclude_var_txt, regions_txt
        self.filtus = filtus
        self.model = model if model is not None else 'Dominant'
        self.columnfilters_original = columnfilters
        self.benignPairs = None if benignPairs is True and not controls else benignPairs
        self.regions = regions
        self.closePairLimit = closePairLimit
        self.controls = controls
        self.exclude_variants,     self.exclude_variants_original     = exclude_variants,     exclude_variants
        self.restrict_to_variants, self.restrict_to_variants_original = restrict_to_variants, restrict_to_variants
        self.exclude_genes,        self.exclude_genes_original        = exclude_genes,        exclude_genes
        self.restrict_to_genes,    self.restrict_to_genes_original    = restrict_to_genes,    restrict_to_genes
        
        self.cleanup()
        
        
    def cleanup(self):
        self.columnfilters = self.prepareColumnfilters(self.columnfilters_original)
        self.regionsChromdic = self.prepareRegions(self.regions)
        
        exclvar, inclvar = self.exclude_variants_original, self.restrict_to_variants_original
        exclgen, inclgen = self.exclude_genes_original, self.restrict_to_genes_original
        controls, model = self.controls, self.model
        
        if exclvar and inclvar:
            self.restrict_to_variants = inclvar.difference(exclvar)
            self.exclude_variants = None
        
        if exclgen and inclgen:
            self.restrict_to_genes = inclgen.difference(exclgen)
            self.exclude_genes = None
        
        if controls:
            #Dominant: Remove all variants in controls. Both recessive: Remove all homozygous variants. Compound heterozygous: Remove benign pairs.
            if model == 'Dominant':
                self.exclude_variants = set.union(*[VF.getUniqueVariants() for VF in controls])
            else:
                self.exclude_variants = set.union(*[VF.getUniqueVariants(alleles=2) for VF in controls])
            if model == 'Recessive' and self.benignPairs is True: # i.e. this is only included if 'benignPairs = True'
                self.benignPairs = self.extractBenignPairs(controls)

    def prepareRegions(self, regions):
        if not regions: return None
        regions.sort()
        chromdic = collections.defaultdict(list)
        for reg in regions:
            chromdic[reg[0]].extend([float(reg[1]), float(reg[2])+1]) # +1 ensures endpoint in included in region, see filtusUtils.firstGreater().
        return chromdic
     
     
    def __str__(self):
        txtlist = self.details()
        return '## ' + '\n## '.join(txtlist) if txtlist else "## No filters"

    def details(self):
        txtlist = []
        if self.closePairLimit: txtlist.append('Close variants removal:' + str(self.closePairLimit))
        if self.restrict_genes_txt: 
            txtlist.append('Gene restriction: %s (%d genes)' %(self.restrict_genes_txt, len(self.restrict_to_genes_original)))
        elif self.restrict_to_genes:
            n = len(self.restrict_to_genes)
            txtlist.append('Gene restriction: ' + (', '.join(self.restrict_to_genes) if n < 7 else str(n) + ' genes'))
        
        if self.exclude_genes_txt: 
            txtlist.append('Gene exclusion: %s (%d genes)' %(self.exclude_genes_txt, len(self.exclude_genes_original)))
        elif self.exclude_genes:
            n = len(self.exclude_genes)
            txtlist.append('Gene exclusion: ' + (', '.join(self.exclude_genes) if n < 7 else str(n) + ' genes'))
        
        if self.exclude_var_txt: txtlist.append('Variant exclusion: %s (%d variants)' %(self.exclude_var_txt, len(self.exclude_variants_original)))
        if self.restrict_var_txt: txtlist.append('Restrict to variants: %s (%d variants)' %(self.restrict_var_txt, len(self.restrict_to_variants_original)))
        if self.regions_txt: txtlist.append('Restrict to regions: %s (%d regions)' %(self.regions_txt, len(self.regions)))
        
        if self.columnfilters_original:
            txtlist.append("Column filters:")
            txtlist.extend(['  ' + ' ::: '.join(map(str, cf)) for cf in self.columnfilters_original])
        if self.model != "Dominant":
            txtlist.append("Genetic model: " + self.model)
        return txtlist
        
    @classmethod # making it available without a Filter instance: Used in geneSharingFamily.
    def extractBenignPairs(self, VFcontrols):
        comb = itertools.combinations
        benignPairs = collections.defaultdict(set)
        for i, VF in enumerate(VFcontrols):
            for gene, vars in VF.geneDict(addIndex=i).iteritems():
                pairs = comb(vars.getUniqueVariants(alleles=1), 2) #tuples
                benignPairs[gene].update(map(frozenset, pairs)) #need frozenset to make set of sets.
        return benignPairs

    def readFromFile(self, filename):
        with open(filename, 'rU') as f:
            lines = [line.strip() for line in f if line.strip()]
        L = len(lines)
        if L < 6:
            raise IOError("Format error (too few lines).")
        for s in lines[:4]:
            if not (s.startswith('0') or s.startswith('1')):
                raise ValueError("Expected string with initial character 0 or 1, but recieved\n'%s'" % s)
        Ncolfilt = int(lines[4])
        if 5 + Ncolfilt != L:
            raise ValueError("Wrong number of column filters:\n" \
            "Number indicated is %d, but the file contains %d" % (Ncolfilt, L - 5))

        restrict_genes_txt = lines[0][2:] if lines[0].startswith('1') and lines[0][2:].strip() else None
        exclude_genes_txt  = lines[1][2:] if lines[1].startswith('1') and lines[1][2:].strip() else None
        exclude_var_txt    = lines[2][2:] if lines[2].startswith('1') and lines[2][2:].strip() else None
        regions_txt        = lines[3][2:] if lines[3].startswith('1') and lines[3][2:].strip() else None
        columnfilters = []
        for k in range(5, 5 + Ncolfilt):
            col, rel, val, keep = lines[k].split(' ::: ')
            if ('less' in rel or 'greater' in rel):
                try:
                    val = float(val)
                except ValueError:
                    FiltusUtils.warningMessage("Column filter ignored:\n\n'%s  %s  %s'\n\nNumerical value needed."%(col, rel, val))
                    continue
            columnfilters.append((col, rel, val, int(keep)))
        return restrict_genes_txt, exclude_genes_txt, exclude_var_txt, regions_txt, columnfilters

    def prepareColumnfilters(self, cfs):
        if not cfs: return None
        # litt tungvint aa gjenta operatornavn her. Grunner: reltext har prioritert rekkefoelge, samt at Filter-klassen ikke har referanse til filtus.FM.
        reltexts = ['equal to', 'missor_equal to', 'not equal to', 'missor_not equal to', 'greater than', 'missor_greater than',
                        'less than', 'missor_less than', 'starts with', 'missor_starts with', 'does not start with', 'missor_does not start with',
                        'contains', 'missor_contains', 'does not contain', 'missor_does not contain']
        opnames = [getattr(FiltusUtils, op) for op in ['myeq', 'missor_eq', 'myne', 'missor_ne', 'floatgt', 'missor_floatgt', 'floatlt',
                        'missor_floatlt', 'mystartswith', 'missor_mystartswith', 'not_mystartswith', 'missor_not_mystartswith', 'contains',
                        'missor_contains', 'not_contains', 'missor_not_contains']]
        opDic = dict(zip(reltexts, opnames))

        def rel2op(rel, keepMissing):
            if keepMissing: rel = 'missor_' + rel
            return opDic[rel]

        #Sort: Always do regex last; AND/OR filters before those, and all others first, sorted according to order in reltexts.
        cfs_sorted = sorted(cfs, key = lambda x: (str(x[2]).startswith('REGEX'), ' AND ' in str(x[2]) or ' OR ' in str(x[2]), reltexts.index(x[1]))) # val may have been converted to float

        # translate from relation text to operator
        return [(col, rel2op(rel, keep), val, keep) for (col, rel, val, keep) in cfs_sorted]


    def checks(self, VF):
        errors = []
        if VF.geneGetter is None and (self.exclude_genes or self.restrict_to_genes or self.model == "Recessive"):
            errors.append("I don't know which column contains gene names.")
        if VF.genotypeGetter is None and self.model != 'Dominant':
            errors.append("I don't know which column contains genotype information.")
        if errors:
            raise ValueError('\n\n'.join(["Filtering error in the following sample:\n%s"% VF.longName] + errors))

        if self.columnfilters:
            col = [cf[0] for cf in self.columnfilters_original if not cf[3] and cf[0] not in VF.columnNames]
            if col:
                plur = len(col)>1
                cols = '\n'.join(col)
                xx = 'heading %s does'%cols if len(col) == 1 else 'headings %s do'%cols
                message = "'%s' does not have column%s '%s'. No variants will survive filtering." %(VF.shortName, 's' if plur else '', "', '".join(col)) \
                             +"\n\nTip: If you check the 'keep if missing' box next to a column filter, it will be ignored for all files lacking the corresponding column."
                raise ValueError(message)

    def apply(self, VF, checks = True, inplace=False):
        if checks:
            try:
                self.checks(VF)
            except ValueError as message:
                FiltusUtils.warningMessage(message)
                if inplace: 
                    VF.setVariants([])
                    VF.appliedFilters = self
                    return
                else:
                    return VF.copyAttributes(variants=[], appliedFilters=self)

        headers = VF.columnNames
        columnfilters = self.columnfilters
        exclude_variants = self.exclude_variants
        restrict_to_variants = self.restrict_to_variants
        exclude_genes = self.exclude_genes
        restrict_to_genes = self.restrict_to_genes
        regionsChromdic = self.regionsChromdic
        res = VF.variants[:]

        ### 0. if removeClosePairs - do this first
        if self.closePairLimit > 0:
            res[:] = removeClosePairs(VF, minDist = self.closePairLimit, variants_only=inplace) #TODO

        ### 1. restriction filters. Usually not used ###
        if restrict_to_genes: #Doing this first because: Usually a very small set, thus reducing the variant set substantially
            annGenes = VF.annotatedGenes
            res[:] = [v for v in res if any(g in restrict_to_genes for g in annGenes(v))]
        if restrict_to_variants: #Not yet implemented in the GUI, but used in Family Gene Sharing ((stemmer det??)
            varDef = VF.varDefGetter
            res[:] = [v for v in res if varDef(v) in restrict_to_variants]

        ### 2. exclude variants. Usually present, and usually gives large reduction.
        if exclude_variants:
            varDef = VF.varDefGetter
            res[:] = [v for v in res if varDef(v) not in exclude_variants]

        ### 3. column filters: These are already sorted w.r.t. speed ###
        if columnfilters:
            for col, op, entry, keep in columnfilters:
                if col in headers:
                    getcol = itemgetter(headers.index(col))
                    res[:] = [v for v in res if op(getcol(v), entry)]
                elif keep: continue
                # else: This is dealt with in checks()

        ### 4. regions: ###
        if regionsChromdic is not None:
            firstGreater = FiltusUtils.firstGreater
            chrom, pos = VF.chromGetter, VF.posGetter
            chrom_vars = collections.defaultdict(list)
            for v in res:
                chrom_vars[chrom(v)].append(v)
            res = []
            for chr in sorted(chrom_vars.keys(), key=FiltusUtils.chromInt):
                reg_startstops = regionsChromdic[chr]
                vars = chrom_vars[chr]
                vars.sort(key=lambda v: float(pos(v))) # possible to avoid doing float(pos()) twice??
                interv = firstGreater((float(pos(v)) for v in vars), reg_startstops) # --> odd/even if inside/outside some region!
                res.extend(v for v, intv in zip(vars, interv) if intv % 2 == 1)
     
     
        ### 5. exclude genes. Might be relatively slow because of annGenes. Haven't checked this though. ###
        if exclude_genes:
            annGenes = VF.annotatedGenes
            res[:] = [v for v in res if not any(g in exclude_genes for g in annGenes(v))]

        ### 6. model filter and controls (benign pairs in compound rec models) ###
        comb = itertools.combinations
        model = self.model
        if model == 'Recessive homozygous':
            isHom = VF.isHomALT()
            res[:] = [v for v in res if isHom(v)]
        elif model == 'Recessive':
            varDef = VF.varDefGetter
            annGenes = VF.annotatedGenes
            isHom = VF.isHomALT()
            benignPairs = self.benignPairs
            if benignPairs:
                remov = set()
                heteroDict = collections.defaultdict(set) #dict of heterozygous variants
                for v in res:
                    if not isHom(v):
                        vdef = varDef(v)
                        for g in annGenes(v):
                            heteroDict[g].add(vdef)
                # the following ignores that variants can be annotated with multiple gene. TODO
                for g, vars in heteroDict.iteritems():
                    if g in benignPairs:
                        nonBenign = set(frozenset(pair) for pair in comb(vars, 2)) - benignPairs[g]
                        remov.update(vars - set(vdef for pair in nonBenign for vdef in pair))
                res[:] = [v for v in res if varDef(v) not in remov]
            geneCount = collections.Counter(g for w in res for g in annGenes(w))
            res[:] = [v for v in res if isHom(v) or any(geneCount[g] > 1 for g in annGenes(v))]
        if inplace: 
            VF.setVariants(res)
            VF.appliedFilters = self
            return
        else:
            return VF.copyAttributes(variants=res, appliedFilters=self)


    
    def prepareTablePrompt(self):
        names = ['Candidate genes', 'Exclude genes', 'Exclude variants', 'Restrict to regions']
        argnames = ['restrict_to_genes', 'exclude_genes', 'exclude_variants', 'regions']
        variables = [self.restrict_to_genes_original, self.exclude_genes_original, self.exclude_variants_original, self.regions]
        filterlist = [(name, filter) for name, filter in zip(argnames, variables) if filter] + [('columnfilters', cf) for cf in self.columnfilters_original]

        def filterString(argname, filter):
            if argname == 'columnfilters':
                col, rel, val, keep = filter
                if keep: rel = 'missing or ' + rel
                return '%s %s "%s"' %(col, rel, val)
            else:
                return names[argnames.index(argname)]

        filtertext = [filterString(argname, filter) for argname, filter in filterlist]
        prompt_text = 'Available filters:\n\n' + '\n'.join('     %d: ' %(i + 1,) + string for i, string in enumerate(filtertext))
        prompt_text += '\n\nUse the numbers above to indicate the order of filtering steps, separated\nby comma. Use "&" to combine several filters in a single step. If you start\nwith a comma, the first step is unfiltered counts.\n\nExample: If you enter "3,  1 & 2,  4" then the first step consists of filter 3, the \nsecond of filters 1 and 2 combined, and in the last step filter 4.\n\n\nFiltering steps:'
        return [filterlist, filtertext, prompt_text]

def trioRecessiveFilter(VFch, VFfa, VFmo, model):
    '''Input: 
    VFch, VFfa, VFmo: VCFtypeData objects which must have been loaded from a single joint calling variant file, with keep00=1.
    model: Either 'Recessive homozygous' or 'Recessive' (which includes both homozygous and compound heterozygous).
    '''
    
    ### Setup
    # Local function referances for speed
    vardef = VFch.chromPosRefAlt
    annGenes = VFch.annotatedGenes
    XminusPAR = FiltusUtils.XminusPAR
    
    # Utility function for extracting the genotype (as a tuple (a0, a1)) of a single variant
    # Numerical codes (as in GT field) are used: Important that all family members are in the same VCF file.
    alleles = VFch.alleles(actualAlleles=False)
    
    # Store observed genotypes of all variants in each parents
    fatherAL = {vardef(v) : alleles(v) for v in VFfa.variants} 
    motherAL = {vardef(v) : alleles(v) for v in VFmo.variants} 
    
    # The set of all variants where at least one parent is homozygous. Entries: ((chrom, pos, ref, alt), A) where A is the observed homozygous allele.
    parHOM = {(vdef,a[0]) for vdef, a in itertools.chain(fatherAL.iteritems(), motherAL.iteritems()) if a[0]==a[1]!='.'}
    
    ### Loop through the variants of the child and check compatibility with AR or XR inheritance
    if model == 'Recessive homozygous':
        # homAllele is a function which takes a single variant as input and returns A if homozygous A/A, or -1 if heterozygous.
        homAllele = VFch.homAllele()
        res = []
        for v in VFch.variants:
            a = alleles(v)
            a0 = a[0]
            
            # if not homozygous or REF/REF or missing, skip to next
            if (not a[1] == a0) or a0 == '0' or a0 == '.':
                continue
                
            vdef = vardef(v)
                
            # Check if homoz in either parent (this step is for speed, could be skipped)
            if (vdef, a0) in parHOM: 
                continue
            
            # Collect the parental alleles. Error imply missing data. If so: dont include.
            try:
                father_alleles = fatherAL[vdef]
                mother_alleles = motherAL[vdef]
            except KeyError:
                continue
            
            # Check first if mother is not heterozygous (this is a requirement for both AR and XR inheritance)
            if not sum(mo == a0 for mo in mother_alleles):
                continue
            
            # Now check the number of copies carried by the father
            father_num = sum(fa == a0 for fa in father_alleles)
            if father_num == 0 and XminusPAR(vdef[:2]):
                # X-linked hemizygosity: Include if child homozyous & mother heterozygous & father is non-carrier
                # NB: Only "correct" XR for male child. If child is female, this situation requires a de novo event.
                res.append(v)
            elif father_num == 1:             
                # Autosomal homozygosity: Include if both parents are heterozygous
                res.append(v)
        
                
    elif model == 'Recessive':
        
        # Storage for variants compatible with simple recessive (homozygous) inheritance
        homoz = []
        
        # Dictionary for gene-wise storing of potential comp. het. variants. (Needed to apply rule 4.)
        chGene = collections.defaultdict(list)
        
        # Dictionary for gene-wise storing of parental origins of the comp.het variants. (Needed to apply rule 5.)
        fromParent = collections.defaultdict(set)
        
        for v in VFch.variants:
            a = alleles(v)
            
            # If REF/REF or missing, skip to next variant.
            if a == ('.','.') or a == ('0','0'):
                continue
            
            vdef = vardef(v)
            
            # Collect the parental alleles. Error imply missing data. 
            try:
                father_alleles = fatherAL[vdef]
                mother_alleles = motherAL[vdef]
            except KeyError:
                continue
                
            ### Part 1: Homozygous ###
            
            # If homozygous in child: store if heterozygous in both parents (AR) or heteroz in mother and missing in father (XR).
            # In any case, continue to next variant.
            if a[0] == a[1]:
                a0 = a[0]
                if (vdef, a0) in parHOM: #check if homoz in either parent (this step is for speed, could be skipped)
                    continue
                
                # If mother is not heterozygous, continue (this is a requirement for both AR and XR inheritance)
                if not sum(mo == a0 for mo in mother_alleles) == 1:
                    continue
            
                # Now check the number of copies carried by the father
                father_num = sum(fa == a0 for fa in father_alleles)
                if XminusPAR(vdef[:2]):
                    # X-linked hemizygosity: Include if mother heterozygous & father is non-carrier
                    # NB: Only "correct" XR for male child. If child is female, this situation requires a de novo event.
                    if father_num == 0:
                        homoz.append(v)
                else:
                    # Autosomal homozygosity: Include if both parents are heterozygous
                    if father_num == 1:             
                        homoz.append(v)
                    
                continue
        
            ### Part 2: Compound heterozygous. Use the 5 rules from Kamphans & al. ###
            
            # Rule 1: the child must be heterozygous. This is automatically fulfilled by now.
            # Rule 2: unaffected indivs must not be homozygous. In the trio case, this is covered by Rule 3 below.
            
            for a1 in a:
                if a1 == '0': 
                    continue
                
                # count allele copies in each parent
                father_num = sum(fa == a1 for fa in father_alleles)
                mother_num = sum(mo == a1 for mo in mother_alleles) 
                
                #Rule 3: the variant should be heterozygous in exactly one of the parents.
                if not father_num + mother_num == 1:
                    continue
                
                # Add v to the list of potential comphet variants of its associated gene(s).
                for g in annGenes(v):
                    chGene[g].append(v)
                    fromParent[g] |= {mother_num} # for rule 5 below (add 0 if paternal, 1 if maternal)
                
        # Rule 4, at gene level: at least two variants in the same gene.
        # Rule 5, at gene level: at least one variant from each parent.
        comphets = [v for gene,vars in chGene.iteritems() if len(vars)>1 and len(fromParent[gene])==2 for v in set(vars)]
        res = homoz + comphets            
    
    # Store resulting variants in a VCFtypeData object with the same attributes as VFch.
    VFres = VFch.copyAttributes(variants=res, appliedFilters=None)
    return VFres

    
def removeClosePairs(VF, minDist, variants_only):
    headers = VF.columnNames
    getChr = VF.chromGetter
    getPosStr = VF.posGetter

    def getPos(v): return int(getPosStr(v))

    res = []
    for i in [str(chri) for chri in range(1, 22)]+['X','Y']:
        chr_vars = sorted([v for v in VF.variants[:] if getChr(v) == i], key=getPos)
        L = len(chr_vars)
        if L == 0: continue
        positions = [getPos(v) for v in chr_vars]
        i = 0
        while(i < L-2):
            if positions[i + 1] - positions[i] > minDist:
                res.append(chr_vars[i])
                i += 1
            elif positions[i + 2] - positions[i] > minDist:
                i += 2
            else:
                i += 3
        if i == L-2:
            if positions[L-1] - positions[L-2] > minDist:
                res.extend(chr_vars[-2:])
        elif i == L-1:
            res.append(chr_vars[-1])

    return res if variants_only else VF.copyAttributes(variants=res)

    
if __name__ == "__main__":    
    reader = VariantFileReader.VariantFileReader()
    
    def test_comphetTrio():
        test = "C:\\Projects\\Filtus\\debug_ying.vcf"; ch_fa_mo=[0,1,2]
        input = dict(sep="\t", chromCol="CHROM", posCol="POS", geneCol="CLINVAR__GENEINFO_INFO", formatCol="FORMAT", splitAsInfo="INFO")
        vflist = reader.readVCFlike(test, keep00=1, **input)
        VFch, VFfa, VFmo = [vflist[i] for i in ch_fa_mo] 
        res = trioRecessiveFilter(VFch, VFfa, VFmo, "Recessive")
        for v in res.variants:
            print v
    
    def test_comphetTrio2():
        test = "C:\\Projects\\Filtus\\filtus.testRecessive.vcf"; ch_fa_mo=[0,2,1]
        input = dict(sep="\t", chromCol="CHROM", posCol="POS", geneCol="GENE_SYMBOL_INFO", formatCol="FORMAT", splitAsInfo="INFO")
        vflist = reader.readVCFlike(test, keep00=1, **input)
        VFch, VFfa, VFmo = [vflist[i] for i in ch_fa_mo] 
        res = trioRecessiveFilter(VFch, VFfa, VFmo, "Recessive")
        pos = [v[1] for v in res.variants]
        if not pos == ['22005167', '22017410']: 
            print "Error:\n", pos, str([22005167, 22017410])
            
    def test_regionfilter():
        #filt = Filter(regions_txt="C:/Projects/EXOME_ANALYSIS/Addison/Info/Design_Annotation_files/Target_Regions/SeqCap_EZ_Exome_v2_target_utenCHR.bed")
        filt = Filter(regions_txt="1:3-6")
        print filt
        #input = dict(sep="\t", chromCol="chr", posCol="pos", geneCol="")
        #vflist = reader.readVCFlike(test, keep00=0, **input)
        vf = reader.readNonVCF("C:\\Users\\magnusdv\\Desktop\\kast.txt", sep="\t", chromCol="chr", posCol="pos", geneCol="", gtCol="")
        vf2 = filt.apply(vf)
        print vf.length, vf2.length
            
    test_comphetTrio()
    test_comphetTrio2()
