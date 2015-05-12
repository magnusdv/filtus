import os.path
import re
import collections
from operator import itemgetter

import Filter
import FiltusUtils

#import memory_profiler
#import psutil

class ColumnData(object):
    def __init__(self, columnNames, variants, columnDescriptions=None, meta=''):
        self.columnNames = columnNames
        self.variants = variants if variants is not None else []
        self.columnDescriptions = {} if columnDescriptions is None else columnDescriptions
        self.meta = meta
        self._mainAttributes = ['columnNames', 'variants', 'columnDescriptions', 'meta']
        self.length = len(self.variants)

    def copyAttributes(self, **attrs):
        args = {atr : getattr(self, atr) for atr in self._mainAttributes}
        args.update(attrs)
        if not 'variants' in attrs: 
            args['variants'] = self.variants[:]
        if not 'columnNames' in attrs: 
            args['columnNames'] = self.columnNames[:]
        if "nGenes" in args and 'variants' in attrs:
            args['nGenes'] = None 
        return type(self)(**args)
        
    def columnGetter(self, *columns):
        colNames = self.columnNames
        if not columns: return None
        if any(col not in colNames for col in columns): 
            return None
        ind = [colNames.index(col) for col in columns]
        return itemgetter(*ind)
        
    def addData(self, newObj): #
        for i, a, b in zip(range(3), [self.varDefColNames, self.geneCol, self.gtCol],
                                   [newObj.varDefColNames, newObj.geneCol, newObj.gtCol]):
            if a!= b and False: #remove False...
                print '%s column: %s vs %s' %(['varDefColNames', 'geneCol', 'gtCol'][i], a, b)

        newvars = newObj.variants[:]
        heads, newheads = self.columnNames, newObj.columnNames
        commonheads = set([h for h in newheads if h in heads])
        
        def emptystring(x): return ''
        old_getters = [itemgetter(newheads.index(h)) if h in commonheads else emptystring for h in heads]

        newh = [h for h in newheads if not h in heads]
        if newh:
            getNew = itemgetter(*[newheads.index(h) for h in newh])
            if len(newh) == 1:
                newv = [tuple(getOld(v) for getOld in old_getters) + (getNew(v),) for v in newvars]
            else:
                newv = [tuple(getOld(v) for getOld in old_getters) + getNew(v) for v in newvars]
            addon = ('',)*len(newh)
            self.setVariants([v + addon for v in self.variants] + newv)
            self.columnNames.extend(newh)
        else:
            newv = [tuple(getOld(v) for getOld in old_getters) for v in newvars]
            self.setVariants(self.variants + newv)

    def setVariants(self, variants): #
        self.variants[:] = variants
        self.length = len(variants)
        if hasattr(self, 'getUniqueGenes'):
            self.nGenes = len(self.getUniqueGenes())

    def sort(self, column, descending=False, stringsFirst=False): 
        getCol = self.columnGetter(column)
        tryFloat = FiltusUtils.tryFloat
        stringval = float('inf')
        if descending: stringval = -stringval
        if stringsFirst: stringval = -stringval
        
        def getval(x):  
            return tryFloat(getCol(x), stringval=stringval)
            
        self.variants.sort(key = getval, reverse=descending)
        
        
    def _check(self, h, body):
        n = len(h)
        mi = min(len(v) for v in body)
        ma = max(len(v) for v in body)
        i = 0
        if mi < n or ma > n:
            print h
            print n, mi, ma
            for v in body:
                i += 1
                if len(v)!= n:
                    print i
                    print v
                    break

    def printData(self, trunc=None, pad=3):
        if self.length == 0:
            return '    '.join(self.columnNames), []
        heads = self.columnNames
        bodylist = self.variantStringList()
        #self._check(heads, bodylist) #include this to debug errors in mainDisplay
        if trunc is None: trunc = 50 #TODO: self.filtus.truncate
        if trunc > 0:
            _getWidth=self._getWidth
            wids = [_getWidth(i, heads, bodylist, trunc) for i in range(len(heads))]
            padded = [w + pad for w in wids]
            headerLine = ''.join('%-*s' % (pad, h) for pad, h in zip(padded, heads))
            body = '\n'.join(''.join('%-*.*s' % (pad, wid, entry) for pad, wid, entry in zip(padded, wids, line)) for line in bodylist)
        else:
            wids = [max(len(h), max(len(line[i]) for line in bodylist)) + pad for i, h in enumerate(heads)]
            headerLine = ''.join('%-*s' % (wid, h) for wid, h in zip(wids, heads))
            body = '\n'.join(''.join('%-*s' % (wid, entry) for wid, entry in zip(wids, line)) for line in bodylist)
        return headerLine, body

    def _getWidth(self, i, heads, bodylist, trunc):
        hlen = len(heads[i])
        if hlen >= trunc:
            return hlen
        elif any(len(line[i]) >= trunc for line in bodylist[:1000]):
            return trunc
        else:
            return max(max(len(line[i]) for line in bodylist[:1000]), hlen)

    def variantStringList(self):
        return [[str(elem) for elem in var] for var in self.variants]

    def save(self, outfilename, sep="\t", colnames=True, preamblePos=('Top',)): 
        meta = self.meta
        with open(outfilename, 'w') as utfil:
            if meta and 'Top' in preamblePos: 
                utfil.write(meta)
            if colnames and self.columnNames: 
                utfil.write(sep.join(self.columnNames) + '\n')
            utfil.write('\n'.join(sep.join(line) for line in self.variantStringList()) + '\n')
            if meta and 'Bottom' in preamblePos: 
                utfil.write(meta)

class GeneSharingResult(ColumnData):
    def __init__(self, data, nSamples, geneMaster, shareCounts, minSampleCount=1, meta=''):
        heads = ['Gene', 'SampleCount', 'SampleIndex', 'Variants', 'Unique', 'Length']
        if data and len(data[0]) == 8: heads += ['P_raw', 'P_bonf' ]
        
        columnDescriptions = {'Gene':'Gene name', 
            'SampleCount':'Number of affected samples whose remaining variants in this gene match the given model', \
            'SampleIndex':'Indices of the affected samples whose remaining variants in this gene match the given model', \
            'Variants':'Total number of variants in all affected samples matching the model criteria in this gene', \
            'Unique':'Number of unique variants matching the model criteria in this gene', \
            'Length':'Gene length (available if a file with gene lengths is indicated in Settings -> Gene lengths)'}
        ColumnData.__init__(self, columnNames=heads, variants=data, columnDescriptions=columnDescriptions, meta=meta)
        self.geneMaster = geneMaster
        self.shareCounts = shareCounts
        self.nSamples = nSamples
        self.minSampleCount = minSampleCount
        
    @classmethod
    def geneMaster(cls, geneMaster, nSamples, minSampleCount=1, genelengths={}, pValues=False, meta=''):
        intlist2string = FiltusUtils.intlist2string
        shareCounts, data = [0]*nSamples, []
        for gene in geneMaster.keys():
            geneData = geneMaster[gene]
            samplecount = geneData.nFiles()
            if samplecount < minSampleCount: 
                del geneMaster[gene]
                continue
            shareCounts[samplecount-1] += 1
            samples = intlist2string(geneData.getFiles())
            nvars = geneData.length
            nuniqvars = geneData.nUniqVars()
            length = genelengths[gene] if gene in genelengths else '-'
            _info = [gene, samplecount, samples, nvars, nuniqvars, length]
            data.append(_info)

        return cls(data, nSamples, geneMaster, shareCounts, minSampleCount=minSampleCount, meta=meta)
        
    def addPvalues(self, gene, genelengths, m_aver, length, totL, model): #TODO: fix this!
        if gene in genelengths:
            pval = pValue(m_aver, length/totL, nSamples, samplecount, model=model)
            pval_bonf = min(pval * M, 1)
            return ['{:.3g}'.format(pval), '{:.3g}'.format(pval_bonf)]
        else:
            return ['-', '-']
            
    def copyAttributes(self, data=None):
        if data is None: data=self.variants[:]
        return GeneSharingResult(data=data, nSamples=self.nSamples, geneMaster=self.geneMaster, 
               shareCounts=self.shareCounts, minSampleCount=self.minSampleCount, meta=self.meta)
        
    def variantsInGenes(self, genes=None):
        gD = self.geneMaster
        if not gD: 
            return ColumnData(columnNames=["No content"], variants=[], meta=self.meta)
        if genes is None: genes = list(gD)
        allVars = MultiFileData(parent=gD[genes[0]], variants = [], genes = [], meta=self.meta)
        for g in genes:
            allVars.addData(gD[g])
        allVars.genes = genes #TODO: avoid this by adding addData function in multifiledata class.
        return allVars  
        
    #def remove(self, other):
    #    gM, gM2 = self.geneMaster, other.geneMaster
    #    new = {g:v for g,v in gM.iteritems() if not (g in gM2 and len(gM[g]) == len(gM2[g]))}
    #    return GeneSharingResult.geneMaster(new, nSamples=nSamples, minSampleCount=self.minSampleCount, 
    #                genelengths={}, pValues=False, meta=''):
        
class NgData(ColumnData):
    def __init__(self, data, title, subtitle='', meta=''):
        ColumnData.__init__(self, columnNames=[], variants=data, columnDescriptions=None, meta=meta)
        bord = '='*max(len(title), len(subtitle))
        complete_title = title + subtitle if subtitle else title
        if subtitle:
            self.intro = '\n'.join([bord, title, subtitle, bord]) + '\n'
            analysis_text = '\n## '.join(["STEP-BY-STEP FILTERING TABLE:", title, subtitle])
        else:    
            self.intro = '\n'.join([bord, title, bord]) + '\n'
            analysis_text = '\n## '.join(["STEP-BY-STEP FILTERING TABLE:", title])
    
    def copyAttributes(self):
        return self
        
    def printData(self):
        wid1_all = [map(len, line[0].split('\n')) for line in self.variants]
        wid1_max = max(i for t in wid1_all for i in t)
        bodyaslist = [line[0] + ''.join([' ']*(wid1_max + 1 - wid1_all[i][-1])) + ''.join(s.rjust(8) for s in line[1:]) for i, line in enumerate(self.variants)]
        body = self.intro + '\n\n'.join(bodyaslist)
        header = ''
        return header, body

 
class VariantData(ColumnData):
    def __init__(self, filename, columnNames, variants, chromCol, posCol, geneCol, gtCol, homSymbol=None, 
                 columnDescriptions=None, appliedFilters=None, startupFilter=None, nGenes=None, meta=''):
        ColumnData.__init__(self, columnNames, variants, columnDescriptions, meta=meta)

        self.filename = filename
        self.longName = filename
        self.shortName = self.getShortName()
        self.chromCol = chromCol
        self.posCol = posCol
        self.geneCol = geneCol
        self.gtCol = gtCol
        self.homSymbol = homSymbol
        self.varDefColNames = [chromCol, posCol]
        self.isVCFtype = False # overridden to True in VCF
        self.keep00 = False # ad hoc. Overridden in VCF files.
        
        self.chromGetter   = self.columnGetter(chromCol)
        self.posGetter     = self.columnGetter(posCol)
        self.geneGetter    = self.columnGetter(geneCol)
        self.genotypeGetter = self.columnGetter(gtCol)
        self.varDefGetter  = self.columnGetter(*self.varDefColNames)
        self.chromPosRefAlt = None # used in VCF files
        
        self.appliedFilters = appliedFilters ## NB! Will be overwritten by other filters!
        if startupFilter is not None:
            if isinstance(startupFilter, basestring):
                startupFilter = Filter.Filter(filterFile=startupFilter)
            startupFilter.checks(self)
            startupFilter.apply(self, checks = False, inplace=True)
        self.nGenes = len(self.getUniqueGenes()) if nGenes is None else nGenes
        self._mainAttributes += ['filename', 'chromCol', 'posCol', 'geneCol', 'gtCol', 'homSymbol', 'nGenes']
        
        
    def getShortName(self):
        if self.filename is None: return None
        return os.path.splitext(os.path.basename(self.filename))[0]
    
    def index(self, longFileNameList):
        return longFileNameList.index(self.longName)
            
    def getVar(self, values, firstHit, columns=None):
        '''return first (or all) variant(s) whose (chrom,pos) equals the string pair chromPos'''
        if columns is None:
            getcols = self.varDefGetter if len(values)==2 else self.chromPosRefAlt if len(values)==4 else None
        elif len(values)==len(columns):
            getcols = self.columnGetter(*columns)
        else:
            getcols = None
        if getcols is None: return None
        if firstHit:
            return next((v for v in self.variants if getcols(v)==values), None)
        else:
            return [v for v in self.variants if getcols(v)==values]
        
    def confirmGTcolumn(self):
        if self.genotypeGetter is None:
            raise ValueError("Unknown genotype column in the following sample:\n\n%s" %self.longName)
    
    def confirmGeneColumn(self):
        if self.geneGetter is None:
            raise ValueError("Unknown gene column in the following sample:\n\n%s" %self.longName)
    
    def noHomozygotes(self):
        isHom = self.isHomALT()
        return isHom and not any(isHom(v) for v in self.variants)

    def pedrow(self, total_set_sorted): # input: sorted list of variants (chr, pos) from ALL samples. Out: vector of length 2 * len(total_set), containing the alleles (with 0's where self does not have the variant).
        chrom_getter = self.chromGetter
        pos_getter = self.posGetter
        
        res = [0]*2 * len(total_set_sorted)
        if self.length == 0:
            return res
        chromInt = FiltusUtils.chromInt
        def chrom_pos(v):
            return (chromInt(chrom_getter(v)), int(pos_getter(v)))
        
        GTnum = self.GTnum()
        GTdict = {0: (1,1), 1: (1,2), 2:(2,2)}
        def gt(v):
            return GTdict[GTnum(v)]

        variants_sorted_rev = sorted(self.variants, key=chrom_pos, reverse=True) #reversing to prepare for pop() below
        pos_reversed = [chrom_pos(v) for v in variants_sorted_rev]
        gt_reversed = [gt(v) for v in variants_sorted_rev]
        
        for i, data in enumerate(total_set_sorted):
            pos = data[:2]
            if pos == pos_reversed[-1]:
                pos_reversed.pop()
                res[2 * i:2*(i + 1)] = gt_reversed.pop()
                if not gt_reversed:
                    break
                if pos == pos_reversed[-1]: # ad hoc solution to repeated positions.
                    pos_reversed.pop()
                    gt_reversed.pop()
        return res

    def allChromPos(self, what="set", reverse=False):
        '''Returns a set or a sorted list of pairs (chrom, position) for all variants in the sample'''
        chrom, pos=self.chromGetter, self.posGetter
        chromInt = FiltusUtils.chromInt
        if what == "set":
            return set((chromInt(chrom(v)), int(pos(v))) for v in self.variants)
        if what == "list":
            return sorted(((chromInt(chrom(v)), int(pos(v))) for v in self.variants), reverse=reverse)
        
    def allDataSet(self, columns):
        '''Returns a set of unique combinations of the specified columns for all variants in the sample'''
        cols = self.columnGetter(*columns)
        if cols is None: 
            raise KeyError("Unknown columns selected")
        return set(cols(v) for v in self.variants)
        
    def nVars(self, diploid=False): #If diploid, homozygous variants are counted twice
        if diploid:
            nAlleles = self.GTnum()
            if nAlleles is None: return 
            return sum(nAlleles(v) for v in self.variants)
        else:
            return self.length # TODO: not correct with REFREFs

    def nUniqVars(self): #
        return len(self.getUniqueVariants())

    def isHomALT(self):
        ''' Returns a function determining whether a variant is ALT/ALT '''
        if not self.gtCol: return
        gt = self.columnGetter(self.gtCol)
        homSymbol = self.homSymbol
        def _f(v):
            return gt(v)==homSymbol
        return _f
        
    def GTnum(self, noREFREF=True):
        ''' Returns a function taking v --> 0 if GT=0/0, 1 if a/b (a!=b) and 2 if GT=c/c (c>0).'''
        if not self.gtCol: return
        isHom = self.isHomALT()
        def _f(v):
            return 2 if isHom(v) else 1
        return _f
        
    def collapse(self): #?
        collapsedList = []
        uniq = set()
        varDefGetter = self.varDefGetter
        for var in self.variants:
            v_string = varDefGetter(var)
            if v_string not in uniq:
                collapsedList.append(var)
                uniq.add(v_string)
        self.setVariants(collapsedList)

    def annotatedGenes(self, v): #
        g = self.geneGetter(v)
        if '(' in g:
            g = re.sub('\([^\)]*\)', '', g) #remove all substrings in parenthesis, e.g. A1(dist=2), A2(dist=3) --> A1, A2
        if ';' in g or ',' in g:
            g = set(re.split('[;,]', g))
        else: g = set([g])
        return g

    def getUniqueGenes(self): #
        if self.length == 0 or self.geneGetter is None: return set()
        annGenes = self.annotatedGenes
        return set.union(*[annGenes(v) for v in self.variants])

    def getUniqueGenes_homoz(self): # not in use
        isHom = self.isHomALT()
        annGenes = self.annotatedGenes
        if self.length == 0 or not isHom or not self.geneGetter: return set()
        return set.union(*[annGenes(v) for v in self.variants if isHom(v)])

    def getUniqueGenes_comphet(self): #Not in use
        isHom = self.isHomALT()
        if self.length == 0 or not isHom or not self.geneGetter: return set()
        return set(g for g, vars in self.geneDict().iteritems() if len(vars)>1 or any(isHom(v) for v in vars))

    def getUniqueVariants(self, alleles=None, homoz=False): #
        varDefGetter = self.varDefGetter
        if not homoz and alleles is None:
            if self.keep00:
                nAlleles = self.GTnum()
                return set(varDefGetter(v) for v in self.variants if nAlleles(v)>0)
            else:
                return set(varDefGetter(v) for v in self.variants)
        elif homoz:
            isHom = self.isHomALT()
            return set(varDefGetter(v) for v in self.variants if isHom(v))
        else:
            nAlleles = self.GTnum()
            return set(varDefGetter(v) for v in self.variants if nAlleles(v) == alleles)

    def geneVars(self, genes, index): #
        if self.geneGetter is None: return {}
        annGenes = self.annotatedGenes
        fileno = (str(index + 1), )
        res = [fileno + v for v in self.variants if any(gene in annGenes(v) for gene in genes)]
        return MultiFileData(parent=self, variants=res, genes=genes)

    def geneDict(self, addIndex=None): #
        if self.geneGetter is None: return {}
        geneDict = collections.defaultdict(list)
        annGenes = self.annotatedGenes
        if addIndex is not None: 
            fileno = (str(addIndex + 1), )
        GT = self.GTnum()
        for v in self.variants:
            if GT(v) == 0: continue
            newv = v if addIndex is None else fileno + v
            for gene in annGenes(v):
                geneDict[gene].append(newv)
        for gene, vars in geneDict.iteritems():
            geneDict[gene] = MultiFileData(parent=self, variants=vars[:], genes=gene)
        return geneDict



class VCFtypeData(VariantData):
    def __init__(self, filename, columnNames, columnDescriptions, variants, chromCol, posCol, geneCol, 
                 formatHeads, splitFormat, splitInfo,  startupFilter=None, appliedFilters=None, keep00=None, nGenes=None, meta=''):    
        gtCol = 'GT' if splitFormat else columnNames[-1]
        homSymbol = '<vcf format>'
        VariantData.__init__(self, filename, columnNames, variants, chromCol, posCol, geneCol, gtCol, homSymbol, 
                             columnDescriptions=columnDescriptions, startupFilter=startupFilter, appliedFilters=appliedFilters, nGenes=nGenes, meta=meta)
        self.splitFormat = splitFormat 
        self.splitInfo = splitInfo
        self.isVCFtype = True
        self.keep00 = keep00
        self.formatHeads = formatHeads
        try:
            self.refCol = next(h for h in columnNames if h.lower() in ['ref', 'vcf_ref', 'ref_vcf'])
            self.altCol = next(h for h in columnNames if h.lower() in ['alt', 'vcf_alt', 'Obs'])
            self.chromPosRefAlt = self.columnGetter(chromCol, posCol, self.refCol, self.altCol) 
        except:
            FiltusUtils.warningMessage("Unknown REF/ALT columns in %s" %filename)
            self.chromPosRefAlt = self.varDefGetter
        self._mainAttributes = [a for a in self._mainAttributes if not a in ['gtCol', 'homSymbol']] + ['formatHeads', 'splitFormat', 'splitInfo',  'keep00']
                     
                     
    def variantDetailsDict(self, vdef):
        res = collections.OrderedDict(Sample=self.shortName) # preserves order of elements, makes sorting the headers easier 
        v = self.getVar(vdef, firstHit=True)
        if not v: 
            return res
        colnames = self.columnNames
        filterCol = next((h for h in colnames if h.lower() in ['filter', 'vcf_filter']), None)
        if filterCol:
            res[filterCol] = v[colnames.index(filterCol)]
        res.update(zip(['REF','ALT'], self.chromPosRefAlt(v)[2:]))
        
        if self.splitFormat:
            for h in self.formatHeads: 
                if h in colnames: res[h] = v[colnames.index(h)]
        else:
            res.update(zip(v[-2].split(':'), v[-1].split(':')))
        return res
        
    def getShortName(self):
        if self.filename is None: return None
        base = os.path.basename(self.filename)
        if base.rfind('.') > base.rfind('_'):
            return os.path.splitext(base)[0]
        else:
            return base.split('_')[-1]
        
    
    def GTnum(self):
        ''' Returns 0 if GT=0/0, 1 if a/b (a!=b) and 2 if GT=c/c (c>0).'''
        gt = self.columnGetter(self.gtCol)
        def _f(v):
            g = gt(v)
            return 1 if g[0] != g[2] else 0 if g[0]=='0' else 2
        return _f
            
    def isHomALT(self):
        gt = self.columnGetter(self.gtCol)
        def _f(v):
            g = gt(v)
            return 0 != g[0] == g[2]
        return _f    
        
    

class MultiFileData(VariantData):
    def __init__(self, parent, variants, genes, meta=''):
        columnNames = parent.columnNames[:]
        if columnNames[0] != 'Sample': columnNames = ['Sample'] + columnNames
        if isinstance(genes, basestring): genes = [genes]
            
        VariantData.__init__(self, filename=parent.filename, columnNames=columnNames, variants=variants, 
            columnDescriptions=parent.columnDescriptions, chromCol=parent.chromCol, posCol=parent.posCol, 
            geneCol=parent.geneCol, gtCol=parent.gtCol, homSymbol=parent.homSymbol, nGenes=len(genes), meta=meta)
        
        if parent.isVCFtype:
            self.splitFormat = parent.splitFormat 
            self.splitInfo = parent.splitInfo
            self.isVCFtype = True
            self.keep00 = parent.keep00
            self.formatHeads = parent.formatHeads
            self.GTnum = self.GTnum_vcfFIX
            self.isHomALT = self.isHomALT_vcfFIX
        self.genes = genes
        
        
    def GTnum_vcfFIX(self):
        ''' Returns 0 if GT=0/0, 1 if a/b (a!=b) and 2 if GT=c/c (c>0).'''
        gt = self.columnGetter(self.gtCol)
        def _f(v):
            try:
                g = gt(v)
                return 1 if g[0] != g[2] else 0 if g[0]=='0' else 2
            except: print g; print v[-6:]
        return _f
            
    def isHomALT_vcfFIX(self):
        gt = self.columnGetter(self.gtCol)
        def _f(v):
            g = gt(v)
            return 0 != g[0] == g[2]
        return _f    
        
    def copyAttributes(self, variants=None):
        if variants is None:
            variants = self.variants[:]
        return MultiFileData(parent=self, variants=variants, genes=self.genes[:], meta=self.meta)
        
    def intersectData(self, sgd):
        varDefGetter = self.varDefGetter
        shared = self.getUniqueVariants() & sgd.getUniqueVariants()
        self.addData(sgd)
        self.setVariants([v for v in self.variants if varDefGetter(v) in shared])

    def collapse(self): #TODO?
        collapsedVars = []
        vUniq = set()
        fileDic = collections.defaultdict(set)
        varDefGetter = self.varDefGetter
        for var in self.variants:
            vdef = varDefGetter(var)
            fileDic[vdef].add(int(var[0]))
            if vdef not in vUniq:
                collapsedVars.append(var)
                vUniq.add(vdef)

        i2s = FiltusUtils.intlist2string
        res = [(i2s(fileDic[varDefGetter(v)]), ) + v[1:] for v in collapsedVars]
        return MultiFileData(self, res, self.genes, meta=self.meta)

    def getFiles(self):
        return set(v[0] for v in self.variants)

    def nFiles(self):
        return len(self.getFiles())
