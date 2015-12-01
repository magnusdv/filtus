import os.path
import csv
import collections
import itertools
import re
import math
import time
from operator import itemgetter
from subprocess import call

import Pmw

import Filter
import FiltusUtils
import DataContainer
import AutEx

def geneLookup(genes, VFlist):
    if all(VF.geneGetter is None for VF in VFlist):
        raise IndexError("None of the variant files have known gene columns")
    
    result = None
    for i, VF in enumerate(VFlist):
        if VF.geneGetter is None: continue
        if result is None: 
            result = VF.geneVars(genes, index=i)
        else:
            result.addData(VF.geneVars(genes, index=i))
    
    result.meta = FiltusUtils.preambleNY(VFlist=VFlist, analysis='GENE LOOKUP\n## Genes: %s' % ', '.join(genes))
    return result
    
def genotypeData(vdef, VFlist):
    VFlist = [VF for VF in VFlist if VF.isVCFtype] # for now: only those with split FORMAT
    details = [VF.variantDetailsDict(vdef) for VF in VFlist]
    headers = FiltusUtils.listUnique([f for dic in details for f in list(dic)]) 
    
    values = []
    for d in details:
        if d: values.append([d.get(h, '') for h in headers]) 
        else: values.append(['']*len(headers))
    
    descr = VFlist[0].columnDescriptions
    meta = FiltusUtils.preambleNY(VFlist=VFlist, analysis='SINGLE VARIANT - UNFILTERED DATA FROM ALL SAMPLES\n## Chromosome: %s\n## Position: %s' %vdef)
    result = DataContainer.ColumnData(columnNames=headers, variants=values, columnDescriptions=descr, meta=meta)
    return result
            
def saveAsPed2(VFlist, pedigree, sampleIndex, dat, map, freq, maptype, genmapfile, freqColumn, defaultFreq, dir, prefix, freqBuffer=0.01): #merlin = True writes datfile and uses cM in mapfile.
     
    nInd = len(VFlist)
    chromInt = FiltusUtils.chromInt
    
    data_set = set()
    for sample in sampleIndex:
        if sample is None: continue
        VF = VFlist[sample]
        if not (VF.chromGetter and VF.posGetter and VF.gtCol):
            raise IndexError("Unknown chromosome/position/genotype column for sample %s" %VF.shortName)
        if freq:
            if freqColumn in VF.columnNames:
                thisVF_data = VF.allDataSet([VF.chromCol, VF.posCol, freqColumn])
            else:
                thisVF_data = {(chr, pos, defaultFreq) for chr, pos in VF.allDataSet([VF.chromCol, VF.posCol])}
            data_set.update(thisVF_data)
        else:   
            data_set.update(VF.allDataSet([VF.chromCol, VF.posCol]))
    
    if freq:
        def _freq(x):
            try: 
                fr = float(x)
                if fr < freqBuffer: return freqBuffer
                if fr > 1-freqBuffer: return 1-freqBuffer
                return fr
            except: 
                return defaultFreq
        
        # reduce (in case of multiple freq info) and sort
        seen = set() 
        all_sorted = sorted((chromInt(chr), int(pos), _freq(fr)) for (chr, pos, fr) in data_set \
                            if chromInt(chr) < 24 and (chr, pos) not in seen and not seen.add((chr, pos)))
    else: 
        all_sorted = sorted((chromInt(chr), int(pos)) for (chr, pos) in data_set)
    
    nMark = len(all_sorted)
    if not nMark:
        raise ValueError("The selected samples contain no variants")
    
    pedigreeLines = pedigree.strip().split('\n')
    pedname, mapname, datname, freqname = [os.path.join(dir, prefix + s) for s in ('.ped', '.map', '.dat', '.freq')]
    
    with open(pedname, 'w') as pedfile:
        for sample, pedline in zip(sampleIndex, pedigreeLines):
            pedfile.write(pedline + ' ')
            if sample is not None:
                VF = VFlist[sample]
                genotypes = ' '.join(str(x) for x in VF.pedrow(total_set_sorted=all_sorted))
            else:
                genotypes = ' '.join('0'*(2*nMark))
            pedfile.write(genotypes)
            pedfile.write('\n')
    files = [pedname]
    
    
    if map:
        mapfile = open(mapname, 'w')
        files.append(mapname)
        if 'cm' in maptype:
            genmap = GeneticMap(mapfilename=genmapfile)
            phys2cm = genmap.phys2cm
        maptype_paste = ''.join(maptype)
    if dat:
        datfile = open(datname, 'w')
        files.append(datname)
        datfile.write('A disease\n')
    if freq:
        freqfile = open(freqname, 'w')
        files.append(freqname)
    
    try:
        for i, vardat in enumerate(all_sorted): # vardat depends on the freq flag: Either triples (chrom, pos, fr) or pairs (chrom, pos)
            chrom, pos = vardat[:2]
            marker = 'var%d' %(i+1)
            if map: 
                if maptype_paste=='cm': mapline = '%d\t%s\t%f\n' %(chrom, marker, phys2cm(chrom, pos))
                elif maptype_paste=='phys': mapline = '%d\t%s\t%d\n' %(chrom, marker, pos)
                else: mapline = '%d\t%s\t%d\t%f\n' %(chrom, marker, pos, phys2cm(chrom, pos))
                mapfile.write(mapline)
            if dat: datfile.write('M %s\n' %marker)
            if freq: freqfile.write('M %s\nF %f %f\n' %(marker, 1-vardat[2], vardat[2])) 
    finally:
        if map: mapfile.close()
        if dat: datfile.close()
        if freq: freqfile.close()
        
    return files
  
  
class GeneticMap(object):
    def __init__(self, mapfilename=None):
        '''Map file must have headings, and 3 columns: CHROM, MB, CM'''
        conv = FiltusUtils.convertType
        chromInt = FiltusUtils.chromInt
        self.has_map = bool(mapfilename)
        if self.has_map:
            with open(mapfilename, 'rU') as mapfil:
                header = mapfil.next()
                mapMatrix = [map(conv, line.split('\t')) for line in mapfil]
            
            chromGroups = itertools.groupby(mapMatrix, key=itemgetter(0))
            self.chromMaps = {chromInt(chr) : zip(*chrmap)[1:] for chr, chrmap in chromGroups}
    
    def phys2cm(self, chr, phys):
        if not self.has_map:
            try:
                return phys/1.0e6
            except TypeError:
                return [p/1.0e6 for p in phys]
        pwl = FiltusUtils.piecewise_linear
        chr = FiltusUtils.chromInt(chr)
        chrmap = self.chromMaps[chr]
        try:
            return pwl(chrmap, phys/1.0e6) 
        except TypeError:
            return [pwl(chrmap, p/1.0e6) for p in phys]

class AlleleFreq(object):
    def __init__(self, VF, defaultFreq=0.5, altFreqCol=None, MAFcolumns=None, minmax=None):
        self.VF = VF
        self.default = float(defaultFreq)
        self.minmax = minmax
        self.minmaxRestrict = self._doNothing if minmax is None else self._doMinmax(minmax)
        self.altFreqCol = altFreqCol
        self.MAFcolumns = MAFcolumns
        
        for c in ([altFreqCol] if altFreqCol else []) + (list(MAFcolumns) if MAFcolumns else []):
            if c not in VF.columnNames: raise IndexError("Non-existent column name: %s" % c)
        if altFreqCol:
            self.type = "ALT"
            self.freqGetter = self._freqGetterALT()
        elif MAFcolumns:
            self.type = "MAF"
            self.needsplit = len(MAFcolumns)==3
            self.freqGetter = self._freqGetterMAF(self.needsplit)
        else:
            self.type = "BASIC"
            self.freqGetter = self._freqGetterBASIC()
    
    def __call__(self, v):
        return self.freqGetter(v)
        
    def __str__(self):
        if self.type == "ALT": return 'ALT allele frequency: Column "%s"; missing entries set to %g' % (self.altFreqCol, self.default)
        elif self.type == "MAF" and self.needsplit: 
            return 'Minor allele frequency: Column "%s" (format = "A:0.123"); missing entries set to %g' % (self.MAFcolumns[0], self.default)
        elif self.type == "MAF" and not self.needsplit: 
            return 'Minor allele/minor allele frequency: Columns %s/%s; missing entries set to %g' % (self.MAFcolumns[0], self.MAFcolumns[1], self.default)
        else: return 'ALT allele frequency: Constant value %g' % self.default
        
    def _doNothing(self, x):
        return x
    
    def _doMinmax(self, minmax):
        min, max = minmax
        def _mm(x):
            if x>max: return max
            if x<min: return min
            return x
        return _mm
        
    def _freqGetterBASIC(self):
        default = self.default
        def _f(v): return default
        return _f  
    
    def _freqGetterALT(self):
        freqcol = self.VF.columnGetter(self.altFreqCol)
        minmaxRestrict = self.minmaxRestrict
        default = self.default
        def _f(v):
            try: return minmaxRestrict(float(freqcol(v)))
            except: return default
        return _f
    
    def _freqGetterMAF(self, needsplit):
        '''
        MAFcolumns should be either:
            a tuple of 3 column names: (MAF, REF, ALT), where MAF are as "MAF_FREQ:MAF_ALLELE" (ex: "A:0.34")
            a tuple of 4 column names: (MAF_FREQ, MAF_ALLELE, REF, ALT)
        The freqGetter returns MAF_freq if MAF_allele == ALT, or 1-MAF_freq if MAF_allele == REF.
        ''' 
        MAFdata = self.VF.columnGetter(*self.MAFcolumns)
        default = self.default
        minmaxRestrict = self.minmaxRestrict
        maf_fix = self._maf_fix
        
        if needsplit:
            def _f(v):
                try:
                    mafstring, ref, obs = MAFdata(v)
                    allele, maf = mafstring.split(":")
                    return minmaxRestrict(maf_fix(allele, maf, ref, obs))
                except: 
                    return default
        else:
            def _f(v):
                try:
                    return minmaxRestrict(maf_fix(*MAFdata(v)))
                except: 
                    return default
        return _f 
    
    def _maf_fix(self, allele, maf, ref, obs):
        maf = float(maf)
        if allele==obs: return maf
        elif allele==ref: return 1 - maf
        elif len(ref)>len(obs):
            return maf if allele=='-' else 1 - maf
        elif len(ref)<len(obs):
            return 1 - maf if allele=='-' else maf
        
  
def pairwiseSharing(seleci, filtus):
    n = len(seleci)
    VFlist = [filtus.filteredFiles[i] for i in seleci]
    uniques = [VF.getUniqueVariants() for VF in VFlist]
    
    res1 = [['']+[i + 1 for i in seleci]]
    for i in range(n):
        res1.append([seleci[i]+1] + ['']*i + [len(set.intersection(uniques[i], uniques[j])) for j in range(i, n)])
    emptyhead, body1 = DataContainer.ColumnData(variants=res1, columnNames=['']*(n + 1)).printData()

    matr = [[float(x) if x else 0.0 for x in resline[1:]] for resline in res1[1:]]
    res2 = [['']+[i + 1 for i in seleci]]
    for i in range(n):
        diag = matr[i][i]
        if diag != 0:
            res2.append([seleci[i]+1] + [round(matr[j][i]/diag * 100, 1) for j in range(i)] +['-'] + [round(matr[i][j]/diag * 100, 1) for j in range(i + 1, n)])
        else:
            res2.append([seleci[i]+1] + ['-']*n)
    emptyhead, body2 = DataContainer.ColumnData(variants=res2, columnNames=['']*(n + 1)).printData()

    text = filtus.text
    text.clearAll()
    text.tag_configure("bold", font = filtus.monobold)
    
    text.settext("Variant sharing matrix:\n\n" + body1)
    text.tag_add('bold', '3.0', '3.end')
    for i in range(4, 4 + n):
        text.tag_add('bold', '%d.0'%i, '%d.0 wordend'%i)

    text.insert('end', "\n\n\nVariant sharing in percentages.\n(Row x, column y) = % of variants in x also present in y.\n\n")
    line = int(text.index('end').split('.')[0])
    text.insert('end', body2)
    text.tag_add('bold', '%d.0'%(line - 1,), '%d.end'%(line-1,))
    for i in range(line, line + n):
        text.tag_add('bold', '%d.0'%i, '%d.0 wordend'%i)
    text.meta = FiltusUtils.preambleNY(VFlist=VFlist, VFindex=seleci, analysis="PAIRWISE VARIANT SHARING")

    
def merge(VFlist, collapse):
    meta = FiltusUtils.preambleNY(VFlist=VFlist, analysis="MERGED SAMPLES - UNIQUE VARIANTS" if collapse else "MERGED SAMPLES")
    resVF = VFlist[0].copyAttributes(meta=meta)
    for VF in VFlist[1:]:
        resVF.addData(VF)
    if collapse: 
        resVF.collapse()
    return resVF

def makeIntersectionGeneDict(VFlist, VFindex, recessive): # used in GeneSharing classes and NgTable
    geneDic = VFlist[0].geneDict(addIndex=VFindex[0])
    for VF, i in zip(VFlist[1:], VFindex[1:]):
        for gene, vars in VF.geneDict(addIndex=i).iteritems():
            if gene in geneDic:
                gd = geneDic[gene]
                gd.intersectData(vars)
                if gd.length == 0:
                    del geneDic[gene]
                elif recessive and len(gd.getUniqueVariants())==1:
                    GTnum = gd.GTnum()
                    if any(GTnum(v)<2 for v in gd.variants): 
                        del geneDic[gene]
    return geneDic
    
def makeGeneDict(VFlist, VFindex, family=False): # used in GeneSharing classes and NgTable
    geneDic = VFlist[0].geneDict(addIndex=VFindex[0])
    for VF, i in zip(VFlist[1:], VFindex[1:]):
        for gene, vars in VF.geneDict(addIndex=i).iteritems():
            if gene in geneDic:
                geneDic[gene].addData(vars)
            else:
                geneDic[gene] = vars
    return geneDic


class ColumnSummary(object):
    def __init__(self):
        self.categorical_limit = 25
        self._NA = ('', 'NA', '-', '.')

    def summarize(self, VFlist, column):
        N = len(VFlist)
        if N==0: return
        hasCol = [column in VF.columnNames for VF in VFlist]
        getCols = [VF.columnGetter(column) for VF in VFlist]

        stringvecs = [[getCol(v) for v in VF.variants] if column in VF.columnNames else None for VF, getCol in zip(VFlist, getCols)]

        limit = self.categorical_limit
        many = any(len(set(vec)) > limit for vec in stringvecs if vec) or \
                 len(set.union(*[set(vec) if vec else set() for vec in stringvecs])) > limit #needed tweak here to ensure union() gets an argument

        allheaders = ['Sample', 'Total']
        numerical = False
        if many:
            _NA = self._NA
            try:
                numvecs = [[float(x) for x in vec if x not in _NA] if vec is not None else None for vec in stringvecs]
                numerical = True
            except ValueError: pass

        if numerical:
            summaries = [self.numSummary(has, vec, numvec) for has, vec, numvec in zip(hasCol, stringvecs, numvecs)]
            headers = ['min', 'mean', 'max', '<missing>', '<nonmissing>']
            listlist = [[str(i+1), VF.length] + [dic[h] for h in headers] for i, VF, dic in zip(range(N), VFlist, summaries)]
            allheaders += headers
        else:
            summaries = [self.catSummary(has, vec, many) for has, vec in zip(hasCol, stringvecs)]
            if many:
                headers = ['<missing>', '<nonmissing>']
            else:
                headers = sorted(set.union(*[set(cnt.keys()) for cnt in summaries if cnt is not None]), key = FiltusUtils.convertType)
                headers.sort(key=lambda x: x in ['<missing>', '<nonmissing>'])  #put missing/nonmissing last
            allheaders += headers + ['%%(%s)'%(h,) for h in headers]

            listlist = []
            for i, VF, dic in zip(range(N), VFlist, summaries):
                info = [str(i+1), VF.length]
                if dic is None:
                    x = info + ['-' for h in headers]*2
                elif VF.length == 0:
                    x = info + [0 for h in headers]*2
                else:
                    percents = ['%.2f' %(dic[h]/float(VF.length)*100,) if h in dic else 0 for h in headers]
                    x = info + [dic[h] if h in dic else 0 for h in headers] + percents
                listlist.append(x)
        meta = FiltusUtils.preambleNY(VFlist=VFlist, analysis="SUMMARY OF COLUMN '%s'"% column)
        res = DataContainer.ColumnData(allheaders, variants = listlist, meta=meta)
        return res

    def numSummary(self, has_header, vec, numvec):
        if not has_header:
            return {'min':'-', 'max':'-', 'mean':'-', '<nonmissing>':'-', '<missing>':'-'}
        N = len(numvec)
        if N == 0:
            return {'min':'-', 'max':'-', 'mean':'-', '<nonmissing>':0, '<missing>':len(vec)}
        return {'min':min(numvec), 'max':max(numvec), 'mean':'%.2f'%(sum(numvec)/N,), '<nonmissing>':N, '<missing>':len(vec) - N}

    def catSummary(self, has_header, vec, many):
        if not has_header:
            return None
        if many:
            missing = sum(vec.count(na) for na in self._NA)
            return {'<missing>' : missing, '<nonmissing>' : len(vec)-missing}
        elif len(vec) == 0:
            return {'<missing>' : 0, '<nonmissing>' : 0}
        else:
            res = collections.Counter(vec)
            for na in self._NA:
                if na in res:
                    res['<missing>'] += res[na]
                    del res[na]
        return res


class GeneSharingComputer(object):
    def __init__(self):
        pass
        #self.genelengths = {}
        
    def analyze(self, VFcases, VFcontrols, model, family, VFcases_index=None, VFcontrols_index=None, minSampleCount=1, genelengths=None):
        intlist2string = FiltusUtils.intlist2string
        
        filter = Filter.Filter(model=model, controls=VFcontrols, benignPairs=None if family else True)
        VFcases_filt = [filter.apply(VF) for VF in VFcases]
        
        VFcases_index = list(VFcases_index) if VFcases_index else range(len(VFcases))
        if VFcontrols_index:
            VFcontrols_index = list(VFcontrols_index)
        else:
            VFcontrols_index = range(len(VFcases), len(VFcases)+len(VFcontrols)) if VFcontrols else []
        
        recessive = model=="Recessive"
        if not family:
            gD = makeGeneDict(VFcases_filt, VFindex=VFcases_index, family=False)
        else:
            gD = makeIntersectionGeneDict(VFcases_filt, VFindex=VFcases_index, recessive=recessive)
            
            if recessive and len(VFcontrols) > 0:
                benignPairs = Filter.Filter.extractBenignPairs(VFcontrols)
                comb = itertools.combinations
                for g in gD.keys():
                    SGDat = gD[g]
                    GTnum = SGDat.GTnum()
                    if SGDat.length==1 and GTnum(SGDat.variants[0])<2:
                        del gD[g]
                    elif g in benignPairs and benignPairs[g]:
                        varDef = SGDat.varDefGetter
                        heterovars = SGDat.getUniqueVariants(alleles=1)
                        nonBenign = set(frozenset(pair) for pair in comb(heterovars, 2)) - benignPairs[g]
                        remov = heterovars - set(vdef for pair in nonBenign for vdef in pair)
                        keep = [v for v in SGDat.variants if not varDef(v) in remov]
                        if len(keep)>1 or (keep and GTnum(keep[0])==2): 
                            SGDat.setVariants(keep)
                        else:
                            del gD[g]
        
        # shareCounts, datalist = [0]*n, []
        # for gene in gD.keys():
            # geneData = gD[gene]
            # samplecount = geneData.nFiles()
            # if samplecount < minSampleCount: 
                # del gD[gene]
                # continue
            # shareCounts[samplecount-1] += 1
            # samples = intlist2string(geneData.getFiles())
            # nvars = geneData.length
            # nuniqvars = geneData.nUniqVars()
            # length = genelengths[gene] if gene in genelengths else '-'
            # _info = [gene, samplecount, samples, nvars, nuniqvars, length]
            # if includePvals:
                # if gene in genelengths:
                    # pval = pValue(m_aver, length/totL, n, samplecount, model=model)
                    # pval_bonf = min(pval * M, 1)
                    # _info += ['{:.3g}'.format(pval), '{:.3g}'.format(pval_bonf)]
                # else:
                    # _info += ['-', '-']
            # datalist.append(_info)

        analys_txt = "GENE SHARING - %s\n## Model: %s\n## Cases: %s\n## Controls: %s" \
                    %('FAMILY' if family else 'CASE/CONTROL', model, intlist2string([i + 1 for i in VFcases_index]), \
                    intlist2string([i + 1 for i in VFcontrols_index]))
        if minSampleCount > 1: 
            analys_txt += "\n## Minimum number of affected: %d" % minSampleCount
        meta = FiltusUtils.preambleNY(VFlist=VFcases+VFcontrols, VFindex=VFcases_index+VFcontrols_index, analysis=analys_txt)
        result = DataContainer.GeneSharingResult.geneMaster(gD, nSamples=len(VFcases), minSampleCount=minSampleCount, 
                                                            genelengths=genelengths, model=model, meta=meta)
        result.sort(column='SampleCount', descending=True)
        result.sort(column='P_raw', ascending=True)
        return result
    
    # def readGenelengths(self, file):
        # file = os.path.normpath(file)
        # res = {}
        # if not os.path.isfile(file):
            # raise IOError("Gene length file does not exist: %s."% file)
        # try:
            # with open(file, 'rU') as f:
                # for line in f:
                    # try:
                        # gen, length = line.split('\t')[:2]
                        # res[gen] = float(length)
                    # except ValueError:
                        # continue
            # if not res:
                # raise IOError("Empty file: %s."%file)
        # except ValueError as e:
            # raise ValueError("Gene length file does not have correct format (tab-separated columns: Gene - length in bp).\n\n%s" % e)

class VariantSharingComputer(object):
    def __init__(self):
        pass

    def analyze(self, VF_inputs, VF_index_list):
        field_alleles = (None, 1, 2, 2, None) #corresponding to f0, f1, f2, f01, f12
        adic = {VF:field_alleles[k] for k, field in enumerate(VF_inputs) for VF in field}

        shared_vdefs = set.intersection(*[VF.getUniqueVariants(alleles=adic[VF]) for VF in VF_inputs[1]+VF_inputs[2]+VF_inputs[4]])
        if VF_inputs[0] or VF_inputs[3]:
            shared_vdefs.difference_update(*[VF.getUniqueVariants(alleles=adic[VF]) for VF in VF_inputs[0]+VF_inputs[3]])
        
        heads = (VF_inputs[1]+VF_inputs[2]+VF_inputs[4])[0].varDefColNames

        analys_txt = "VARIANT SHARING\n## " + '\n## '.join('%s alleles: %s' %(a, ', '.join(str(i + 1) for i in field))  for a, field in zip(['0', '1', '2', '0/1', '1/2'], VF_index_list))
        meta = FiltusUtils.preambleNY(VFlist=[VF for input in VF_inputs for VF in input], VFindex=[i for field in VF_index_list for i in field], analysis = analys_txt)
        
        result = DataContainer.ColumnData(heads, variants = list(shared_vdefs))
        result.sort(column=heads[1], descending=False)
        result.sort(column=heads[0], descending=False)
        return result
         

class DeNovoComputer(object):
    def __init__(self):
        pass
        
    def _transm(self, m):
        if not 0 < m < 1: raise ValueError("Prior mutation rate must be between 0 and 1")
        n = 1-m
        T = [[[n**2, 2*m*n, m**2], # AA*AA
        [.5*n, .5, .5*m], # AA*AB
        [n*m, 1-2*n*m, n*m]], # AA*BB
        
        [[.5*n, .5, .5*m], # AB*AA
        [.25, .5, .25], # AB*AB
        [.5*m, .5, .5*n]], # AB*BB
        
        [[n*m, 1-2*n*m, n*m], # BB*AA
        [.5*m, .5, .5*n], # BB*AB
        [m**2, 2*m*n, n**2]]] # BB*BB
        return T
    
    def _REFREFchecker(self, GTonly, includeMissing, GT=None, AD=None, minREFperc=None):
        _REFREF = ("0/0", "0|0", "./.") if includeMissing else ("0/0", "0|0")
        if GTonly:
            def _f(v): return GT(v) in _REFREF
        else:
            default = 100 if includeMissing else 0
            _ADperc = self._ADperc
            def _f(v): return (GT(v) in _REFREF) or (_ADperc(AD(v), alleleNum=0, default=default) > minREFperc)
        return _f
    
    def _ADperc(self, ADfield, alleleNum, default=0.0):
        '''Input: AD field (string), alleleNum (integer or list of integers). 
        Output: Percentage of reads with the specified allele numbers.'''
        try: 
            reads = map(int, ADfield.split(','))
            if isinstance(alleleNum, list):
                return sum(float(reads[num]) for num in alleleNum)/sum(reads)*100
            else:
                return float(reads[alleleNum])/sum(reads)*100
        except (ValueError, ZeroDivisionError):
            return default
        except Exception as e:
            print e
            return default
            
    def _incompatibleAllele(self, v_child, v_fa, v_mo, alleleGetter):
        '''
        Example 1: Child=2/2, Father=0/1, mother=1/2 --> 2
        Example 2: Child=1/2, Father=0/0, mother=0/0 --> [1,2]
        '''
        al_ch, al_fa, al_mo = alleleGetter(v_child), alleleGetter(v_fa), alleleGetter(v_mo)
        if al_ch in [al_fa, al_mo]:
            return None #if child's genotype equals either parent, then defined to be benign
        compatible = (al_ch[0] in al_fa and al_ch[1] in al_mo) or (al_ch[0] in al_mo and al_ch[1] in al_fa)
        if compatible:
            return None
        unseen = [x for x in sorted(set(al_ch)) if x not in al_fa and x not in al_mo]
        if len(unseen) == 1: 
            return int(unseen[0])
        elif len(unseen) == 2: 
            return [int(unseen[0]), int(unseen[1])] # coding as single numeric (assuming < 10 alleles!)
        
        if (al_ch[0] in al_fa and al_ch[1] not in al_mo) or (al_ch[0] in al_mo and al_ch[1] not in al_fa):
            return int(al_ch[1])
        return int(al_ch[0])
            
    def PL2postprob(self, PLch, PLfa, PLmo, TRmatrix, bfrq):    
        '''Only works for diallelic variants (easy to generalize, but need frequencies for all alleles)'''
        log10, pow = math.log10, math.pow
        try:
            glc, glf, glm = map(self.pl2loglik, (PLch, PLfa, PLmo))  #convert PL's to logliks
            if len(glc) > 3: raise RuntimeError("Only works for diallelic variants")
            loga, logb = log10(1-bfrq), log10(bfrq)
            logHW = (2*loga, log10(2)+loga+logb, 2*logb) #(1-b)^2,2*(1-b)*b, b^2)
            logliks = [glf[f] + glm[m] + glc[c] + logHW[f] + logHW[m] + log10(TRmatrix[f][m][c]) for f in range(3) for m in range(3) for c in range(3)]
            postprob = pow(10, logliks[1])/sum(pow(10, ll) for ll in logliks) #logliks[1] corresponds to (0,0,1), i.e. de novo.
        except Exception as e:
            postprob = None
        return postprob
        
    def pl2loglik(self, PL):
        '''converts PL string to true logliks of AA, AB, BB'''
        pl10 = [-float(a)/10 for a in PL.split(',')]
        logRdenom = math.log1p(sum(math.pow(10, x) for x in pl10 if x < 0))/math.log(10)
        return [y - logRdenom for y in pl10]
        
    def analyze(self, VFch, VFfa, VFmo, trioID, mut, defaultFreq, altFreqCol=None, MAFcolumns=None, threshold=None, GTonly=None, minALTchild=None, maxALTparent=None, includeMissing=False):
        #### Checking input data
        if not 0 < defaultFreq < 1: raise ValueError("Default allele frequency must be between 0 and 1")
        if threshold is not None and not 0 <= threshold <= 1: raise ValueError("Posterior probability threshold must be between 0 and 1")
        if minALTchild is not None and not 0 <= minALTchild <= 100: raise ValueError("Minimum child ALT percentage must be between 0 and 100")
        if maxALTparent is not None and not 0 <= maxALTparent <= 100: raise ValueError("Maximum parent ALT percentage must be between 0 and 100")
        if GTonly is True and (minALTchild is None or maxALTparent is None):
            raise RunTimeError("When 'GTonly' is False, both parameters 'minALTchild' and 'maxALTparent' must be specified (between 0 and 100).") 
        if VFch.chromPosRefAlt is None: raise ValueError("The input file format is not VCF-like.") 
        for col in ['GT', 'AD', 'PL']:
            for VF in [VFch, VFfa, VFmo]:
                if not col in VF.columnNames: 
                    raise ValueError("This analysis requires a '%s' column. This is missing from:\n\n%s" %(col, VF.longName))
        
        #### Setup
        if GTonly is None: GTonly = (minALTchild is None or maxALTparent is None)
        
        TRmatrix = self._transm(mut)
        freq = AlleleFreq(VFch, defaultFreq=defaultFreq, altFreqCol=altFreqCol, MAFcolumns=MAFcolumns, minmax=(0.001, 0.999))
        
        vardef = VFch.chromPosRefAlt 
        GT, AD, PL = [VFch.columnGetter(x) for x in ('GT', 'AD', 'PL')]
        _isREFREF = self._REFREFchecker(GTonly, includeMissing, GT=GT, AD=AD, minREFperc=(100 - maxALTparent if maxALTparent else None))
        
        item02 = itemgetter(0,2) # used to extract alleles from genotype, e.g. '0/1' -> ('0', '1')
        def _alleles(v): return sorted(item02(GT(v)))
        def _multiallelic(v): return ',' in vardef(v)[3]
        
        pow, log10 = math.pow, math.log10
        _ADperc, _incompatibleAllele = self._ADperc, self._incompatibleAllele
        PL2postprob = self.PL2postprob
        
        #### Actual computations
        # Variants where both parents are REF/REF - or variant is multiallelic
        fa00 = {vardef(v):v for v in VFfa.variants if _isREFREF(v) or _multiallelic(v)}
        mo00 = {vardef(v):v for v in VFmo.variants if _isREFREF(v) or _multiallelic(v)}
        both00 = set(fa00) & set(mo00)
        denovo = []
        
        for v in VFch.variants:
            vdef = vardef(v)
            if vdef not in both00 or _isREFREF(v): continue
            v_fa, v_mo = fa00[vdef], mo00[vdef]
            
            if ',' in vdef[3]: #multiallelic!
                alleleNum = _incompatibleAllele(v, v_fa, v_mo, alleleGetter=_alleles)
                if not alleleNum: continue
            else:
                alleleNum = 1
            ALTch = _ADperc(AD(v), alleleNum, default=0)
            if minALTchild and ALTch < minALTchild: continue
            ALTfa = _ADperc(AD(v_fa), alleleNum, default=0)
            ALTmo = _ADperc(AD(v_mo), alleleNum, default=0)
            if maxALTparent and (ALTfa > maxALTparent or ALTmo > maxALTparent): continue
            
            post = PL2postprob(PLch=PL(v), PLfa=PL(v_fa), PLmo=PL(v_mo), TRmatrix=TRmatrix, bfrq=freq(v))
            if post is not None:
                if threshold and post < threshold: continue
                post_txt = '%.4f'%post
            else:
                post_txt = '-'
            
            denovo.append((post_txt, '%.1f'% ALTch, '%.1f'% ALTfa, '%.1f'% ALTmo) + v)
        
        heads = ['P(de novo|data)', '%ALT child', '%ALT father', '%ALT mother'] + VFch.columnNames
        
        analys_txt = "DE NOVO\n## Child: %d\n## Father: %d\n## Mother: %d\n" % tuple(i+1 for i in trioID)
        analys_txt += "## Mutation rate: %g\n## %s" % (mut, str(freq))
        meta = FiltusUtils.preambleNY(VFlist=[VFch, VFfa, VFmo], VFindex=trioID, analysis = analys_txt)
        
        resultVF = VFch.copyAttributes(columnNames=heads, variants = denovo, filename=None, meta=meta)
        resultVF.sort(column=heads[0], descending=True)
        return resultVF
        
   
    #def PL(liks):
    #    '''Convert true likelihoods (sum 1) to normalized phred-scaled PL values'''
    #    x,y,z = [-10*math.log10(y) for y in liks]
    #    mm = min(x,y,z)
    #    return (x-mm, y-mm,z-mm)

    
class NgTable(object):
    #TODO: This needs cleaning up. Move GUI stuff out of here.
    def __init__(self, sharingPage, count_what): #if what = 'genes': count genes. If what = 'variants', count variants
        self.sharingPage = sharingPage
        self.filtus = sharingPage.filtus
        self.count_what = count_what
        self.nFilters = None
        self.filter_order = None
        self.title = "%s counts after step-wise filtering." % ('Gene' if count_what == 'genes' else 'Variant',)
        self.prompt = Pmw.PromptDialog(self.filtus.parent, title = 'Sharing table setup', defaultbutton = 0,
                        buttons = ('OK', 'Cancel'), command=self.execute, entryfield_labelpos='n', label_justify = 'left')
        self.prompt.withdraw()

    def validateFilter(self):
        VFlist = self.filtus.files
        page = self.sharingPage
        if self.count_what == 'genes':
            files_missing_genecol = [VFlist[i].filenames for i in page.getCases() if VFlist[i].geneGetter is None]
            if files_missing_genecol:
                FiltusUtils.warningMessage("Missing gene column in:\n\n%s" %'\n'.join(files_missing_genecol))
                return False

    def table(self):
        page = self.sharingPage
        VFall = self.filtus.files
        try:
            entries = page.validateEntries(VFall)
        except ValueError as e:
            FiltusUtils.warningMessage(e)
            return
        complete_filter = self.filtus.FM.getFilter()
        filterlist, filtertextlist, prompt_text = complete_filter.prepareTablePrompt()

        if self.filter_order is None:
            self.filter_order = ', '+', '.join(map(str, range(1, len(filterlist) + 1)))

        self.nFilters = len(filterlist)
        self.prompt.configure(label_text = prompt_text)
        self.prompt.component('entryfield').setentry(self.filter_order)

        userinput = FiltusUtils.activateInCenter(self.filtus.parent, self.prompt)
        if not userinput:
            return

        f_steps, indiv_step, clean_start = userinput

        filterSteps = [(Filter.Filter(), 'No filters')] if clean_start else []
        for step in f_steps:
            step_filter = [filterlist[i-1] for i in step]
            step_text = '\n'.join(filtertextlist[i-1] for i in step)
            filter_args = {argname:filter for argname, filter in step_filter if argname != 'columnfilters'}
            cfs = [filter for argname, filter in step_filter if argname == 'columnfilters']
            if cfs:
                filter_args.update({'columnfilters':cfs})
            filterSteps.append((Filter.Filter(**filter_args), step_text))

        if self.count_what == 'genes':
            model, cases, controls = entries['model'], entries['VFcases_index'], entries['VFcontrols_index']
            if model == 'Dominant': subtitle = 'Dominant model.'
            elif model == 'Recessive': subtitle = 'Recessive model (homozygous or compound heterozygous).'
            else: subtitle = 'Recessive model (homozygous only.)'
            coldat = self.computeGeneTable(VFall, cases, controls, model, filterSteps, indiv_step)
        elif self.count_what == 'variants':
            subtitle = "The column of individual 'k' shows the numbers of remaining variants after including this individual."
            fields = [self.sharingPage.getentry(field=k) for k in range(5)] #list of VF's for each of the 5 allele fields 0, 1, 2, 01, 12
            coldat = self.computeVariantTable(VFall, fields, filterSteps, indiv_step)

        data = DataContainer.NgData(coldat, self.title, subtitle)
        meta = FiltusUtils.preambleNY(VFlist=VFall, analysis=data.analysis_text)
        self.filtus.text.prettyPrint(data, meta=meta, label="Filter table")

    def execute(self, button):
        if button is None or button == 'Cancel':
            self.prompt.deactivate(None)
        else:
            order = self.prompt.get()
            if bool(re.search('[^0-9,&* ]', order)) or ('*' in order and '***' not in order):
                FiltusUtils.warningMessage("Syntax error")
                return

            if '***' in order:
                order, indiv_step = order.split('***')
            else:
                indiv_step = 1

            steps = [map(int, s.strip().split('&')) for s in order.strip().split(',') if s.strip()]
            if not all(i>0 and i<= self.nFilters for s in steps for i in s):
                FiltusUtils.warningMessage("Non-existent filter number")
                return

            self.filter_order = order
            clean_start = order.strip().startswith(',')
            self.prompt.deactivate((steps, indiv_step, clean_start))

    def computeGeneTable(self, VFlist, cases, controls, model, filterSteps, indiv_step):
        n = len(cases)
        ind_order = range(0, n, int(indiv_step))
        if ind_order[-1] != n-1: ind_order.append(n-1)

        coldat = [["Shared by at least:"] + [str(i + 1) for i in ind_order]]
        for filter, stepText in filterSteps:
            VFlist = [filter.apply(VF, checks=False) for VF in VFlist]
            filter2 = Filter.Filter(model=model, controls = [VFlist[i] for i in controls], benignPairs=True)
            affVFlist = [filter2.apply(VFlist[i], checks=False) for i in cases]
            Ngrow = self._Ngrow(makeGeneDict(affVFlist, VFindex=range(len(affVFlist))), ind_order)
            coldat.append([stepText] + [str(i) for i in Ngrow])
        return coldat

    def _Ngrow(self, dicDict, ind_order):
        shareCounts = [0]*(ind_order[-1] + 1) # in other words: [0]*len(cases)
        for dic in dicDict.itervalues():
            shareCounts[dic.nFiles()-1] += 1
        return [sum(shareCounts[i:]) for i in ind_order]

    def computeVariantTable(self, VFlist, fields, filterSteps, indiv_step):
        # TODO: indiv_step??
        field_alleles = (None, 1, 2, 2, None) #corresponding to f0, f1, f2, f01, f12
        adic = {i:field_alleles[k] for k, field in enumerate(fields) for i in field}

        carriers = sorted(fields[1])+sorted(fields[2])+sorted(fields[4])
        remove = sorted(fields[0])+sorted(fields[3])

        def countVariants(VFlist):
                first = carriers[0]
                variantset = VFlist[first].getUniqueVariants(alleles=adic[first])
                resvec = [len(variantset)]
                for i in carriers[1:]:
                    variantset.intersection_update(VFlist[i].getUniqueVariants(alleles=adic[i]))
                    resvec.append(len(variantset))
                for i in remove: #fjerne alle varianter i f0-individer, og alle homozygote i f01-inidivder.
                    variantset.difference_update(VFlist[i].getUniqueVariants(alleles=adic[i]))
                    resvec.append(len(variantset))
                return resvec

        coldat = [["Individual:"] + [str(i + 1) for i in carriers + remove]]
        for filter, stepText in filterSteps:
            VFlist = [filter.apply(VF) for VF in VFlist]
            Ngrow = countVariants(VFlist)
            coldat.append([stepText] + [str(i) for i in Ngrow])
        return coldat


class PlinkComputer(object):
    #def __init__(self): pass
        
    def runPlink(self, VF, pedprefix="filtus2plink", outprefix='filtus2plink', dir='', verbose=True, strictROH = False, **args):
        params = dict(file=pedprefix, out=outprefix, noweb='')
        if strictROH:
            strict = {'homozyg-window-kb' : 0, 'homozyg-window-snp' : 1, 'homozyg-window-het' : 0,
                      'homozyg-window-missing' : 0, 'homozyg-window-threshold' : 1, 'homozyg-snp' : 2,
                      'homozyg-kb' : 0, 'homozyg-density' : 10000, 'homozyg-gap' : 1000000}
            params.update(strict)
        params.update(args) # overwrites strict params
        
        param_strings = ['--%s %s'%(key,val) for key, val in params.iteritems()]
        plinkcommand = 'plink ' + ' '.join(param_strings)
        if verbose: 
            print 'PLINK command:\n' + plinkcommand
        try:
            files = self.saveAsPedSingle(VF, dir=dir, prefix=pedprefix)
            if verbose:
                print "Created files:\n%s\n%s"%tuple(files)
        except Exception as e:
            raise IOError('Failed writing ped/map files:\n%s' %e)
        try:
            call(plinkcommand)
        except Exception as e:
            raise OSError('PLINK execution failed!\n\nError message:\n%s' %e)
        
        with open(outprefix+'.hom', "rU") as homfile:
            heads = [h.strip() for h in homfile.next().split(' ') if h]
            goodcols = itemgetter(*[heads.index(H) for H in heads if H not in ['FID', 'PHE', 'SNP1', 'SNP2']])
            heads = goodcols(heads)
            homReader = csv.reader(homfile, delimiter = ' ', skipinitialspace=True)
            homlist = [tuple(goodcols(line)) for line in homReader]
        return DataContainer.ColumnData(heads, variants = homlist)
            
    def saveAsPedSingle(self, VF, dir, prefix):
        if VF.length == 0: 
            raise IndexError("The selected samples contain no variants")
        if not VF.chromGetter and VF.posGetter and VF.gtCol:
            raise KeyError("Unknown chromosome/position/genotype column for some of the selected samples")
        pedpath = os.path.abspath(os.path.join(dir, prefix + '.ped'))
        mappath = os.path.abspath(os.path.join(dir, prefix + '.map'))

        all_sorted = sorted(VF.allChromPos())
        with open(pedpath, 'w') as pedfile:
            pedfile.write('1 1 0 0 1 1 ')
            pedfile.write(' '.join(str(x) for x in VF.pedrow(total_set_sorted=all_sorted)))
            pedfile.write('\n')
        
        with open(mappath, 'w') as mapfile:
            # format (X, var10, 0, 1234567) standard plink. Could have used merlin format, but require extra parameters.
            mapfile.write('\n'.join('%s\tvar%d\t0\t%d' %(str(chrom), i + 1, pos) for i, (chrom, pos) in enumerate(all_sorted)))
        
        return [pedpath, mappath]
        
    def summary(self, segmentsColDat):
        n = segmentsColDat.length
        if n==0:
            totMB, maxMB, fraction = [0]*3
        else:
            mb = [float(s[4])/1000 for s in segmentsColDat.variants]
            totMB = sum(mb)
            maxMB = max(mb)
            fraction = totMB/3000.0
        result = {"Segments":n, "Fraction":fraction, "Total (MB)":totMB,"Longest (MB)":maxMB, }
        return result

        
class AutExComputer(object):
    def __init__(self, genmapfile):
        self.genmap = GeneticMap(mapfilename=genmapfile)
        self.genmapfile = genmapfile
        
    def ranges(self, boo):
        ''' Makes a list of index ranges [start, stop] where boo is True.''' 
        res = []
        prev = 0
        for i,b in enumerate(boo):
            if prev == b: pass
            elif b: start = i
            else: res.append([start, i-1])
            prev = b
        if prev:
            res.append([start, i])
        return res
    
    def singleChrom_scores(self, VF, chrom, f, a, error, defaultFreq, altFreqCol=None, MAFcolumns=None):
        chr, pos, gtnum = VF.chromGetter, VF.posGetter, VF.GTnum()
        chrom = str(chrom)
        freq = AlleleFreq(VF, defaultFreq=defaultFreq, altFreqCol=altFreqCol, MAFcolumns=MAFcolumns, minmax=(0.01, 0.99))
        chromDataSorted = sorted((int(pos(v)), gtnum(v), freq(v)) for v in VF.variants if chr(v) == chrom)
        pos_phys, obs, freqs = zip(*chromDataSorted)
        pos_cm = self.genmap.phys2cm(chrom, pos_phys)
        scores = AutEx.FwdBwd(obs, pos_cm, freqs, f=f, a=a, error=error)
        return pos_phys, obs, freqs, pos_cm, scores
    
    def autex_segments(self, VF, f, a, error, defaultFreq, altFreqCol=None, MAFcolumns=None, threshold=0.5, minlength=0.0, unit='cM*', mincount=0, overrule_count=100):
        chr, pos, gtnum = VF.chromGetter, VF.posGetter, VF.GTnum()
        freq = AlleleFreq(VF, defaultFreq=defaultFreq, altFreqCol=altFreqCol, MAFcolumns=MAFcolumns, minmax=(0.01, 0.99))
        sortedVariants = sorted(VF.variants, key=VF.chromGetter)
        chromGroups = itertools.groupby(sortedVariants, key=chr)
        segments = []
        for chr, vars in chromGroups:
            if FiltusUtils.chromInt(chr) > 22: continue
            chromDataSorted = sorted((int(pos(v)), gtnum(v), freq(v)) for v in vars)
            pos_phys, obs, freqs = zip(*chromDataSorted)
            
            pos_cM  = self.genmap.phys2cm(chr, pos_phys)
            scores = AutEx.FwdBwd(obs, pos_cM, freqs, error=error, a=a, f=f)
            segs = self.scores2segments(scores, pos_phys=pos_phys, pos_cM=pos_cM, threshold=threshold, minlength=minlength, unit=unit, 
                                        mincount=mincount, overrule_count=overrule_count)
            chrres = [[chr] + seg for seg in segs]
            segments.extend(chrres)
        
        meta = FiltusUtils.preambleNY(VFlist=[VF], analysis = self.meta_string(f,a,error,freq,threshold=threshold, minlength=minlength, unit=unit, mincount=mincount, overrule_count=overrule_count))
        res = DataContainer.ColumnData(variants=segments, columnNames=['CHR', 'FROM*', 'TO*', 'MB*', 'CM*', 'FROM', 'TO', 'MB', 'CM', 'N'], meta=meta)
        res.sort(column="FROM*")
        res.sort(column="CHR")
        return res
    
    def meta_string(self, f, a, error, freq, threshold, minlength, unit, mincount, overrule_count=100):
        txt = "AUTOZYGOSITY MAPPING\n## Input parameters:\n##    f = %g\n##    a = %g\n## %s\n" % (f,a, str(freq))
        txt += "## Recombination map: %s\n" %("Uniform map with 1 cM = 1 Mb" if not self.genmapfile else "DecodeMap_thin.txt (adapted from Kong et al. (2010)).", )
        txt += "## Posterior probability threshold: %s\n" % threshold
        txt += "## Minimum region size: %s %s and %s variants" %(minlength, unit, mincount)
        if overrule_count:
            txt += "\n## (But always include regions with %d or more variants)" % overrule_count
        return txt
        
    def scores2segments(self, scores, pos_phys, pos_cM, threshold, minlength=0.0, mincount=0, unit='cM', overrule_count=100):
        '''converting vector of scores to list of segments of the form [start_ext, end_ext, mb_ext, cm_ext, start, end, mb, cm, count]'''
        rr = self.ranges(p >= threshold for p in scores)
        if not rr: 
            return []
        pos_physEXT = list(pos_phys) + [pos_phys[-1], 0] # so that (1) pos_phys[-1]=0 and (2) end+1 don't give indexError
        pos_cmEXT = list(pos_cM) + [pos_cM[-1], 0] #
        
        digits = 2
        spans = [(pos_physEXT[start-1], pos_physEXT[end+1], pos_phys[start], pos_phys[end], end-start+1, pos_cmEXT[start-1], pos_cmEXT[end+1], pos_cM[start], pos_cM[end]) for start, end in rr]
        segs = [[s[0], s[1], round((s[1]-s[0])/1.0e6, digits), round(s[6]-s[5], digits), 
                 s[2], s[3], round((s[3]-s[2])/1.0e6, digits), round(s[8]-s[7], digits), s[4]] for s in spans]
        
        useLengthColumn = {'Mb*':2, 'cM*':3, 'Mb':6, 'cM':7}[unit]
        segs = [s for s in segs if s[-1]>=overrule_count or (s[useLengthColumn] >= minlength and s[-1] >= mincount)]
        return segs
            
    def summary(self, segmentsColDat):
        n = segmentsColDat.length
        if n==0:
            totMB_ext, totCM_ext, maxMB_ext, maxCM_ext, totMB, totCM, maxMB, maxCM, fraction_ext, fraction = [0]*10
        else:
            mb_ext = [s[3] for s in segmentsColDat.variants]; totMB_ext = sum(mb_ext); maxMB_ext = max(mb_ext)
            cm_ext = [s[4] for s in segmentsColDat.variants]; totCM_ext = sum(cm_ext); maxCM_ext = max(cm_ext)
            mb = [s[7] for s in segmentsColDat.variants]; totMB = sum(mb); maxMB = max(mb)
            cm = [s[8] for s in segmentsColDat.variants]; totCM = sum(cm); maxCM = max(cm)
            fraction_ext = totCM_ext/3200.0
            fraction = totCM/3200.0
        result = {"Number of segments":n, "Fraction*":fraction_ext, "Fraction":fraction, "Total MB*":totMB_ext, "Total CM*":totCM_ext,
                      "Total MB":totMB, "Total CM":totCM, "Longest MB*":maxMB_ext, "Longest CM*":maxCM_ext, "Longest MB":maxMB, "Longest CM":maxCM}
        return result
        
       

        
def _printVF(VF, truncate=10):        
    pretty = VF.printData(trunc=truncate)
    print pretty[0]
    print pretty[1]
   
   
if __name__ == "__main__":
    
    import VariantFileReader
    reader = VariantFileReader.VariantFileReader()
    
    def test_denovo():
        test = "example_files\\multiallelic_tests.vcf"; frqCol=""; ch_fa_mo=[0,1,2]
        vflist = reader.readVCFlike(test, sep="\t", chromCol="CHROM", posCol="POS", geneCol="", splitAsInfo="INFO", keep00=1)
        VFch, VFfa, VFmo = [vflist[i] for i in ch_fa_mo] 
        dn = DeNovoComputer()
        res1 = dn.analyze(VFch, VFfa, VFmo, trioID=ch_fa_mo, mut=1e-8, defaultFreq=.1, altFreqCol=frqCol, GTonly=None, minALTchild=30, maxALTparent=10, includeMissing=0)
        _printVF(res1)
        res2 = dn.analyze(VFch, VFfa, VFmo, trioID=ch_fa_mo, mut=1e-8, defaultFreq=.1, altFreqCol=frqCol, GTonly=1, minALTchild=30, maxALTparent=10, includeMissing=1)
        _printVF(res2)
        
    test_denovo()
    