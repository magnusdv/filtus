import FiltusUtils
import DataContainer
import csv
import re

class VariantFileReader(object):
    def __init__(self):
        pass
        
    def readNonVCF(self, filename, sep, splitAsInfo=None, split_general=None, skiplines=0, **params): # gtCol, homSymbol:
        with open(filename, "rU") as ifile:
            for i in range(skiplines): ifile.next()
            reader = csv.reader(ifile, delimiter = sep, skipinitialspace = True, strict=True)
            headers = reader.next()
            if splitAsInfo or split_general:
                variants = [row for row in reader]
            else:
                variants = [tuple(row) for row in reader]
        
        if splitAsInfo or split_general: # this is allowed also in nonVCF files
            if splitAsInfo:
                infoInd = headers.index(splitAsInfo)
                headers[:], variants[:] = self._splitINFO(headers, variants, infoInd)
            for col, sep in split_general:  
                headers[:], variants[:] = self._splitGeneral(headers, variants, col, sep)        
            variants[:] = map(tuple, variants)    
        
        if len(variants)>0 and len(headers) > len(variants[0]):
            raise csv.Error("Format error: The number of headings is greater than the number of entries in the first row.")
        
        ## Nijmegen fix: If a blank line appears, delete it everything after it. (For speed: Only last 200 lines checked)
        blank_tail = [all(len(entry.strip()) == 0 for entry in line) for line in reversed(variants[-200:])]
        if any(blank_tail):
            variants[:] = variants[:-blank_tail.index(True)-1]

        #create instance
        return DataContainer.VariantData(filename=filename, columnNames=headers, variants=variants, **params)

    def readVCFlike(self, filename, sep, chromCol, posCol, geneCol, keep00=0, split_general=[],
                    splitAsInfo="", formatCol="FORMAT", splitFormat=1, commentChar="##", startupFilter=None):
        preamble = []
        with open(filename, "rU") as ifile:
            line = ifile.next()
            while commentChar and line.startswith(commentChar):
                preamble.append(line)
                line = ifile.next()
            headers = [h.strip() for h in line.split(sep)]
            
            data = [row for row in csv.reader(ifile, delimiter = sep, skipinitialspace = False, strict=True)] # saving time (?) with skip.. = False
        descriptions = self._parseDescriptions(preamble)
        
        if headers[-1] == "Otherinfo": ### Annovar fix: Re-inserting VCF column names.
            headers = self._fixAnnovarOtherinfo(headers, firstvar=data[0])
        if headers[0] == '#CHROM': ### VCF tweak
            headers[0] = 'CHROM'
        if splitAsInfo:
            infoInd = headers.index(splitAsInfo)
            headers[:], data[:] = self._splitINFO(headers, data, infoInd)
        
        split_general_common = [(x,y) for x,y in split_general if x in headers]
        split_general_special = [(x,y) for x,y in split_general if x not in headers]
        
        for col, sep in split_general_common:  
                headers[:], data[:] = self._splitGeneral(headers, data, col, sep)
        
        formatIndex = headers.index(formatCol)    
        if splitFormat:
            formatHeads, allEqual = self._getFormatHeads(data, formatIndex)
        else: 
            formatHeads, allEqual = None, None
            
        if keep00:
            def missing(val): return len(val) < 3 or (val[0]==val[2]=='.')
        else:
            def missing(val): return len(val) < 3 or (val[2] in ['0', '.'] and val[0] in ['0', '.'])
        
        paramsVCF = dict(chromCol=chromCol, posCol=posCol, geneCol=geneCol, startupFilter=startupFilter,
                 splitFormat=splitFormat, splitInfo=bool(splitAsInfo), keep00=keep00, formatHeads=formatHeads)
        VFlist = []
        Nsamples = len(data[0]) - (formatIndex+1) 
        for i in range(formatIndex+1, len(data[0])):
            samplename = headers[i]
            fname = filename + '_' + samplename if samplename!="<sample>" else filename
            if splitFormat:
                h, variants = self._splitFORMAT(headers[:], data, formatIndex, i, formatHeads, allEqual, missingFunc=missing)
                for col, sep in split_general_special:  
                    h[:], variants[:] = self._splitGeneral(h, variants, col, sep)
                    
                variants[:] = map(tuple, variants)
            else:
                common = formatIndex + 1
                variants = [tuple(x[:common] + [x[i]]) for x in data if not missing(x[i])]
                h = headers[:common] + [samplename]
            VFlist.append(DataContainer.VCFtypeData(fname, h, descriptions, variants, **paramsVCF))
            #print 'sample finished: %s'%str(time.time()-st); st=time.time()
        
        return VFlist

    def _getFormatHeads(self, data, formatIndex):
        formats = set(x[formatIndex] for x in data)
        allEqual = len(formats) == 1
        if allEqual:
            formatHeads = list(formats)[0].split(':')
        else:
            formatHeads = FiltusUtils.listUnique([field for F in formats for field in F.split(':')])
        return formatHeads, allEqual

    def _splitFORMAT(self, h, data, formatIndex, sampleIndex, formatHeads, allEqual, missingFunc):
        L = formatIndex + len(formatHeads)

        if allEqual:
            variants = [x[:formatIndex] + x[sampleIndex].split(':') for x in data if not missingFunc(x[sampleIndex])]
            [x.extend(['']*(L-len(x))) for x in variants if len(x) < L]
        else:
            variants = []
            for x in data:
                if missingFunc(x[sampleIndex]): continue
                dic = dict(zip(x[formatIndex].split(':'), x[sampleIndex].split(':')))
                variants.append(x[:formatIndex] + [dic.get(tag, '') for tag in formatHeads])

        h[formatIndex:] = formatHeads
        return h, variants


    def _splitINFO(self, h, variants, infoInd):
        info = [dict(s.split('=') if '=' in s else [s,'1'] for s in r[infoInd].split(';')) for r in variants]
        all_tagheads = sorted(set(tag for dic in info for tag in dic))
        h[infoInd:(infoInd + 1)] = [tag + '_INFO' for tag in all_tagheads]
        for x, infodic in zip(variants, info):
            x[infoInd:(infoInd + 1)] = [infodic.get(tag, '') for tag in all_tagheads]
        return h, variants

    def _splitGeneral(self, h, variants, splitCol, sep):
        if splitCol not in h:
            err = "Something is wrong: Column '%s' to be split, is not among the current column names:\n\n%s" %(splitCol, ', '.join(h))
            raise RuntimeError(err)
        ind = h.index(splitCol)
        counts = set(v[ind].count(sep) for v in variants)
        maxL = max(counts) + 1
        if len(counts) == 1:
            for v in variants:
                v[ind:(ind+1)] = v[ind].split(sep)
        else:
            for v in variants:
                s = v[ind].split(sep)
                if len(s) < maxL: s.extend(['']*(maxL-len(s)))
                v[ind:(ind+1)] = s
        h[ind:(ind+1)] = [h[ind]+'_%d' %i for i in range(1, maxL+1)]
        return h, variants
       
    def _fixAnnovarOtherinfo(self, headers, firstvar):
        print "FIXING OTHERINFO!!"
        lenh = len(headers)
        lenfst = len(firstvar)
        if (lenfst - lenh) == 9 and firstvar[lenh+7].startswith('GT'):
            headers[-1:] = ['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '<sample>'] 
        if (lenfst - lenh) == 10 and firstvar[lenh+8].startswith('GT'):
            headers[-1:] = ['zygosity', 'CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '<sample>'] 
        return headers
        
    def _parseDescriptions(self, lines):
        '''Parses the intro part of a VCF file. Returns a dict with elements like: {(FORMAT, GT): "Description"}.'''
        
        PATTERN = re.compile(r'''((?:[^,"']|"[^"]*"|'[^']*')+)''')
        def process(line):
            col = line[2:line.index('=')]
            start, stop = line.index('<')+1, line.rindex('>')
            dat = dict([b.split('=', 1) for b in PATTERN.split(line[start:stop])[1::2]])
            return ((col, dat['ID']), dat['Description'])

        descr = [process(line) for line in lines if all(x in line for x in ['=<', '>', 'ID', 'Description'])]
        if not descr:
            return None
        descriptions = {key[1] : val for key, val in descr if key[0] == "FORMAT" or key[0] == "INFO" }
        descriptions.update({key[1]+'_INFO' : val for key, val in descr if key[0] == "INFO"})
        filters = [x for x in descr if x[0][0] == "FILTER"]
        if filters:
            max_width = max(len(x[0][1]) for x in filters)
            descriptions['FILTER'] = '\n'.join('%-*s %s' %(max_width, x[0][1] + ':', x[1]) for x in filters)
        
        return descriptions


if __name__ == "__main__":
    import Filter
    import FiltusAnalysis
    autex = FiltusAnalysis.AutExComputer(genmapfile="C:\\Projects\\FILTUS\\DecodeMap_thin.txt")
    dn = FiltusAnalysis.DeNovoComputer()
    
    reader = VariantFileReader()
    #testfile = "C:\\Projects\\FILTUS\\FILTUS_EXE\\example_files\\test_file.csv"
    testfile = "C:\\Projects\\testfiles_all\\Eirik\\Frengen_PE_exome_batch5\\famPE15_PE_batch5_140206.variantsOnly.targetsPad50.ug.vqsrAndHard.SNPEFFtopImpact.vepAllAndTopImpacts.phased.vcf"
    vflist = reader.readVCFlike(testfile, sep="\t", chromCol="CHROM", posCol="POS", geneCol="", splitAsInfo="INFO", keep00=1)
    
    filter = Filter.Filter(filterFile="C:\\Projects\\FILTUS\\PASSfilter.fconfig")
    vflist = map(filter.apply, vflist)
    
    for vf in vflist:
        hom = autex.autex_segments(vf, f=0.01, a=0.5, error=0.005, altFreqCol='VEP_ENSEMBL_ALLELE_FREQ_INFO', defaultFreq=0.5, minlength=1.0, mincount=0)
        print autex.summary(hom, txt=1)
    
    #res = hom.printData(trunc=12)
    #print res[0]
    #print res[1]
    #print dn.analyze(vflist[0], vflist[1], vflist[2], mut=1e-8, freq=('', 0.1)).length
    #print dn.analyze(vflist[1], vflist[2], vflist[0], mut=1e-8, freq=('', 0.1)).length
    #print dn.analyze(vflist[2], vflist[0], vflist[1], mut=1e-8, freq=('', 0.1)).length
    
    #import sys
    #infile = sys.argv[1]
    