import re
import bisect  
import tkMessageBox
import datetime

def warningMessage(message):
    tkMessageBox.Message(icon='warning', type = 'ok', message = message, title = "Warning").show()

def infoMessage(message):
    tkMessageBox.showinfo(title='Info', message=message)

def yesnoMessage(message, **kw):
    answer = tkMessageBox.askyesno(title="Warning", message=message, **kw)
    return answer

def activateInCenter(parent, widget):
    width, height, offx, offy = parent.winfo_width(), parent.winfo_height(), parent.winfo_rootx(), parent.winfo_rooty()
    x = offx + width/2 - widget.winfo_reqwidth()/2
    y = offy + height/8
    return widget.activate(geometry='+%d+%d' %(x, y))

def composeMeta(VFlist=None, VFindex=None, sort=True, analysis=None, appendTo=None, version='1.0.2'):
    txt = '###############################################################\n' \
          '####  FILTUS %s, analysis performed %s  ####\n' \
          '###############################################################\n' \
          '##\n' %(version, datetime.datetime.now().strftime("%c"))
    
    if VFlist:
        if sort and VFindex is not None:
            VFindex, VFlist = zip(*sorted(zip(VFindex, VFlist)))
        elif VFindex is None:
            VFindex = range(len(VFlist))
    
        samplenames = ['%d: %s' %(i + 1, VF.longName) for i, VF in zip(VFindex, VFlist)]
        txt += '\n## '.join(["## SAMPLES:"] + samplenames + ['\n'])

        ### Prefilters
        all_prefilters = [VF.prefilter for VF in VFlist]
        prefilters = list({f for f in all_prefilters if f})
        
        if len(prefilters) == 0: 
            txt += "## PREFILTERS: None\n##\n"
        elif all(f == prefilters[0] for f in all_prefilters):
            filtertext = 'Keep lines which %s %r' %prefilters[0]
            if len(VFlist) == 1: 
                txt += "## PREFILTER APPLIED:\n## %s\n##\n" % filtertext
            else:
                txt += "## PREFILTER APPLED TO ALL SAMPLES:\n## %s\n##\n" % filtertext
        else:
            txt += "## PREFILTERS:\n"
            # list of all pairs (ids, fi) where ids are the samples having the filter fi
            IDs = [([i + 1 for i, f in enumerate(all_prefilters) if f==pref], pref) for pref in prefilters]
            IDs.sort(key=lambda x: x[0][0])
            
            for ids, pref in IDs:
                filtertext = 'Keep lines which %s %r' %pref
                txt += "## Applied to %s: %s\n" %(', '.join(map(str,ids)), filtertext)
            txt += "##\n"
                
        ### Filters
        all_filters = [VF.appliedFilters for VF in VFlist]
        filters = list({f for f in all_filters if f and f.details()})
        
        if len(filters) == 0:
            txt += "## FILTERS: None\n##\n"
        elif len(filters) == 1:
            f0 = all_filters[0]
            if len(VFlist) == 1:
                txt += "## FILTERS APPLIED:\n%s\n##\n" %str(f0)
            elif all(f is f0 for f in all_filters):
                txt += "## FILTERS APPLIED TO ALL SAMPLES:\n%s\n##\n" %str(f0)
        else:
            for filt in filters:
                IDs = [str(i + 1) for i, VF in zip(VFindex, VFlist) if VF.appliedFilters is filt]
                if len(IDs)==1:
                    txt += "## FILTERS APPLIED TO SAMPLE %s:\n%s\n##\n" %(IDs[0], str(filt))
                else:
                    txt += "## FILTERS APPLIED TO SAMPLES %s:\n%s\n##\n" %(', '.join(IDs), str(filt))
    
    if analysis:
        txt += "## ANALYSIS: %s\n##\n" %analysis
    if appendTo:
        txt = appendTo + '\n' + txt
    return txt
        
        
        
def interval_union(x):
    '''
    Union of intervals.
    x: a sequence of intervals of the form (start, end) or [start, end]
    output: list of union intervals
    '''
    z = sorted(x, reverse=1) # to use pop()
    res = []
    a,b = z.pop()
    while z:
        c,d = z.pop()
        if c <= b:
            b = max(b,d)
        else: 
            res.append([a,b])
            a,b = c,d
    res.append([a,b])
    return res

def interval_set_intersection(x, y):
    '''
    Intersection of two interval sets.
    x, y: sequences of intervals of the form (start, end) or [start, end]
    output: list of intersection intervals
    '''
    x = interval_union(x) # the algorithm below assumes that the intervals in x (and y) are disjoint
    y = interval_union(y)
    z = sorted(x+y, reverse=1) # to use pop()
    res = []
    a,b = z.pop()
    while z:
        c,d = z.pop()
        if c <= b: # non-empty intersection (assumes each of x and y are disjoint!)
            res.append([c, min(b,d)])
            a,b = min(b,d), max(b,d)
        else:
            a,b = c,d
    return res

def region_intersection(*x):
    commonChr = sorted(set.intersection(*[set(r[0] for r in regs) for regs in x]))
    res = []
    for chr in commonChr:
        chrRes = [r[1:] for r in x[0] if r[0]==chr]
        for xx in x[1:]:
            y = [r[1:] for r in xx if r[0]==chr]
            chrRes = interval_set_intersection(chrRes, y)
            if not chrRes: break
        res.extend([[chr]+r for r in chrRes])
    return res
    
def inGeneList(entry, genelist):
    if entry in genelist:
        return True
    elif any(sep in entry for sep in [',', ';', '(']):
        return any(entr in genelist for entr in re.split('[,;(]', entry))
    else: return False

def ignore_pass(event):
    pass

def ignore_break(event):
    return "break"

def convertType(s):
    for func in (int, float):
        try:
            n = func(s)
            return n
        except: pass
    return s

def tryFloat(s, stringval=float('inf')):
    try:
        num = float(s)
    except: 
        num = stringval
    return (num, s)
    
def checkDuplicates(fields):
    '''given list of numeric vectors: raise error if any number appears twice'''
    flatlist = set(i for f in fields for i in f)
    for i in flatlist:
        if sum(i in f for f in fields) > 1:
            raise ValueError("Same sample number indicated in several fields: %d" %(i+1,))
      
def chromInt(s):
    try:
        return int(s)
    except ValueError:
        if s == 'X': return 23
        if s == 'Y': return 24
        if s.startswith('chr'): return chromInt(s[3:])
        else: return 25

def XminusPAR(vdef):
        '''input: (chr, pos), both strings'''
        if chromInt(vdef[0]) == 23:
            pos = int(vdef[1])
            return not (60000<pos<2699521 or 154931043<pos<155270561) #X:0-60000; X:2699521-154931043; X:155270561-1e9
        return False
        
def string2intlist(txt, idlist=None):
    #example: '1, 4-6, 8-10' --> [1, 4, 5, 6, 8, 9, 10]
    txt = txt.strip()
    if not txt: return []

    res = []
    a = [i.split('-') for i in txt.split(',')]
    for s in a:
        if len(s) == 1: res.append(int(s[0]))
        elif len(s) == 2:
            start = int(s[0])
            stop = int(s[1])
            if not start < stop:
                raise ValueError("Invalid range: %s"%'-'.join(s))
            res.extend(range(start, stop + 1))
    return res

def convert2indices(txt, idlist):
    exclude = False
    txt = txt.strip()
    if txt.startswith("NOT"):
            exclude = True
            txt = txt[3:].strip()
    res = []

    if txt.startswith("ID"):
        nomatch = []
        uniq = [id.strip() for id in txt[2:].split(',') if id.strip()]

        for id in uniq:
            tf = [id in name for name in idlist]
            if sum(tf) == 0: nomatch.append(id)
            elif sum(tf) > 1: raise ValueError("Identifier '%s' matches several sample names."%id)
            else:
                res.append((i for i, bool in enumerate(tf) if bool).next())

        if len(nomatch)>0:
            raise ValueError("Identifiers not matching any sample name: %s" %', '.join(nomatch))

        nonexist = [str(j) for j in res if j not in range(len(idlist))]
    else:
        res = [i-1 for i in string2intlist(txt)]
        nonexist = [str(j + 1) for j in res if j not in range(len(idlist))]

    if len(nonexist) > 0:
        raise ValueError("Nonexistent sample number(s): %s" %', '.join(nonexist))
    if len(res) > len(set(res)):
        raise ValueError("Repeated number(s) in entered text.")

    if exclude:
        res = [i for i in range(len(idlist)) if i not in res]

    return res


def intlist2string(x, space=False):
    if len(x) == 0: return ''
    xs = sorted(map(int, x))
    N = len(xs)
    i, intstart, res = 0, xs[0], str(xs[0])
    while True:
        k = 1
        while i + k < N and xs[i + k] == intstart + k:  k += 1
        if k == 2: res += ', %d' %(intstart + 1,) if space else ',%d' %(intstart + 1,)
        elif k > 2: res += '-%d' %(intstart + k-1,)
        i += k
        if i < N:
                intstart = xs[i]
                res += ', %d' %intstart if space else ',%d' %intstart
        else: break
    return res

#### Binomial test NB!! Naive implementation, only for small sample size
def choose(n, k):
    if not 0 <= k <= n: return 0
    A = 1
    for t in xrange(min(k, n - k)):
        A = (A * (n - t)) // (t + 1)
    return A

def binomTest(n, k, p): #Prob(X>= k), where X~Bin(n, p)
    return sum(choose(n, i)*p**i*(1-p)**(n-i) for i in range(k, n + 1))

def cummax(v):
    a = [v[0]]
    for i in range(1, len(v)):
        a.append(max(a[-1], v[i]))
    return a

def cummin(v):
    a = [v[0]]
    for i in range(1, len(v)):
        a.append(min(a[-1], v[i]))
    return a

def pValue(m, Lrel, n, k, model, approx=False): #m = average vars after filtr; Lrel = relative length of gene in question; n = #samples; k = #successes
    if model == 'Dominant':
        bernouP = max(m * Lrel, 1.0) if approx else 1-(1.0-Lrel)**m
    elif model.startswith('Recessive'):
        bernouP = m*(m-1.0)*Lrel * Lrel if approx else 1-(1.0-Lrel)**m -m * Lrel*((1.0-Lrel)**(m-1))
    else:
        return -1
    pval = binomTest(n, k, bernouP)
    if pval < -1: pass #print n, k, bernouP
    return binomTest(n, k, bernouP)

def holm(pvs, n=None):
    if n is None: n = len(pvs)
    pvs_decor = [(p, i) for i, p in enumerate(pvs)]
    pvs_decor.sort()
    adj = cummax([min(1, pi[0]*(n-j)) for j, pi in enumerate(pvs_decor)])
    return [adj[pi[1]] for j, pi in enumerate(pvs_decor)]


def hoch(pvs, n=None):
    if n is None: n = len(pvs)
    pvs_decor = [(p, i) for i, p in enumerate(pvs)]
    pvs_decor.sort()
    adj1 = [min(1, pi[0]*(n-j)) for j, pi in enumerate(pvs_decor)]
    adj = list(reversed(cummin(list(reversed(adj1)))))
    return [adj[pi[1]] for j, pi in enumerate(pvs_decor)]

#### Operators:
def myeq(a, b):
    if not a: return False
    if ' OR ' in b:
        return any(bb.strip() == a for bb in b.split(' OR '))
    return b == a

def myne(a, b):
    if not a: return False
    if ' OR ' in b:
        return all(bb.strip() != a for bb in b.split(' OR '))
    return b != a

def contains(a, b):
    if not a: return False
    if ' AND ' in b:
        return all(bb.strip() in a for bb in b.split(' AND '))
    if ' OR ' in b:
        return any(bb.strip() in a for bb in b.split(' OR '))
    if b.startswith('REGEX'):
        return bool(re.search(b[6:], a))
    return b in a

def not_contains(a, b):
    if not a: return False
    if ' AND ' in b:
        return not all(bb.strip() in a for bb in b.split(' AND '))
    if ' OR ' in b:
        return not any(bb.strip() in a for bb in b.split(' OR '))
    if b.startswith('REGEX'):
        return not bool(re.search(b[6:], a))
    return not b in a

def mystartswith(a, b):
    if not a: return False
    if ' OR ' in b:
        return any(a.startswith(bb.strip()) for bb in b.split(' OR '))
    return a.startswith(b)

def not_mystartswith(a, b):
    if not a: return False
    if ' OR ' in b:
        return not any(a.startswith(bb.strip()) for bb in b.split(' OR '))
    return not a.startswith(b)

def myendswith(a, b):
    if not a: return False
    if ' OR ' in b:
        return any(a.endswith(bb.strip()) for bb in b.split(' OR '))
    return a.endswith(b)

def not_myendswith(a, b):
    if not a: return False
    if ' OR ' in b:
        return not any(a.endswith(bb.strip()) for bb in b.split(' OR '))
    return not a.endswith(b)

def floatgt(a, b):
    try:
        return float(a) > b
    except:
        if not a: return False
        elif a == 'X': return 23 > b
        elif a == 'Y': return 24 > b
        else: return False

def floatlt(a, b):
    try:
        return float(a) < b
    except:
        if not a: return False
        elif a == 'X': return 23 < b
        elif a == 'Y': return 24 < b
        else: return False

def missor_eq(a, b):
    if not a or a=="NA": return True
    if ' OR ' in b:
        return any(bb.strip() == a for bb in b.split(' OR '))
    return a == b

def missor_ne(a, b):
    if not a or a=="NA": return True
    if ' OR ' in b:
        return all(bb.strip() != a for bb in b.split(' OR '))
    return a!= b

def missor_contains(a, b):
    if not a or a=="NA": return True
    if ' AND ' in b:
        return all(bb.strip() in a for bb in b.split(' AND '))
    if ' OR ' in b:
        return any(bb.strip() in a for bb in b.split(' OR '))
    if b.startswith('REGEX'):
        return bool(re.search(b[6:], a))
    return b in a

def missor_not_contains(a, b):
    if not a or a=="NA": return True
    if ' AND ' in b:
        return not all(bb.strip() in a for bb in b.split(' AND '))
    if ' OR ' in b:
        return not any(bb.strip() in a for bb in b.split(' OR '))
    if b.startswith('REGEX'):
        return not bool(re.search(b[6:], a))
    return not b in a

def missor_mystartswith(a, b):
    if not a or a=="NA": return True
    if ' OR ' in b:
        return any(a.startswith(bb.strip()) for bb in b.split(' OR '))
    return a.startswith(b)

def missor_not_mystartswith(a, b):
    if not a or a=="NA": return True
    if ' OR ' in b:
        return not any(a.startswith(bb.strip()) for bb in b.split(' OR '))
    return not a.startswith(b)

def missor_myendswith(a, b):
    if not a or a=="NA": return True
    if ' OR ' in b:
        return any(a.endswith(bb.strip()) for bb in b.split(' OR '))
    return a.endswith(b)

def missor_not_myendswith(a, b):
    if not a or a=="NA": return True
    if ' OR ' in b:
        return not any(a.endswith(bb.strip()) for bb in b.split(' OR '))
    return not a.endswith(b)

def missor_floatgt(a, b):
    if a=="" or a=="NA": return True
    try:
        return float(a) > b
    except:
        #if not a: return True
        if a == 'X': return 23 > b
        elif a == 'Y': return 24 > b
        else: return True

def missor_floatlt(a, b):
    if a=="" or a=="NA": return True
    try:
        return float(a) < b
    except:
        #if not a: return True
        if a == 'X': return 23 < b
        elif a == 'Y': return 24 < b
        else: return True

def listUnique(x):
    seen = set()
    return [y for y in x if y not in seen and not seen.add(y)]

def makeUniqueNames(x, remove_underscores=False):
    uniq = []
    for e in x:
        ename = e
        if '_' in ename and remove_underscores:
            ename = ''.join(ename.split('_'))
        if ename in uniq:
            i = 2
            while ename + str(i) in uniq:
                i += 1
            ename += str(i)
        uniq.append(ename)
    return uniq
    
def wrapFilename(fname, lim=70, joinsep='\n'):
    '''Wraps long file names (tries to break at either \, /, ., _, or - as close to "lim" as possible)'''
    
    def wrap(txt):
        if len(txt) <= lim: return txt, ''
        s = txt[lim:30:-1].find('\\') #find last instance of \ between 30 and 70. Returns -1 if none.
        if s<0:
            s = txt[lim:30:-1].find('/')
        if s<0 and '.' in txt[lim:]:    #don't break at last '.' (which is probably the file extension)
            s = txt[lim:30:-1].find('.')
        if s<0:
            s = txt[lim:30:-1].find('_')
        if s<0:
            s = txt[lim:30:-1].find('-')
        if s<0:
            if len(txt) > lim + 10:
                s = 0      #break at 'lim' if no better options
            else:         #if txt is just a bit too long - let it go.
                return txt, ''
        splt = lim - s + 1
        return txt[:splt], txt[splt:]

    parts = []
    remainder = fname
    while(remainder):
        beforesplit, remainder = wrap(remainder)
        parts.append(beforesplit)
    
    return joinsep.join(parts)

def piecewise_linear(map, aval):
    ''' 
    map: sortert matrise med verdier for to variabler (a og b)
    aval: verdi for 1. variabel
    output: verdi for b tilsvarende aval ved stykkevis linear avbildning.. Maksimalt lik siste b-verdi. 
    Kan nok vektoriseres og bli raskere med numpy.
    '''
    a,b = map
    i = bisect.bisect(a, aval) - 1
    ai, bi = a[i], b[i]
    try:
        return bi + (float(aval) - ai) / (a[i+1] - ai) * (b[i+1] - bi)
    except IndexError:
        return b[-1]