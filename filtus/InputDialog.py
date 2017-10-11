import Tkinter
import Pmw
import csv
import FiltusWidgets
import FiltusUtils

import os.path
import VariantFileReader

class InputDialog(object):
    def __init__(self, filtus):
        self.filtus = filtus
        self.parent = filtus.parent
        self.reader = VariantFileReader.VariantFileReader()
        self.currentfile = None
        self.commentChar = '##'
        self.sep = None
        self.currentHeaders = []
        self.originalHeaders = []
        self.firstvariants = []
        
        self.chromCol = None
        self.posCol = None
        self.geneCol = None
        self.gtCol = None
        self.homSymbol = None
        self.infoCol = None
        self.formatCol = None
        
        self.vcf = False
        self._INFOheaders = []
        self._FORMATheaders = []
        self._sampleNames = []
        self.splitFormat = False
        self.hasCSQ = False
        self._CSQheaders = []
        self.splitCsq = False
        self.keep00 = False
        self.split_general = False
        self.splitFormatVar = Tkinter.IntVar(self.parent)
        self.splitCsqVar = Tkinter.IntVar(self.parent)
        self.keep00Var = Tkinter.IntVar(self.parent)
        
        self.prefilter = None
        self.skiplines = 0
        self.prompt = True
        self.guess = False
        self.skipFile = False
        self.stopLoading = False
        self.checkHomozygosity = False

        self._separators = ["comma", "tab", "semicolon", "space"]
        self._sepDic = dict(comma=',', tab = '\t', semicolon = ';', space = ' ')
        self._sepDicInv = {v:k for k, v in self._sepDic.items()}
        self._createDialog()


    def _createDialog(self):
        self.dialog = Pmw.Dialog(self.parent, buttons = ('Use for all files', 'Use for this file', 'Skip this file', 'Cancel'),
                            defaultbutton = 0, title = 'Input file settings', command=self._executeDialogButton,
                            dialogchildsite_pady=0, buttonbox_pady=10)
        self.dialog.withdraw()
        interior0 = self.dialog.interior()
        fr = Tkinter.Frame(interior0) #self.dialog.interior()
        fr.grid(row=0, column=0, pady=10, sticky='news')
        FiltusWidgets.HelpButton(interior0, filtus=self.filtus, page="loading", bookmark="inputsettings").grid(row=0, column=0, sticky="ne")
        
        OM = FiltusWidgets.OptionMenuExt
        filename_group = Pmw.Group(fr, tag_text = 'File name')
        self.fileLabel = Tkinter.Label(filename_group.interior(), justify = "left", font = self.filtus.textfont)

        grid_nw = dict(sticky='nw', padx=10, pady=2)
        grid_nw_right = dict(sticky='nw', padx=(20, 10), pady=2)
        
        self.fileLabel.grid(row=0, column=1, **grid_nw)

        basic_group = Pmw.Group(fr, tag_text = 'Basic settings')
        basic_interior = basic_group.interior()
        pmw_OPTIONS = dict(labelmargin=10, labelpos='w')
        width=15 # if self.filtus.windowingsystem != 'aqua' else 12
        OM_OPTIONS = dict(labelmargin = 10, labelpos='w', menubutton_anchor = 'w', menubutton_padx=5, menubutton_pady=1, 
                          menubutton_width=width, menu_font=self.filtus.defaultfont)

        self.sepInputOM = OM(basic_interior, label_text = "Column separator:", items = self._separators,
                                        command=self._readAndSetHeaders, **OM_OPTIONS)
        self.commentEntry = Pmw.EntryField(basic_interior, label_text = "Preamble lines start with:", value = self.commentChar,
                                modifiedcommand=self._noDefButton, command=self._readAndSetHeaders, entry_width=12, **pmw_OPTIONS)
        self.sepInputOM.grid(row=0, column=0, **grid_nw)
        self.commentEntry.grid(row=0, column=1,  **grid_nw_right)
        
        ### variant settings
        variant_group = Pmw.Group(fr, tag_text = 'Variant settings')
        variant_interior = variant_group.interior()
        self.chromColMenu = OM(variant_interior, label_text = "Chrom column:", **OM_OPTIONS)
        self.posColMenu = OM(variant_interior, label_text = "Position column:", **OM_OPTIONS)
        self.geneColMenu = OM(variant_interior, label_text = "Gene name column:", **OM_OPTIONS)
        
        self.chromColMenu.grid(row=0, column=0, **grid_nw)
        self.posColMenu.grid(row=0, column=1, **grid_nw_right)
        self.geneColMenu.grid(row=1, column=0, **grid_nw)
        
        Tkinter.Frame(variant_interior, height=2, borderwidth=2, relief="sunken").grid(sticky='ew', pady=(8,4), columnspan=2)
        self.vcfChooser = Pmw.RadioSelect(variant_interior, buttontype="radiobutton",
                        labelpos="w", labelmargin=0, label_text = "Genotype format:", command=self._vcfCallback)
        self.vcfChooser.add("VCF")
        self.vcfChooser.add("Other")
        self.vcfChooser.setvalue('Other')
        self.vcfChooser.grid(sticky='w', padx=10, columnspan=2)
        
        self.VCFframe = Tkinter.Frame(variant_interior)
        self.formatColMenu = OM(self.VCFframe, label_text = "FORMAT column:",  **OM_OPTIONS)
        self.formatColMenu.grid(**grid_nw)
        Tkinter.Checkbutton(self.VCFframe, text="Keep 0/0   (only for autozygosity/de novo)", 
            variable=self.keep00Var, anchor='w').grid(row=0, column=1, **grid_nw_right)
            
        self.nonVCFframe = Tkinter.Frame(variant_interior)
        self.gtColMenu = OM(self.nonVCFframe, label_text = "Genotype column:",  **OM_OPTIONS)
        self.homSymbolEntry = Pmw.EntryField(self.nonVCFframe, label_text = "Homozygosity symbol:", entry_width=12, **pmw_OPTIONS)
        self.gtColMenu.grid(**grid_nw)
        self.homSymbolEntry.grid(row=0, column=1, **grid_nw_right)
        
        self.VCFframe.grid(sticky='news', columnspan=2)
        self.VCFframe.grid_remove()
        self.nonVCFframe.grid(sticky='news', columnspan=2)
        
        ### split settings
        split_group = Pmw.Group(fr, tag_text = 'Column splits')#, tag_pyclass = Button, tag_relief = 'raised', tag_command=self.vcfToggle, collapsedsize=10)
        split_interior = split_group.interior()
        self.splitFormatButt = Tkinter.Checkbutton(split_interior, variable=self.splitFormatVar, text="  Split FORMAT/genotypes",
            command=self._splitFORMAT_update)
        self.splitFormatButt.grid(sticky='w', padx=(10,5), pady=2)
        
        self.infoColMenu = OM(split_interior, label_text='Split as "INFO":', command=self._splitINFO_update, **OM_OPTIONS)
        self.infoColMenu.grid(**grid_nw)
        
        self.splitCsqButt = Tkinter.Checkbutton(split_interior, variable=self.splitCsqVar, text="and split CSQ",
            command=self._splitCsq_update)
        self.splitCsqButt.grid(row=1, column=1, **grid_nw)
        
        self.splitcol1Menu = OM(split_interior, label_text="Split column:", **OM_OPTIONS)
        self.splitcol1_sep = Pmw.EntryField(split_interior, label_text = "by separator", entry_width=4, **pmw_OPTIONS)
        self.splitcol2Menu = OM(split_interior, label_text="Split column:", **OM_OPTIONS)
        self.splitcol2_sep = Pmw.EntryField(split_interior, label_text = "by separator", entry_width=4, **pmw_OPTIONS)
        
        self.splitcol1Menu.grid(**grid_nw)
        self.splitcol1_sep.grid(row=2, column=1, **grid_nw)
        self.splitcol2Menu.grid(**grid_nw)
        self.splitcol2_sep.grid(row=3, column=1, **grid_nw)
        
        prefilter_group = Pmw.Group(fr, tag_text = "Prefilter")
        prefilter_interior = prefilter_group.interior()
        prefilter_interior.columnconfigure(1, weight=1)
        self.prefilter_operatorOM = OM(prefilter_interior, label_text = "Keep only lines which ", items = ['', 'contain', 'do not contain', 'start with', 'do not start with'], **OM_OPTIONS)
        self.prefilter_valueEntry = Pmw.EntryField(prefilter_interior, entry_width=10, **pmw_OPTIONS)
        
        self.prefilter_operatorOM.grid(row=0, column=0, **grid_nw)
        self.prefilter_valueEntry.grid(row=0, column=1, sticky='nwe', padx=(0, 10), pady=2)
        
        for g in (filename_group, basic_group, variant_group, split_group, prefilter_group):
            g.configure(ring_borderwidth=1, tag_font = self.filtus.smallbold)
            g.grid(sticky='news', pady=6, padx=10, ipady=2)
        self.align()

    def _vcfCallback(self, button):
        self.vcf = button=='VCF'
        if self.vcf:
            self.splitFormatButt.configure(state="normal")
            self.nonVCFframe.grid_remove()
            self.VCFframe.grid()
        else:
            self.splitFormatButt.deselect()
            self._splitFORMAT_update()
            self.splitFormatButt.configure(state="disabled")
            self.VCFframe.grid_remove()
            self.nonVCFframe.grid()
        
        
    def align(self):
        Pmw.alignlabels([self.sepInputOM, self.chromColMenu, self.vcfChooser, self.geneColMenu, self.formatColMenu, self.gtColMenu, 
            self.prefilter_operatorOM, self.infoColMenu, self.splitcol1Menu, self.splitcol2Menu], sticky='w') 
        Pmw.alignlabels([self.commentEntry, self.homSymbolEntry], sticky='w')
        
        
    def read(self, filename, **kwargs):
        self.skipFile = False
        self.stopLoading = False
        new_ext = self.currentfile is None or (os.path.splitext(filename)[1] != os.path.splitext(self.currentfile)[1])
        self.prompt = kwargs.pop('prompt', self.prompt or new_ext)
        self.guess = kwargs.pop('guess', self.guess or (self.prompt and new_ext))
        
        try:
            self._guessAndPrepare(filename, kwargs)
            if self.prompt or any(OM.inconsistent for OM in self._activeMenus()):
                FiltusUtils.activateInCenter(self.parent, self.dialog)
            else:
                self._setParameters()
            if self.stopLoading or self.skipFile:
                return
            
            self.filtus.busy()
            common_params = dict(filename=filename, sep=self.sep, chromCol=self.chromCol, posCol=self.posCol, geneCol=self.geneCol, 
                                splitAsInfo=self.infoCol, splitCsq=self.splitCsq, split_general=self.split_general, prefilter=self.prefilter)
            
            if self.vcf:
                VF = self.reader.readVCFlike(formatCol=self.formatCol, splitFormat=self.splitFormat, keep00=self.keep00, **common_params)
            else:
                VF = self.reader.readNonVCF(skiplines=self.skiplines, gtCol=self.gtCol, homSymbol=self.homSymbol, **common_params)
            self.filtus.notbusy()
            
        except (ValueError, RuntimeError) as e:
            self.filtus.notbusy()
            FiltusUtils.warningMessage(e)
            return self.read(filename, guess = False, prompt=True)
        except Exception as e:
            raise
            self.filtus.notbusy()
            typ = type(e).__name__
            FiltusUtils.warningMessage("An error occured while reading this file:\n%s\n\n%s: %s\n\nPlease try again or skip file." %(filename, typ, e))
            return self.read(filename, guess=False, prompt=True)
        
        if self.checkHomozygosity and VF.noHomozygotes():
            tryagain = FiltusUtils.yesnoMessage('The file %s has no homozygous variants. Go back to settings dialog?'%filename)
            if tryagain:
                VF = self.read(filename, guess = False, prompt=True)

        return VF
    
    ############# Functions to follow prepare the prompt
        
    def _guessAndPrepare(self, filename, kwargs):
        self.currentfile = filename
        self.fileLabel.configure(text=FiltusUtils.wrapFilename(filename, joinsep='\n     '))
        preamble, headerline, firstline = self._getFirstLines(filename)
        self.__dict__.update(kwargs)
        
        sep = self.sep
        if sep is None or sep not in headerline:
            sep = next((char for char in ['\t', ',', ';', ' '] if char in headerline and char in firstline), '\t')
        
        self.sepInputOM.invoke(self._sepDicInv[sep])
        headers = self.currentHeaders
        
        def _doGuess(col):
            '''Dont guess if specified in arguments, or if the current value is consistent.'''
            return self.guess and col not in kwargs and (getattr(self, col) is None or getattr(self, col+'Menu').inconsistent)
        
        lowheaders = [h.lower() for h  in headers]
        def _matchHeader(alts):
            for h in alts:
                if h in lowheaders: return headers[lowheaders.index(h)]
            return ''
            
        if 'chromCol' in kwargs: self.chromColMenu.setAndCheck(kwargs['chromCol'])
        elif _doGuess('chromCol'):
            chromCol = _matchHeader(['#chrom', 'vcf_chrom', 'vcf_chr', 'chrom', 'chr', 'chromosome'])
            if chromCol: self.chromColMenu.setAndCheck(chromCol)
        
        if 'posCol' in kwargs: self.posColMenu.setAndCheck(kwargs['posCol'])
        elif _doGuess('posCol'):
            posCol = _matchHeader(['pos', 'vcf_pos', 'vcf_start', 'start', 'position', 'pos_start', 'chromosome_position'])
            if posCol: self.posColMenu.setAndCheck(posCol)
        
        if 'splitAsInfo' in kwargs: 
            self.infoColMenu.setAndCheck(kwargs['splitAsInfo'])
            self. _splitINFO_update()
        
        # If VEP CSQ info present: Store headers
        self.hasCSQ = False
        for line in preamble:
            if 'ID=CSQ,' in line:
                self.hasCSQ = True
                self._CSQheaders = line.split("Format: ")[1].strip().strip('">').split("|")
                break
                
        if 'geneCol' in kwargs: self.geneColMenu.setAndCheck(kwargs['geneCol'])
        elif _doGuess('geneCol'):
            geneCol = _matchHeader(['gene', 'gene.refgene', 'gene symbol'])
            if geneCol =='': 
                genecCol = next((h for h, lowh in zip(headers, lowheaders) if 'gene' in lowh and 'name' in lowh), '')
            if geneCol: self.geneColMenu.setAndCheck(geneCol)
        
        if self.guess:
            vcf, infoCol, formatCol = self._guessVCF(self.originalHeaders, self.firstvariants[0])  # infoCol not used
            self.vcfChooser.invoke(int(not vcf))
            if vcf: self.splitFormatVar.set(1) # Default option: Split FORMAT
            
        if 'formatCol' in kwargs: 
            self.formatColMenu.setAndCheck(kwargs['formatCol'])
        elif self.guess: #from above
            self.formatColMenu.setAndCheck(formatCol)
        
        if 'splitFormat' in kwargs: 
            self.splitFormatVar.set(kwargs['splitFormat'])
            
        self._splitFORMAT_update()
            
        if 'keep00' in kwargs: self.keep00Var.set(kwargs['keep00'])
        
        if 'gtCol' in kwargs: self.gtColMenu.setAndCheck(kwargs['gtCol'])
        elif _doGuess('gtCol'):
            gtCol = '' if vcf else _matchHeader(['genotype', 'gt', 'zygosity', 'homozygous', 'attribute_het'])
            self.gtColMenu.setAndCheck(gtCol)
        
        if 'split_general' in kwargs:
            s = kwargs['split_general']
            split, sep = s[0]
            self.splitcol1Menu.setAndCheck(split)
            self.splitcol1_sep.setvalue(sep)
            if len(s) > 1:
                split, sep = s[1]
                self.splitcol2Menu.setAndCheck(split)
                self.splitcol2_sep.setvalue(sep)
        
        if 'prefilter' in kwargs:
            operatorText, value = kwargs['prefilter']
            self.prefilter_operatorOM.setAndCheck(operatorText)
            self.prefilter_valueEntry.setvalue(value)
            
    def _getFirstLines(self, filename, n=1):
        self.skiplines = 0
        preamble = []
        firstlines = []
        self.commentChar = comment = self.commentEntry.getvalue().strip()
        with open(filename, 'rU') as ifile:
            for line in ifile:
                if comment and line.startswith(comment): 
                    preamble.append(line)
                    self.skiplines += 1
                    continue
                firstlines.append(line)
                if len(firstlines) > n: 
                    break
        if not firstlines or not firstlines[0].strip():
            raise IOError("Skipping empty file: %s" %filename)
        headerline = firstlines[0]
        if n==1:
            first = firstlines[1] if len(firstlines) > 1 else ''
        else:
            first = firstlines[1:]
        return preamble, headerline, first
 
 
    ################### Callback functions
    
    def _readAndSetHeaders(self, sepvalue=None):
        '''Callback for both self.sepInputOM and self.commentEntry'''
        preamble, headerline, firstline = self._getFirstLines(self.currentfile, n=100)
        if sepvalue: 
            self.sep = self._sepDic[sepvalue]
        self.sepInputOM.setColor(test=self.sep in headerline)
        
        top = csv.reader([headerline] + firstline, delimiter=self.sep, skipinitialspace=True)
        h = top.next()
        self.firstvariants = list(top)
        
        if h[-1] == "Otherinfo": ### Annovar fix: Re-inserting VCF column names.
            h[:] = self.reader._fixAnnovarOtherinfo(h, self.firstvariants[0])
        if h[0] == '#CHROM': ### VCF tweak
            h[0] = 'CHROM'
        
        self.originalHeaders = h
        
        self._updateColnameMenus(h, all=True)
        self._splitINFO_update(reset=True)
        self._splitFORMAT_update(reset=True)
        
        
    def _splitINFO_update(self, column=None, reset=False): 
        '''callback for the INFO option menu. Also called from _readAndSetHeaders (with column=None)'''
        if reset:
            self._INFOheaders = []
        if column is None:
            if self.infoColMenu.inconsistent:
               return
            column = self.infoColMenu.getvalue()
        
        self.infoColMenu.setColor(True)
        h = self.currentHeaders[:]
        
        ### Always start by unsplitting everything:
        # If CSQ is split: unsplit this first
        splitCsq = self.hasCSQ and self.splitCsqVar.get()
        if splitCsq:
            self.splitCsqVar.set(0)
            self._splitCsq_update()
            h[:] = self.currentHeaders[:]

        # Unsplit INFO fields
        if self._INFOheaders:
            ind = h.index(self._INFOheaders[0])
            h[ind:(ind + len(self._INFOheaders))] = [self.infoCol]
        self._INFOheaders = []
        self.infoCol = ''
           
        ### If empty selection: Reset and return
        if column == "": 
            self._updateColnameMenus(h)  
            self.splitCsqButt.configure(state="disabled")
            return
        
        ### Otherwise: split selected column as INFO (if possible)
        first_infos = [v[self.originalHeaders.index(column)] for v in self.firstvariants]
        _INFOheaders = sorted(set(s.split('=')[0] + '_INFO' for info in first_infos for s in info.split(';') if '=' in s))
        
        if not _INFOheaders:
            self.infoColMenu.setColor(False)
            self.splitCsqButt.configure(state="disabled")
            self._updateColnameMenus(h)
            FiltusUtils.warningMessage("I don't recognise %s as an INFO column"%column)
            return
        
        ind = h.index(column)
        h[ind:(ind + 1)] = _INFOheaders
        self._updateColnameMenus(h)   
        
        if self.hasCSQ and "CSQ_INFO" in _INFOheaders:
            self.splitCsqButt.configure(state="normal")
            if splitCsq:
                self.splitCsqVar.set(1)
                self._splitCsq_update()
        
        self._INFOheaders = _INFOheaders
        self.infoCol = column
        
        
    def _splitCsq_update(self, reset=False): 
        '''callback for the splitCsq checkbox'''
        h = self.currentHeaders[:]
        csq_heads = self._CSQheaders
        split = self.splitCsqVar.get()
        if split:
            ind = h.index("CSQ_INFO")
            h[ind:(ind + 1)] = csq_heads
        else:
            if all(csq in h for csq in csq_heads): # extra precaution
                ind = h.index(csq_heads[0])
                h[ind:(ind + len(csq_heads))] = ['CSQ_INFO']
        
        self._updateColnameMenus(h)    
        
    def _splitFORMAT_update(self, reset=False): 
        '''callback for the splitFormat checkbox'''
        if reset:
            self._FORMATheaders = []
            self._sampleNames = []
        split = self.splitFormatVar.get()
        column = self.formatColMenu.getvalue()
        if split and (not column or self.formatColMenu.inconsistent):
            self.splitFormatVar.set(0)
            return
        
        h = self.currentHeaders[:]
        
        def unsplit():
            if self._FORMATheaders:
                h[h.index('GT'):] = [self.formatCol] + self._sampleNames
            self._FORMATheaders = []
            self._sampleNames = []
        
        if split: 
            first = self.firstvariants[0][self.originalHeaders.index(column)]
            if not first.startswith('GT'):
                self.formatColMenu.setColor(False)
                FiltusUtils.warningMessage("FORMAT column entries must begin with 'GT'")
                return
            unsplit() # undo possible previous split
            self.formatCol = column
            self._FORMATheaders = first.split(':')
            ind = h.index(column)
            self._sampleNames = h[ind+1:]
            h[ind:] = self._FORMATheaders
        else:
            unsplit()
            
        self._updateColnameMenus(h)
    
    def _updateColnameMenus(self, headers, all=False):
        '''Update various option menus in the dialog. If (all) include those unaffected by splits.'''
        self.currentHeaders = headers[:]
        for OM in self._activeMenus(fixed=all):
            OM.setItems([''] + headers)
        
    def _executeDialogButton(self, button):
        try:
            if button is None or button == 'Cancel':
                self.stopLoading = True
                self.dialog.deactivate()
                return
            elif button == "Skip this file":
                self.skipFile = True
                self.dialog.deactivate()
                return
            self.prompt = button != "Use for all files"  #button is either this or "Use for this file"
            self.guess = False
            
            try:
                self._setParameters()
            except Exception as e:
                FiltusUtils.warningMessage(e)
                return
        
            self.dialog.deactivate()
        except Exception as e:
            FiltusUtils.warningMessage("Something went wrong. Trying to close the input dialog.")
            self.dialog.destroy()
            del self.filtus.fileReader
            return
            
    def _setParameters(self):
        if self.sepInputOM.inconsistent:
            raise RuntimeError('Column separator not found in first line.\n\nPlease check the input settings (including "Skip lines starting with")')
        wrongcols = [OM.getvalue() for OM in self._activeMenus() if OM.inconsistent and OM not in (self.splitcol1Menu, self.splitcol2Menu)]
        if wrongcols: 
            raise RuntimeError("Column(s) not found in file: %s" %', '.join(wrongcols))
        for key in ['chromCol', 'posCol', 'geneCol', 'formatCol', 'gtCol', 'infoCol', 'splitcol1', 'splitcol2']:
            setattr(self, key, getattr(self, key+'Menu').getvalue())
        
        self.splitFormat = self.splitFormatVar.get()
        self.splitCsq = self.splitCsqVar.get()
        self.keep00 = self.keep00Var.get()
        self.split_general = [(self.splitcol1, self.splitcol1_sep.getvalue()), (self.splitcol2, self.splitcol2_sep.getvalue())]
        self.split_general = [(x,y) for x,y in self.split_general if x]
        
        self.prefilter = (self.prefilter_operatorOM.getvalue(), self.prefilter_valueEntry.getvalue())
        if self.prefilter[0] == '' or self.prefilter[1] == '': 
            self.prefilter = None
            
        if not self.chromCol: raise RuntimeError("Please indicate chromosome column.")
        if not self.posCol: raise RuntimeError("Please indicate position column.")
        
        if self.vcf:
            if not self.formatCol: raise RuntimeError("Please indicate vcf-like FORMAT column.")
        else:
            self.homSymbol = self.homSymbolEntry.getvalue()
            if self.gtCol and not self.homSymbol: 
                raise RuntimeError("Missing symbol for homozygous genotype.")
            if self.homSymbol and not self.gtCol: 
                raise RuntimeError("Symbol of homozygosity given, but no genotype column.")
        
        for x in [x for x,y in self.split_general if not y]:
            raise RuntimeError("Please indicate splitting separator for column '%s'." % x)
        if len(self.split_general)==2 and self.split_general[0][0] == self.split_general[1][0]:
            raise RuntimeError("Column cannot be split twice: '%s'." % self.split_general[0][0])
    
        
    def _noDefButton(self):
        '''this is invoked when modifying the commentCharEntry, to stop <Return> from jump to the default button.'''
        self.dialog.component('buttonbox').setdefault(None)

    
    def _guessVCF(self, headers, firstvar):
        '''Returns (vcf [T/F], infoCol, formatCol) Requires "format" in column name AND matching numbers of colons in all remaining cols'''
        n = len(headers)
        lowheads = [h.lower() for h in headers]
        kolon = [x.count(':') for x in firstvar]
        for i in range(n-1, 1, -1): # shorter to go backwards 
            if not 'format' in lowheads[i]: continue
            k = kolon[i]
            if k==0 and not firstvar[i]=='GT': continue
            if all(0 <= kolon[j] <= k or firstvar[j]=='./.' for j in range(i+1, n)):
                infoCol = headers[i-1] if 'info' in lowheads[i-1] else ''
                formatCol = headers[i]
                return (True, infoCol, formatCol)
        return (False,'','')

     
    def _activeMenus(self, fixed=True):
        m = [self.chromColMenu, self.posColMenu, self.geneColMenu, self.gtColMenu, self.splitcol1Menu, self.splitcol2Menu]
        if fixed: m.extend([self.infoColMenu, self.formatColMenu])
        return m
     
        