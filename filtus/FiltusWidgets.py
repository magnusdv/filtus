# -*- coding: latin-1 -*-
import os
import os.path
import re
import gc
from operator import itemgetter

import Tkinter
import tkFileDialog
import webbrowser
import Pmw

import Filter
import FiltusUtils
import FiltusAnalysis
import InputDialog
try:
    import FiltusQC
    #import FiltusRelatedness
    PLOT_available = 1
except ImportError:
    PLOT_available = 0
import pprint

class HelpButton(Tkinter.Button):
    def __init__(self, parent, filtus, page, bookmark=None, text='?', font=None):
        self.filtus = filtus
        if font is None: font=filtus.defaultfont
        Tkinter.Button.__init__(self, parent, text=text, font=font, padx=2, pady=0, bd=1, height=1, highlightthickness=0, command=self.showInBrowser)
        self.pagepath = os.path.join(filtus.manualdir, page + '.html')
        if bookmark is not None:
            self.pagepath += '#%s'%bookmark
            
    def showInBrowser(self):
        self.filtus.openWebPage(self.pagepath)


class PedWriter(Pmw.Dialog):
    def __init__(self, filtus):
        Pmw.Dialog.__init__(self, filtus.parent, buttons = ('OK', 'Cancel'), title = 'File converter', 
                            activatecommand=self.prepare, command=self.execute)
        self.withdraw()
        self.filtus = filtus
        self.VFlist = []
        interior = self.interior()
        interior.columnconfigure(0, weight=1)
        interior.rowconfigure(1, weight=1)
        
        # Title label
        Tkinter.Label(interior, text='Create ped/dat/map/freq files', font=filtus.titlefont).grid(padx=20, pady=10)
        
        entry_opt = dict(entry_width=8, labelpos='w', labelmargin=5)
        grid_opt = dict(padx=10, pady=5, sticky='news')
        
        ht = max(8, len(filtus.files))
        ped_group = Pmw.Group(interior, tag_text = 'Pedigree data')
        ped_interior = ped_group.interior()
        ped_interior.rowconfigure(0, weight=1)
        
        self.pedtext = Pmw.ScrolledText(ped_interior, rowheader = 1, columnheader = 1, rowcolumnheader = 1, labelpos = 'nw', label_text='Write/paste pedigree columns:',
                text_width = 50, columnheader_width=50, text_height = ht, rowheader_height = ht, 
                rowheader_width = 6, rowcolumnheader_width = 6,
                text_wrap='none', Header_wrap='none', text_font = filtus.monofont, Header_font = filtus.monofont, 
                text_padx = 2,  text_pady = 2, Header_padx = 2, Header_pady = 2)
        
        self.pedtext.component('rowcolumnheader').grid(padx=(0,8))
        self.pedtext.component('rowheader').grid(padx=(0,8))
        
        headerLine = '\t'.join(['FAMID','IID','FID','MID','SEX','AFF'])
        self.pedtext.component('columnheader').insert('0.0', headerLine)
        self.pedtext.component('rowcolumnheader').insert('end', 'Sample')
        self.pedtext.configure(columnheader_state = 'disabled', rowcolumnheader_state = 'disabled')
        self.sampleIndex = self.pedtext.component('rowheader')
        
        self.sampleIndex.configure(foreground="red")
        self.pedtext.grid(row=0, column=0, **dict(grid_opt, pady=(5,10)))
        #Tkinter.Button(ped_interior, text="Reset", command=self._setDefaultPed, pady=0, highlightthickness = 0, font = filtus.smallfont).grid(sticky='ne', row=0, column=0, padx=10)
        
        # Frequency input (NB: Must be defined before whatButt in save group)
        freq_group = Pmw.Group(interior, tag_text = 'Allele frequencies')
        freq_interior = freq_group.interior()
        
        self._altFreqMenu = OptionMenuExt(freq_interior, label_text = "ALT frequency column:", labelpos='nw', 
                                        menu_font=filtus.defaultfont, menubutton_anchor = 'w', labelmargin=3,
                                        menubutton_padx=5, menubutton_pady=1, menubutton_width=15)
        self._def_freq_entry = Pmw.EntryField(freq_interior, label_text = "Missing entry value:", entry_justify="center",
                                        value='', entry_width=6, labelpos='nw', labelmargin=5)
        
        self._altFreqMenu.grid(row=0, column=0, sticky='news', padx=(10,25), pady=(5,10))
        self._def_freq_entry.grid(row=0, column=1, sticky='nes', padx=10, pady=(5,11))
        
        # Map input
        map_group = Pmw.Group(interior, tag_text = 'Map options')
        map_interior = map_group.interior()
        
        self.mapCheck = Pmw.RadioSelect(map_interior, buttontype = 'checkbutton', orient = 'vertical', command=self._whatToggle)
        self.mapCheck.add('phys', text="Use physical positions")
        self.mapCheck.add('cm', text="Convert to centiMorgan (using Decode recombination map)")
        self.mapCheck.invoke('cm')
        self.mapCheck.grid(**grid_opt)
        
        ## Output
        save_group = Pmw.Group(interior, tag_text = 'Output')
        save_interior = save_group.interior()
        
        self.whatButt = Pmw.RadioSelect(save_interior, buttontype = 'checkbutton', orient = 'horizontal', command=self._whatToggle)
        for txt in ['ped', 'dat', 'map', 'freq']: 
            self.whatButt.add(txt)
            self.whatButt.invoke(txt)
        self.whatButt.button(0).configure(state="disabled")
        
        self.dirBrowser = FileBrowser(save_interior, filtus=filtus, label="Save files to directory:", 
                checkbutton = False, labelpos='nw', browsesticky='se', entryfield_entry_width=20, browsetitle="")
        self.dirBrowser.browsebutton.configure(command = self._browseDir)
        self.dirBrowser.entryfield.configure(command = None)
        
        self.prefixEntry = Pmw.EntryField(save_interior, entry_width=15, value="filtus2merlin", 
                                label_text = "Prefix:", labelpos='nw', labelmargin=5)
        
        self.dirBrowser.grid(row=1, column=0, sticky='news', padx=(10,25), pady=(5,10))
        self.prefixEntry.grid(row=1, column=1, sticky='nes', padx=10, pady=(5,10))
        self.whatButt.grid(row=2, columnspan=2, **grid_opt)
        
        
        for g in [ped_group, save_group, map_group, freq_group]:
            g.interior().columnconfigure(0, weight=1)
            #g.configure(ring_borderwidth=1)
            g.configure(tag_font = filtus.smallbold)
            g.grid(**grid_opt)
        
    def _whatToggle(self, button, selected):
        if button=="map":
            map_buttons = [self.mapCheck.button(i) for i in range(2)]
            if selected: 
                for b in map_buttons: b.configure(state="normal")
            else:
                for b in map_buttons: 
                    b.deselect()
                    b.configure(state="disabled")
        if button=="freq":
            if selected: 
                self._altFreqMenu.configure(menubutton_state="normal")
                self._def_freq_entry.configure(entry_state="normal")
            else:
                self._altFreqMenu.configure(menubutton_state="disabled")
                self._def_freq_entry.configure(entry_state="disabled")
            
    def _browseDir(self):
        dir = tkFileDialog.askdirectory(initialdir=self.filtus.currentDir, title = "Save database as")
        if dir:
            self.filtus.currentDir = dir
            self.dirBrowser.setvalue(dir)
            
    def checkPed(self, pedtxt):
        lines = [l.strip() for l in pedtxt.split('\n')]
        while not lines[-1]: lines = lines[:-1]
        pedN = len(lines)
        set012 = set([0,1,2])
        rows = [map(int, line.split()) for line in lines]
        if not all(len(r)==6 for r in rows):
            raise ValueError("All rows must have 6 numerical entries.")
        cols = zip(*rows)
        if not all(x in set012 for x in cols[4]):
            raise ValueError("SEX must be 1 (male), 2 (female) or 0( unknown).")
        if not all(x in set012 for x in cols[5]):
            raise ValueError("AFF must be  1 (non-affected), 2 (affected) or 0 (unknown).")
        return pedN
    
    def _setDefaultPed(self):
        n = len(self.filtus.files)
        ped = '\n'.join(['%d\t%d\t0\t0\t1\t1'%(i+1, i+1) for i in range(n)])
        self.pedtext.settext(ped)
        self.sampleIndex.delete('0.0', 'end')
        self.sampleIndex.insert('end', '\n'.join(str(i+1) for i in range(n)))
        
    def prepare(self):
        if not self.pedtext.getvalue().strip():
            self._setDefaultPed()
        cols = [head for VF in self.filtus.files for head in VF.columnNames]
        self._altFreqMenu.setItems([''] + FiltusUtils.listUnique(cols))
    
    def sampleIndexFix(self, sampletxt, pedN):
        splt = [s.strip() for s in sampletxt.split('\n')]
        if len(splt) > pedN:
            if any(splt[pedN:]):
                raise IndexError("Sample column is too long.")
            splt = splt[:pedN]
        elif len(splt) < pedN:
            splt += [''] * (pedN - len(splt))
        if all(s == '' for s in splt):
            raise IndexError("No samples indicated.")
        
        sampleNum = []
        for s in splt:
            try:
                num = int(s)-1 if s else None
            except:
                raise ValueError("Illegal sample index: %s is not an integer." %s)
            if num is not None and (num < 0 or num >= len(self.filtus.files)):
                raise IndexError("Invalid sample index: %s" %s)
            if num in sampleNum:
                raise IndexError("Repeated sample index: %s" %s)
            sampleNum.append(num)
        
        return sampleNum
        
    def execute(self, button):
        if button is None or button=='Cancel':
            return self.deactivate()
        try:
            filtus = self.filtus
            filtus.busy()
            ped = self.pedtext.getvalue().strip()
            pedN = self.checkPed(ped)
            sampletxt = self.sampleIndex.get('0.0', 'end')
            samples = self.sampleIndexFix(sampletxt, pedN)
            what = self.whatButt.getvalue()
            dat, map, freq = 'dat' in what, 'map' in what, 'freq' in what
            maptype = self.mapCheck.getvalue()
            if map and not any(maptype):
                raise RuntimeError('At least one of the "Map options" must be checked')
            
            genmapfile = os.path.join(filtus.datadir, "DecodeMap_thin.txt") if 'cm' in maptype else None
            prefix = self.prefixEntry.getvalue()
            dir = self.dirBrowser.getvalue().strip()
            if not os.path.isdir(dir): raise os.error("Not a valid directory: %s" %dir)
            
            if 'freq' in what:
                frCol = self._altFreqMenu.getvalue()
                try:
                    def_freq = float(self._def_freq_entry.getvalue())
                    if not 0 < def_freq < 1: raise ValueError
                except ValueError:
                    raise ValueError("Please indicate a 'Missing entry value' strictly between 0 and 1. (This is used as the default allele frequency).")
            else:
                frCol, def_freq = None, None
                
            files = FiltusAnalysis.saveAsPed2(VFlist=filtus.filteredFiles, pedigree=ped, sampleIndex=samples, dat=dat, map=map, freq=freq, maptype = maptype, 
                        genmapfile=genmapfile, freqColumn=frCol, defaultFreq=def_freq, dir=dir, prefix=prefix)
            filtus.notbusy()
            FiltusUtils.infoMessage('Created the following files:\n\n' + '\n'.join(files))
        except Exception as e:
            filtus.notbusy()
            FiltusUtils.warningMessage("%s: %s"% (type(e).__name__, e))
            return
   

class AdvancedLoad(object):
    def __init__(self, filtus):
        parent = filtus.parent
        self.parent = parent
        self.filtus = filtus
        self.files = []
        self.dir = filtus.currentFileDir
        self.dialog = Pmw.Dialog(parent, buttons=('Load', 'Cancel'), title='Advanced load', command=self.doCommand, buttonbox_pady=10, 
                dialogchildsite_padx=0, dialogchildsite_pady=0, activatecommand=self.updateFilelist)
        self.dialog.withdraw()
        interior0 = self.dialog.interior()
        interior = Tkinter.Frame(interior0) #self.dialog.interior()
        interior.grid(row=0, column=0, padx=20, pady=10, sticky='news')
        HelpButton(interior0, filtus=filtus, page="loading", bookmark="advanced").grid(row=0, column=0, sticky="ne")
        
        interior.columnconfigure(0, weight=1)
        interior.rowconfigure(3, weight=1)
        
        # Title label
        Tkinter.Label(interior, text='Advanced file selector', font=filtus.titlefont).grid(row=0, padx=20, pady=10)
        dir_group = Pmw.Group(interior, tag_text = 'Directory')
        dir_interior = dir_group.interior()
        
        self.dirEntry = Pmw.EntryField(dir_interior, entry_width=20, entry_state='normal', value = self.dir)
        self.dirEntry.component('entry').bind("<Return>", self.updateFilelist)
        
        buttwidth=8 if filtus.windowingsystem == 'aqua' else 6
        dirBrowse = Tkinter.Button(dir_interior, command=self.askdir, text='Browse', 
                pady=0, padx=1, width=buttwidth, highlightthickness=0, font=filtus.smallfont)

        self.dirEntry.grid(sticky='we', padx=(10,0))
        dirBrowse.grid(row=0, column=1, sticky='w', padx=(5, 10), pady=5)

        ### options
        entry_OPTIONS = dict(labelpos='w', labelmargin=10, entry_width=25)
        
        option_group = Pmw.Group(interior, tag_text = 'Options')
        option_interior = option_group.interior()
        
        self.subdirVar = Tkinter.IntVar(parent, value=1)
        Tkinter.Checkbutton(option_interior, variable = self.subdirVar).grid(row=0, sticky='w', pady=5, padx=(10,0))
        self.subdirLabel = Tkinter.Label(option_interior, text = 'Include subdirectories')
        self.subdirLabel.grid(row=0, column=1, sticky='w', pady=5)

        self.endswithVar, self.containsVar, self.excludeVar = Tkinter.IntVar(parent), Tkinter.IntVar(parent), Tkinter.IntVar(parent)
        self.endswithEntry = Pmw.EntryField(option_interior, label_text='Only files ending in:', **entry_OPTIONS)
        self.containsEntry = Pmw.EntryField(option_interior, label_text='Only path/names containing:', **entry_OPTIONS)
        self.excludeEntry = Pmw.EntryField(option_interior, label_text='Exclude path/names containing:', **entry_OPTIONS)
        
        for i, str in enumerate(['endswith', 'contains', 'exclude']):
            intvar, entry = getattr(self, str+'Var'), getattr(self, str+'Entry')
            Tkinter.Checkbutton(option_interior, variable = intvar).grid(row=i + 1, sticky='w', padx=(10,0), pady=5)
            entry.component('entry').bind("<Return>", self.updateFilelist)
            entry.grid(row=i + 1, column=1, sticky='we', pady=5, padx=(0,10))

        Tkinter.Button(option_interior, text="Update file list", command=self.updateFilelist).grid(columnspan=2, pady=5)
        self.align()
         
        ### file list
        file_group = Pmw.Group(interior, tag_text = 'File list')
        file_interior = file_group.interior()
        file_interior.rowconfigure(0, weight=1)
        self.preview = Pmw.ScrolledListBox(file_interior, labelpos=None, selectioncommand=self.updateFileCount, 
                        listbox_height=6, listbox_activestyle='none', listbox_selectmode='extended', listbox_font=filtus.monofont)
        listbox = self.preview.component('listbox')
        listbox.bind('<Control-a>', self.selectAll)
        listbox.bind('<KeyRelease-Up>', self.updateFileCount)
        listbox.bind('<KeyRelease-Down>', self.updateFileCount)
        
        self.preview.grid(sticky='news', pady=(5,0), padx=10)
        self.previewSize = Tkinter.StringVar(parent)
        Tkinter.Label(file_interior, textvariable = self.previewSize).grid(row=1, column=0, sticky='nw', pady=(0,5), padx=10)
        
        Tkinter.Button(file_interior, text="Select all", command=self.selectAll).grid(row=1, column=0, sticky='ne', pady=(0, 5), padx=10)
        
        for g in [dir_group, option_group, file_group]:
            expandcol = 1 if g is option_group else 0 
            g.interior().columnconfigure(expandcol, weight=1)
            g.configure(tag_font = filtus.smallbold)
            g.grid(sticky='news', padx=10, pady=7)
        
        #self.updateFilelist()
        

    def align(self):
        Pmw.alignlabels([self.subdirLabel, self.endswithEntry, self.containsEntry, self.excludeEntry])

    def selectAll(self, event=None):
        self.preview.setvalue(self.preview.get())
        self.updateFileCount()
        
    def updateFilelist(self, event=None): # event unused, needed when function is used as callback
        newdir = self.dirEntry.getvalue().strip()
        if newdir != self.dir:
            self.dir = newdir
        subdir = self.subdirVar.get()
        endswith = self.endswithEntry.getvalue() if self.endswithVar.get() else None
        contains = self.containsEntry.getvalue() if self.containsVar.get() else None
        exclude = self.excludeEntry.getvalue() if self.excludeVar.get() else None
        
        self.files = self.getFiles(self.dir, subdir, endswith, contains, exclude)
        self.preview.setlist([' ' + os.path.relpath(filpath, self.dir) for fil, filpath in self.files])
        self.updateFileCount()

    def updateFileCount(self, event=None): # event unused, needed when function is used as callback
        n = len(self.preview.getvalue())
        self.previewSize.set("Files selected: %d" % n)
        
    def askdir(self):
        dir = tkFileDialog.askdirectory(initialdir=self.filtus.currentFileDir)
        if dir:
            self.filtus.currentFileDir = self.dir = os.path.normpath(dir)
            self.dirEntry.setvalue(self.dir)
            self.dirEntry.component('entry').xview('end')
            self.updateFilelist()

    def getFiles(self, dir, subdir, endswith, contains, exclude):
        if not dir or not os.path.isdir(dir): return []
        else: self.filtus.currentFileDir = dir
        if subdir:
            files = ((fil, os.path.join(dirpath, fil)) for dirpath, dirnames, files in os.walk(dir) for fil in files)
        else:
            files = ((fil, os.path.join(dir, fil)) for fil in os.listdir(dir))
        if endswith is not None:
            files = ((fil, filwpath) for fil, filwpath in files if FiltusUtils.myendswith(filwpath, endswith))
        if contains is not None:
            files = ((fil, filwpath) for fil, filwpath in files if FiltusUtils.contains(filwpath, contains))
        if exclude is not None:
            files = ((fil, filwpath) for fil, filwpath in files if FiltusUtils.not_contains(filwpath, exclude))

        files = [(fil, os.path.normpath(filpath)) for fil, filpath in files] # convert to list and fix slashes
        return sorted(files, key=itemgetter(1))

    def doCommand(self, result):
        if result != "Load":
            self.dialog.deactivate()
            return
        files = [os.path.join(self.dir, f.strip()) for f in self.preview.getvalue()]
        if not files: 
            FiltusUtils.warningMessage("No files selected.")
            return
        self.dialog.deactivate()
        self.filtus.loadFiles(files)
        

class FiltusText(Pmw.ScrolledText):
    def __init__(self, parent, filtus, **kw):
        defaultopts = dict(borderframe=5, columnheader=True, labelpos=None,
                scrollmargin=0, vscrollmode='static', hscrollmode='static', 
                text_width=50, text_height=20, text_wrap='none', text_font=filtus.monofont, 
                Header_font=filtus.monofont)
        defaultopts.update(kw)
        Pmw.ScrolledText.__init__(self, parent, **defaultopts)
            
        self.filtus = filtus
        self.textfield = self.component('text')
        self.columnheader = self.component('columnheader')
        self.currentColDat = None
        self.headers = None
        self.rightClick = None
        self.pad = 3
        self.meta = ''
        makeReadOnly(self.textfield)
        makeReadOnly(self.columnheader)

        self.sortMenu = Tkinter.Menu(filtus.parent, tearoff=0, font=filtus.defaultfont)
        self.sortMenu.add_command(label="Sort by this column, ascending", command=self._sort_up)
        self.sortMenu.add_command(label="Sort by this column, descending", command=self._sort_down)
        self.sortMenu.add_separator()
        self.sortMenu.add_command(label="Copy column to clipboard", command=self._selectColumn)
        self.sortMenu.add_separator()
        self.sortMenu.add_command(label="Description", command=self._showDescription)
        self.chosenColumn = None
        
        for event in filtus.rightclickevents:
            self.columnheader.bind(event, self._showSortMenu)
    
    def setlabel(self, label):
        self.configure(label_text=label)
        
    def _selectColumn(self):
        self.clipboard_clear()
        VF = self.currentColDat
        colGetter = VF.columnGetter(self.chosenColumn)
        if colGetter is not None:
            txt = '\n'.join(str(colGetter(v)) for v in VF.variants)
            self.clipboard_append(txt)
    
    def clearAll(self):
        self.clear()
        self.columnheader.delete('1.0', 'end')
        self.setlabel('')
        self.currentColDat = None
        self.headers = None
        self._setRightClickMenu(None)
        self.meta = ''
        self.pad = 3

    def _sort_up(self):
        self.currentColDat.sort(self.chosenColumn, False)
        self.prettyPrint(self.currentColDat, rightClick=self.rightClick)

    def _sort_down(self):
        self.currentColDat.sort(self.chosenColumn, True)
        self.prettyPrint(self.currentColDat, rightClick=self.rightClick)

    def _showDescription(self):
        descr = self.currentColDat.columnDescriptions.get(self.chosenColumn, "No description available")
        FiltusUtils.infoMessage('Description of column "%s":\n\n%s' %(self.chosenColumn, descr))

    def _showSortMenu(self, event):
        beg = self.columnheader.search('  ', "current", backwards=True, stopindex='1.0')
        wordstart = '1.0' if beg == '' else beg+'+2c'
        end = self.columnheader.search('  ', "current", backwards=False, stopindex='end')
        wordend = 'end' if end == '' else end
        column = self.columnheader.get(wordstart, wordend)
        if not self.headers or column not in self.headers:
            return
        self.chosenColumn = column
        descr_state = 'normal' if column in self.currentColDat.columnDescriptions else 'disabled'
        self.sortMenu.entryconfigure(5, state=descr_state)
        self.sortMenu.post(event.x_root, event.y_root)

    def firstword(self):
        line = self.textfield.get("current linestart", "current lineend")
        wordend = re.match('[^\\s\(;,]+', line)  # to first whitespace or any of ( ; ,
        if not wordend: return ''
        else: return wordend.group()

    def prettyPrint(self, ColDat, rightClick=None, meta=None, truncate=None, pad=None, label=None):
        self.currentColDat = current = ColDat.copyAttributes() # in case of sorting or other changes.
        self.headers = current.columnNames
        if meta: self.currentColDat.meta = meta
        if truncate is None: truncate=self.filtus.truncate
        if pad is None: pad = self.pad
        headerLine, body = current.printData(trunc=truncate, pad=pad)
        if label is not None: self.setlabel(label)
        self.columnheader.delete('1.0', 'end')
        self.columnheader.insert('1.0', headerLine)
        self.settext(body)
        self._setRightClickMenu(rightClick)
        self.pad = pad
    
    def _setRightClickMenu(self, menu):
        if menu is None: 
            menu = FiltusUtils.ignore_pass
        elif menu == "variantMenu":
            menu = self.variantViewMenu
        for event in self.filtus.rightclickevents:
            self.textfield.bind(event, menu)
        self.rightClick = menu
        
    def save(self):
        filtus = self.filtus
        if self.currentColDat is None and not self.getvalue().strip():
            FiltusUtils.warningMessage("No content to save")
            return
        filename = tkFileDialog.asksaveasfilename(initialdir=self.filtus.currentDir)
        if not filename:
            return
        filtus.currentDir = os.path.dirname(filename)
        try:
            if self.currentColDat:
                self.currentColDat.save(filename, sep = filtus.sepOutput, preamblePos=filtus.includePreamble)
            else:
                self.saveNonVF(filename, meta=self.meta)
        except Exception as e:
            FiltusUtils.warningMessage('%s\n\nFile not saved.' % e)
            return

    def saveNonVF(self, filename, meta):
        includePre = self.filtus.includePreamble
        sep = self.filtus.sepOutput
        with open(filename, 'w') as utfil:
            if meta and 'Top' in includePre: utfil.write(meta)
            utfil.write(self.getvalue())
            if meta and 'Bottom' in includePre: utfil.write('\n' + meta)
            
    def currentvariant(self, event):
        pos = self.textfield.index("@%d,%d" % (event.x, event.y))
        line,col = map(int, pos.split('.'))
        if col > 0 and pos == self.textfield.index("@0,%d" % event.y):
            raise IndexError
        return self.currentColDat.variants[line-1]
    
    def variantViewMenu(self, event):
        filtus = self.filtus
        try:
            currentVar = self.currentvariant(event)
        except IndexError:
            return
        vdef = self.currentColDat.varDefGetter(currentVar)
        
        def _show():
            if not hasattr(filtus, 'variantView'):
                filtus.variantView = VariantView(filtus, 
                    title="Unfiltered data from all samples")
            filtus.variantView.showVariant(vdef)
        
        menu = Tkinter.Menu(filtus.parent, tearoff = 0, font = filtus.defaultfont)
        menu.add_command(label="Display details for this variant", command=_show)
        menu.post(event.x_root, event.y_root)

        
class GeneView(Pmw.Dialog):
    def __init__(self, filtus, title):
        Pmw.Dialog.__init__(self, filtus.parent, title="Gene view", buttons=('Save to file', 'Collapse', 'Close window')) 
        self.withdraw()
        self.configure(command=self._doIt)
        self.filtus = filtus
        interior = self.interior()
        interior.columnconfigure(0, weight=1)
        interior.rowconfigure(1, weight=1)
        Tkinter.Label(interior, text=title, font=filtus.titlefont).grid(padx=20, pady=10)
        self.text = FiltusText(interior, filtus=filtus, labelpos='n')
        self.text.grid(padx=20, pady=20, sticky='news')
        
        self.collapseButton = self.component('buttonbox').button("Collapse")
        self.isCollapsed = None
        self.collapsedVF = None
        self.uncollapsedVF = None
        
    def _doIt(self, button):
        if button == "Save to file":
            self.text.save()
        elif button == "Collapse":
            butt = self.collapseButton
            if not self.isCollapsed:
                if self.collapsedVF is None:
                    self.collapsedVF = self.uncollapsedVF.collapse()
                self._set(self.collapsedVF)#, meta=self.meta)
                self.isCollapsed = True
                butt.configure(text='Uncollapse')
            else:
                self._set(self.uncollapsedVF)#, meta=self.meta)
                self.isCollapsed = False
                butt.configure(text='Collapse')
        else:
            self.deactivate()

    def _set(self, VF, meta=''):
        title = self._geneTitle(VF.genes)
        self.text.prettyPrint(VF, rightClick="variantMenu", label=title)
        
    def display(self, VF):#, meta=''):
        origVF = VF.copyAttributes()
        #self.meta = meta + '## GENE VIEW:\n## Viewing variants in: %s\n##\n' % ', '.join(origVF.genes)
        self._set(origVF)#, meta=self.meta)
        self.uncollapsedVF = origVF
        self.collapsedVF = None
        self.isCollapsed = False
        self.collapseButton.configure(text='Collapse')
        FiltusUtils.activateInCenter(self.filtus.parent, self)

    def _geneTitle(self, genes):
        if len(genes) == 1: 
            genetitle = genes[0]
        elif len(genes) < 15 and len(', '.join(genes)) < 50: 
            genetitle = ', '.join(genes)
        else: 
            genetitle = "%d selected genes"%len(genes)
        return genetitle
        
        
class VariantView(Pmw.Dialog):
    def __init__(self, filtus, title):
        Pmw.Dialog.__init__(self, filtus.parent, title="Variant view", buttons=('Save to file', 'Close window'))
        self.withdraw()
        self.configure(command=self._doIt)
        self.filtus = filtus
        interior = self.interior()
        interior.columnconfigure(0, weight=1)
        interior.rowconfigure(1, weight=1)
        Tkinter.Label(interior, text=title, font=filtus.titlefont).grid(padx=20, pady=10)
        self.text = FiltusText(self.interior(), filtus, labelpos='n')
        self.text.grid(padx=20, pady=20, sticky='news')
        
    def _doIt(self, button):
        if button == "Save to file": self.text.save()
        else: self.deactivate()
            
    def showVariant(self, vdef):
        variantDat = FiltusAnalysis.genotypeData(vdef, self.filtus.files)
        self.text.prettyPrint(variantDat, pad=5, label='Chromosome %s, position %s' %vdef)
        FiltusUtils.activateInCenter(self.filtus.parent, self)

        
class SharingPage(object):
    def __init__(self, notebook, filtus, title, manpage=None):
        self.notebook = notebook
        self.title = title
        self.page = notebook.add(title)
        #self.page.rowconfigure(0, weight=1)
        self.page.columnconfigure(0, weight=1)
        self.filtus = filtus
        self.sharingComputer = None
        
        self.interior = Tkinter.Frame(self.page)
        self.interior.columnconfigure(0, weight=1)
        
        self.buttonFrame = Tkinter.Frame(self.page)
        self.buttonFrame.columnconfigure(1, weight=1)
        self.button = Tkinter.Button(self.buttonFrame, text = "Analyze", command=self.analyze, pady=0)
        self.button.grid(padx=75, pady=(5, 0))
        
        self.summary = Pmw.ScrolledText(self.page, text_height=2, labelpos='nw', label_text='Summary:', text_font=filtus.monofont,
                     text_width=1, borderframe=0, scrollmargin=0, vscrollmode='none', hscrollmode='none', text_wrap='none', text_padx=2, text_pady=5)
        
        self.interior.grid(row=0, sticky='news')
        HelpButton(self.page, filtus=filtus, page=manpage).grid(row=0, sticky="ne")
        self.buttonFrame.grid(row=1, sticky='news', pady=5)
        self.summary.grid(row=2, padx=5, pady=(0, 5), sticky='news')

        makeReadOnly(self.summary.component('text'))

    def focus(self):
        self.notebook.selectpage(self.title)

    def clearAll(self):
        for w in self.fields:
            w.configure(entry_state='normal')
            w.clear()
        self.summary.clear()
        #self.sharingMachine.VFlist = None
        self.clear_more()
        
    def clear_more(self): # overridden in subclassses
        pass
        
    def readIDfield(self, entryfield, loadedFilenames):
        entry = entryfield.getvalue().strip()
        try:
            indices = FiltusUtils.convert2indices(entry, idlist=loadedFilenames)
        except ValueError as e:
            message = "Invalid entry: %s\n\nDetails:\n%s"%(entry, e)
            raise ValueError(message)
        return indices
        
    def confirmColumns(self, fields, gene=False, genotype=False, columnnames=[]):
        '''check if the samples indicated by f (a list of integer lists) have genotype columns'''
        VFlist = self.filtus.files
        for nr in [i for field in fields for i in field]:
            VF = VFlist[nr]
            if gene: VF.confirmGeneColumn()
            if genotype: VF.confirmGTcolumn()
            for col in columnnames:
                if not col in VF.columnNames:
                    raise ValueError("This analysis requires that the sample files have a column called %s. This is missing from:\n\n%s" %(col, VF.longName))
    
    def analyze(self):
        try:
            allFilteredSamples = self.filtus.checkLoadedSamples(select="all")
            input = self.validateEntries(allFilteredSamples)
            resultVF = self.sharingComputer.analyze(**input) 
            self.displayResult(resultVF)
        except ValueError as e:
            FiltusUtils.warningMessage(e)

            
class GeneSharingPage(SharingPage):
    def __init__(self, notebook, filtus, title, manpage=None, family=False):
        SharingPage.__init__(self, notebook, filtus, title, manpage)
        self.sharingComputer = FiltusAnalysis.GeneSharingComputer()
        self.family = family
        self.genelengths = {}
        
        if not family:
            self.tableMaker = FiltusAnalysis.NgTable(self, count_what = 'genes') # not to be used if family = False
            self.ngButton = Tkinter.Button(self.buttonFrame, text = "Table", command=self.tableMaker.table, pady=0)
            self.ngButton.grid(padx=70, row=0, column=1, pady=(5, 0))
            
        self.model = ModelSelect(self.interior)

        frame1 = Tkinter.Frame(self.interior)
        self.cases = Pmw.EntryField(frame1, entry_width=1, label_text = "Affected:", labelpos='w', labelmargin=5)
        self.controls = Pmw.EntryField(frame1, entry_width=1, label_text = "Healthy:", labelpos='w', labelmargin=5)
        self.shortcuts = Pmw.RadioSelect(frame1, buttontype = 'checkbutton', command=self.shortcutCommand, orient = 'vertical', pady=1)
        self.shortcuts.add('All affected')

        self.fields = [self.cases, self. controls]
        self.model.grid(padx=(5, 0), pady=(5, 0), sticky='nw')
        frame1.columnconfigure(0, weight=1)

        self.cases.grid(row=0, column=0, sticky='we', pady=2)
        self.controls.grid(row=1, column=0, sticky='we', pady=2)
        self.shortcuts.grid(padx=(30, 0), row=0, column=1, rowspan=2, pady=0, sticky='e')
        self.align()
        frame1.grid(padx=(5, 0), pady=0, sticky='news')

    def validateEntries(self, loadedVFs):
        loadedFilenames = [VF.longName for VF in loadedVFs]
        cases = self.readIDfield(self.cases, loadedFilenames)
        controls = self.readIDfield(self.controls, loadedFilenames)
        
        if len(cases) == 0: # if no 'cases' indicated
            raise ValueError("'Affected' field cannot be empty.")
        
        FiltusUtils.checkDuplicates([cases, controls])
        model = self.model.getvalue()
        self.confirmColumns([cases, controls], gene=True, genotype=(model!='Dominant'))
        VFcases=[loadedVFs[i] for i in cases]
        VFcontrols=[loadedVFs[i] for i in controls]
        
        # check if changes relevant for genelengths have been made, and if so update self.genelengths
        self.updateGenelengths()
        
        userInput = dict(VFcases=VFcases, VFcontrols=VFcontrols, model=model, family=self.family, genelengths=self.genelengths,
                         VFcases_index=cases, VFcontrols_index=controls)
        return userInput
    
    def updateGenelengths(self):
        genelengthBrowser = self.filtus.genelengthBrowser
        candidBrowser = self.filtus.FM.candidate_genes
        excluBrowser = self.filtus.FM.exclusion_genes
        
        # if no changes have been made, use previously generated data
        if not any(x.modified for x in (genelengthBrowser, candidBrowser, excluBrowser)):
            return self.genelengths
        
        # else, make from scratch
        genelengths = genelengthBrowser.getcontent()
        
        keep = set(genelengths)
        if candidBrowser.on():
            keep.intersection_update(candidBrowser.getcontent())
        if excluBrowser.on():
            keep.difference_update(excluBrowser.getcontent())
        
        self.genelengths = {gene:length for gene,length in genelengths.iteritems() if gene in keep}
        
     
    def rightClickMenu(self, event):
        menu = Tkinter.Menu(self.filtus.parent, tearoff = 0, font = self.filtus.defaultfont)
        text = self.filtus.text
        sharingResult = text.currentColDat
        gene = text.firstword()
        
        if gene in sharingResult.geneMaster:
            menu.add_command(label="Show relevant variants in %s" %gene, command=lambda: self.showVariantsInGenes(gene, choice=1))
            menu.add_command(label="Show relevant variants in %s and all genes above" %gene, command=lambda: self.showVariantsInGenes(gene, choice=2))
        menu.add_command(label="Show relevant variants in all genes", command=lambda: self.showVariantsInGenes(gene, choice=3))
        menu.post(event.x_root, event.y_root)

    def showVariantsInGenes(self, chosenGene, choice):
        filtus = self.filtus
        sharingResult = filtus.text.currentColDat
        if choice==1:
            genes = [chosenGene]
        elif choice==2:
            allgenes = [line[0] for line in sharingResult.variants]
            genes = allgenes[:allgenes.index(chosenGene)+1]
        elif choice==3:
            genes = [line[0] for line in sharingResult.variants]
        if len(genes) == 0: return
        
        varingens = sharingResult.variantsInGenes(genes)
        if not hasattr(self, 'sharingGeneView'):
            self.sharingGeneView = GeneView(filtus, title="Gene view of relevant variants in affected samples")
        self.sharingGeneView.display(varingens)
        
    def clear_more(self):
        self.model.setvalue('Dominant')
        self.shortcuts.setvalue([])
        
    def shortcutCommand(self, button, isTRUE):
        n = len(self.filtus.longFileNameList)
        ALL = FiltusUtils.intlist2string(range(1, n + 1)) if n > 0 else ''
        if isTRUE:
            self.cases.setvalue(ALL)
            self.controls.configure(entry_state='normal')
            self.controls.setvalue('')
            self.controls.configure(entry_state='disabled')
        else:
            self.controls.configure(entry_state='normal')

    def align(self):
        Pmw.alignlabels([self.cases, self.controls])

    def displayResult(self, resultVF):
        shareCounts = resultVF.shareCounts
        n = len(shareCounts)
        cumsum = [sum(shareCounts[i:]) for i in range(n)]
        self.summary.settext("Shared by    : " + "".join('%6d' %(i + 1,) for i in range(n)) + \
                           "\nGenes (cumul): " + "".join('%6d' %i for i in cumsum))    
        self.filtus.text.prettyPrint(resultVF, rightClick = self.rightClickMenu, label="Gene sharing results")
        
                                  
class VariantSharingPage(SharingPage):
    def __init__(self, notebook, filtus, title, manpage=None):
        SharingPage.__init__(self, notebook, filtus, title, manpage)
        self.sharingComputer = FiltusAnalysis.VariantSharingComputer()
        
        self.tableMaker = FiltusAnalysis.NgTable(self, count_what='variants')
        self.ngButton = Tkinter.Button(self.buttonFrame, text = "Table", command=self.tableMaker.table, pady=0)
        self.ngButton.grid(padx=70, row=0, column=1, pady=(5, 0))
        
        page = self.interior
        frame1 = Tkinter.Frame(page)

        Tkinter.Label(frame1, text = 'Allele distribution:').grid(padx=5, pady=10, sticky='nw')
        self.shortcuts = Pmw.RadioSelect(frame1, buttontype = 'checkbutton', command=self.shortcutCommand, orient = 'horizontal', pady=1)
        self.shortcuts.add('All aff dom')
        self.shortcuts.add('All aff rec')

        frame1.columnconfigure(0, weight=1)
        frame1.grid(padx=(5, 0), pady=0, sticky='new')
        self.shortcuts.grid(row=0, column=1, sticky='e')

        self.allelFrame = Tkinter.Frame(page)
        wd = 9
        self.allelFrame.grid(padx=5, pady=0, sticky='nw')
        self.allel0 = Pmw.EntryField(self.allelFrame, entry_width=wd, label_text = "0:", labelpos='w', labelmargin=5)
        self.allel01 = Pmw.EntryField(self.allelFrame, entry_width=wd, label_text = " 0 or 1:", labelpos='w', labelmargin=5)
        self.allel1 = Pmw.EntryField(self.allelFrame, entry_width=wd, label_text = "1:", labelpos='w', labelmargin=5)
        self.allel12 = Pmw.EntryField(self.allelFrame, entry_width=wd, label_text = " 1 or 2:", labelpos='w', labelmargin=5)
        self.allel2 = Pmw.EntryField(self.allelFrame, entry_width=wd, label_text = "2:", labelpos='w', labelmargin=5)
        self.allel0.grid(row=0, column=0, sticky='nw', padx=0)
        self.allel1.grid(row=0, column=1, sticky='nw', padx=7)
        self.allel2.grid(row=0, column=2, sticky='nw', padx=7)
        self.allel01.grid(row=1, column=0, sticky='nw', padx=0)
        self.allel12.grid(row=1, column=1, sticky='nw', padx=7)
        self.fields = [self.allel0, self.allel1, self.allel2, self.allel01, self.allel12]
        self.align()

    def validateEntries(self, loadedVFs):
        loadedFilenames = [VF.longName for VF in loadedVFs]
        fields = [self.readIDfield(w, loadedFilenames) for w in self.fields]
        if len(fields[1]+fields[2]+fields[4]) == 0: # empty 1, 2, 1-2 fields
            raise ValueError("Necessary field(s) empty")
        FiltusUtils.checkDuplicates(fields)
        self.confirmColumns(fields, genotype=True)
        VF_inputs = [[loadedVFs[i] for i in field] for field in fields]
        userInput = dict(VF_inputs=VF_inputs, VF_index_list=fields)
        return userInput
    
    def clear_more(self):
        self.shortcuts.setvalue([])
        
    def shortcutCommand(self, button, isTRUE):
        wlist = set(self.fields)
        if isTRUE:
            if button == 'All aff dom':
                wsel = self.allel12
                self.shortcuts.setvalue(['All aff dom'])
            else:
                wsel = self.allel2
                self.shortcuts.setvalue(['All aff rec'])
            K = len(self.filtus.longFileNameList)
            if K == 0: return
            if K == 1: txt = '1'
            else: txt = '1-%d' %K
            wsel.configure(entry_state='normal')
            wsel.setvalue(txt)
            for w in wlist.difference({wsel}):
                w.configure(entry_state='normal')
                w.clear()
                w.configure(entry_state='disabled')
        else:
            for w in wlist:
                w.configure(entry_state='normal')

    def getentry(self, field):
        return FiltusUtils.convert2indices(self.fields[field].getvalue(), idlist = self.filtus.longFileNameList)

    def align(self):
        Pmw.alignlabels([self.allel0, self.allel01], sticky='e')
        Pmw.alignlabels([self.allel1, self.allel12], sticky='e')

    def displayResult(self, resultVF):
        self.filtus.text.prettyPrint(resultVF, rightClick = None, label="Variant sharing results")
        self.summary.settext("Variants found: %d" % resultVF.length)
        
    
class VCFgenotypeDisplay(Pmw.Dialog):
    def __init__(self, filtus, triodic, vdef, denovoprob):
        VFch, VFfa, VFmo = [filtus.filteredFiles[triodic[x]] for x in ['child', 'father', 'mother']]
        ch_var, fa_var, mo_var = [VF.getVar(vdef, firstHit=True) for VF in (VFch, VFfa, VFmo)]
        columnNames = VFch.formatHeads[:]
        ch_vals = VFch.columnGetter(*columnNames)(ch_var)
        fa_vals = VFfa.columnGetter(*columnNames)(fa_var)
        mo_vals = VFmo.columnGetter(*columnNames)(mo_var)
        
        # expand/explain field names
        descr = dict(GT="Genotype", AD="Allele depths\n0, 1", DP="Read depth", GQ="Genotype\nquality", 
                     PL="Phred likelihoods\n0/0,  0/1,  1/1")
        for key,val in descr.iteritems():
            if key in columnNames: 
                columnNames[columnNames.index(key)] = val
        
        # add space after commas
        ch_vals, fa_vals, mo_vals = [[x.replace(',', ', ') for x in vals] for vals in ch_vals, fa_vals, mo_vals] 
                
        Pmw.Dialog.__init__(self, filtus.parent, buttons = ('OK',), title = 'VCF genotype display', defaultbutton=0)
        self.withdraw()
        interior = self.interior()
        
        # Title label
        Tkinter.Label(interior, text='Genotype details for variant at %s:%s'%vdef[:2], font=filtus.titlefont).grid(padx=20, pady=10)
        
        alleles = 'REF allele (0) = %s\nALT allele (1) = %s' % vdef[2:]
        Tkinter.Label(interior, text=alleles).grid(padx=20, pady=(5,10))
        
        tableframe = Tkinter.Frame(interior)
        margin_labelops = dict(bd=2, relief="raised")
        value_labelops = dict(bg="white", bd=1, relief="groove")
        grid_ops = dict(sticky='news', padx=0, pady=0)
        for i, header in enumerate(columnNames):
            Tkinter.Label(tableframe, text = header, width=17 if i==4 else 11, **margin_labelops).grid(row=0, column=i+1, **grid_ops)
        for row, id in enumerate(['Child', 'Father', 'Mother']):
            Tkinter.Label(tableframe, text = id, **margin_labelops).grid(row=row+1, column=0, ipady=3,  **grid_ops)
        for row, values in enumerate([ch_vals, fa_vals, mo_vals]):   
            for col, val in enumerate(values):
                Tkinter.Label(tableframe, text = str(val), **value_labelops).grid(row=row+1, column=col+1, **grid_ops)
        tableframe.grid(padx=10, pady=10)
        
        denovo = 'Estimated de novo probability = %s%%' %(float(denovoprob)*100,)
        Tkinter.Label(interior, text=denovo).grid(padx=20, pady=(20,30))
        

class LabeledListBox(Pmw.MegaWidget):
    def __init__(self, parent, filtus, **kw):
        # Define the megawidget options.
        optiondefs = (
            ('toptext', '', None),
            ('bottomtext', '', None),
            ('height', 6, Pmw.INITOPT),
            ('width', 36, Pmw.INITOPT),
            ('selectmode', 'extended', Pmw.INITOPT)
        )
        self.defineoptions(kw, optiondefs)

        # Initialise base class (after defining options).
        Pmw.MegaWidget.__init__(self, parent)
        self.filtus = filtus
        self._topString = Tkinter.StringVar(self.winfo_toplevel(), value = self['toptext'])
        self._bottomString = Tkinter.StringVar(self.winfo_toplevel(), value = self['bottomtext'])

        #create components
        interior = self.interior()

        self._scrolledlistbox = self.createcomponent('scrolledlistbox', (('listbox', 'scrolledlistbox_listbox'),), None,
                                Pmw.ScrolledListBox, (interior,), listbox_width=self['width'], listbox_height=self['height'],
                                listbox_activestyle = 'none', labelpos='nw', label_textvariable = self._topString,
                                listbox_selectmode = self['selectmode'],
                                listbox_font = filtus.monofont)
        self._bottomlabel = self.createcomponent('bottomlabel', (), None, Tkinter.Label, (interior,),
                            textvariable = self._bottomString, font = filtus.monofont)
        self.initialiseoptions()

        interior.columnconfigure(0, weight=1)
        interior.rowconfigure(0, weight=1)
        self._scrolledlistbox.grid(sticky='news')
        self._bottomlabel.grid(sticky='w')

        
    def fixselectall(self):
        self.component('scrolledlistbox_listbox').bind('<Control-a>', self.selectall)

    def selectall(self, event=None):
        self.component('scrolledlistbox_listbox').select_set(0, self._scrolledlistbox.size()-1)

    def settoptext(self, text):
        self._topString.set(text)

    def setbottomtext(self, text):
        self._bottomString.set(text)

    def size(self):
        return self._scrolledlistbox.size()

    def clearAll(self):
        self._scrolledlistbox.clear()
        self.setbottomtext('')

    def setAll(self, list, enum=True): # A bit sloppy: This is intended for fileListbox instances, while the SummaryBox subclass has its own version.
        if enum:
            list = ['%2d: %s' %(i + 1, filename) for i, filename in enumerate(list)]
        self.setlist(list)
        self.settoptext('Loaded samples: %d' %len(list))

    def newline(self, text): # adds one line to the listbox and makes it visible. (Used in applyFilters())
        self.insert('end', text)
        self.see('end')
        self.update()
        
Pmw.forwardmethods(LabeledListBox, Pmw.ScrolledListBox, '_scrolledlistbox')


class SummaryBox(LabeledListBox):
    def __init__(self, parent, **kw):
        kw['selectmode'] = 'single'
        kw['scrolledlistbox_hscrollmode'] = 'none'
        kw['scrolledlistbox_dblclickcommand'] = self.showvariants
        apply(LabeledListBox.__init__, (self, parent), kw)
        self.VFlist = []

    def clearAll(self):
        self._scrolledlistbox.clear()
        self.setbottomtext('')
        self.VFlist = []

    def addVF(self, VF):
        nV = VF.length
        if VF.varDefGetter:
            name = 'variant ' if nV == 1 else 'variants'
        else:
            name = 'row ' if nV == 1 else 'rows'
        summary = '%2d: %5d %s' %(self.size() + 1, nV, name)
        if VF.varDefGetter and VF.geneGetter:
            nG = VF.nGenes
            summary += ' in %5d gene%s' %(nG, 's' if nG != 1 else '')
        self.newline(summary)
        self.VFlist.append(VF)

    def setAll(self, list):
        self.clearAll()
        for VF in list:
            self.addVF(VF)
        self.setTotal()

    def setTotal(self):
        if any(VF.varDefGetter is None for VF in self.VFlist):
            return
        tot = len(set([VF.varDefGetter(v) for VF in self.VFlist for v in VF.variants]))
        summary = '%d unique variant%s' % (tot, 's' if tot != 1 else ' ')
        
        VFlistG = [VF for VF in self.VFlist if VF.geneGetter]
        geneGetters = [VF.annotatedGenes for VF in VFlistG] #saving time
        ngenes = len(set(g for VF, geneGetter in zip(VFlistG, geneGetters) for v in VF.variants for g in geneGetter(v)))
        if len(VFlistG)>0:
            summary += ', %d gene%s' %(ngenes, 's' if ngenes != 1 else '')
        self.setbottomtext('Total: ' + summary)
        gc.collect()

    def showvariants(self):
        if len(self.VFlist) == 0 or len(self.curselection())==0: return
        selec = int(self.curselection()[0])
        VF = self.VFlist[selec]
        filtus = self.filtus
        filtus.busy()
        meta = FiltusUtils.preambleNY(VFlist=[VF], sort=False)
        rightClick = "variantMenu" if VF.isVCFtype else None
        labelstart =  'Unfiltered variants: ' if self._topString.get().startswith("Unfilt") else 'Variants after filtering: '
        filtus.text.prettyPrint(VF, rightClick=rightClick, meta=meta, label=labelstart + VF.shortName) # OK to pass on meta in this case
        filtus.notbusy()
        
    
class FilterMachine(Pmw.MegaWidget):
    def __init__(self, parent, filtus, manpage, **kw):
        # Define the megawidget options.
        optiondefs = ()
        self.defineoptions(kw, optiondefs)

        # Initialise base class (after defining options).
        Pmw.MegaWidget.__init__(self, parent)

        self.colnames = [''] + FiltusUtils.listUnique([head for VF in filtus.files for head in VF.columnNames])
        self.operatorNames = ["equal to", "not equal to", "starts with", "does not start with", "contains", "does not contain", "greater than", "less than"]
        self.filtus = filtus
        self.CFlist = []
        self.nColumnFilters = None

        #create components
        group = self.createcomponent('group', (), None, Pmw.Group, self.interior(), tag_text = 'Filters')
        groupi = group.interior()
        common_OPTIONS = dict(filtus=filtus, checkbutton = True, state='disabled')

        self.candidate_genes = self.createcomponent('candidategenes', (), None, GeneFile, groupi, browsetitle = 'Select file with candidate genes', label = "Restrict to genes:", **common_OPTIONS)
        self.exclusion_genes = self.createcomponent('excludegenes', (), None, GeneFile, groupi, browsetitle = "Select file with exclusion genes", label = "Exclude genes:", **common_OPTIONS)
        self.exclusion_variants = self.createcomponent('excludevariants', (), None, VariantFile, groupi, label = "Exclude variants:", **common_OPTIONS)
        self.regions = self.createcomponent('regions', (), None, RegionFile, groupi, browsetitle = 'Select file with genomic regions', label = "Restrict to regions:", **common_OPTIONS)
        self.main_filter_widgs = (self.candidate_genes, self.exclusion_genes, self.exclusion_variants, self.regions)

        cftitle = Tkinter.Frame(groupi)
        Tkinter.Label(cftitle, text = "Column filters:").grid(row=0, sticky='w', pady=(0, 2))
        Tkinter.Label(cftitle, text = "keep if\nmissing", font = filtus.tinyfont, pady=0, bd=0).grid(row=0, column=1, sticky='es')
        self.columnFrame = self.createcomponent('columnframe', (), None, Tkinter.Frame, groupi)
        self.applyButton = Tkinter.Button(groupi, text = "Apply filters", command=self.applyFilters, pady=0)
        self.align()
        for w in (self.interior(), groupi, cftitle, self.columnFrame):
            w.columnconfigure(0, weight=1)
        group.grid(sticky='news')
        self.candidate_genes.grid(row=0, padx=(10, 5), pady=(0, 3), sticky='ew')
        self.exclusion_genes.grid(row=1, padx=(10, 5), pady=3, sticky='ew')
        self.exclusion_variants.grid(row=2, padx=(10, 5), pady=3, sticky='ew')
        self.regions.grid(row=3, padx=(10, 5), pady=3, sticky='ew')
        cftitle.grid(row=4, sticky='ew', padx=(10, 2))
        self.columnFrame.grid(row=5, sticky='ew', pady=(0, 10))
        self.applyButton.grid(row=6, sticky='n', pady=(0, 10))
        
        HelpButton(groupi, filtus=filtus, page=manpage).grid(row=6, sticky="es")
        self.initialiseoptions()

        self._initColFilters(k=3)

    def align(self):
        Pmw.alignlabels([w.component('entryfield') for w  in self.main_filter_widgs])

    def _initColFilters(self, k):
        for w in self.CFlist: w.destroy()
        self.CFlist = [ColumnFilter(parent=self.columnFrame, filtus = self.filtus, colnames = self.colnames, operatorNames = self.operatorNames) for i in range(k)]
        for CF in self.CFlist:
            CF.grid(sticky='new')
        self.nColumnFilters = k

    def setColnames(self, colnames):
        self.colnames = [''] + colnames
        for CF in self.CFlist:
            CF.setcolnames(self.colnames)
    
    def applyFilters(self):
        filtus = self.filtus
        if len(filtus.files) == 0: return
        try:
            filter = self.getFilter(model=None)
        except Exception as e:
            FiltusUtils.warningMessage(e)
            return
        summary2 = filtus.fileSummary2
        summary2.clearAll()
        
        filtus.busy()
        filtus.filteredFiles = []
        for VF in filtus.files:
            VFfilt = filter.apply(VF)
            filtus.filteredFiles.append(VFfilt)
            summary2.addVF(VFfilt)
        summary2.setTotal()
        filtus.notbusy()
        
    def getFilter(self, model = None, closePairLimit=0):
        cfs = [CF.getfilter() for CF in self.CFlist]
        cfs = [cf for cf in cfs if cf] # remove empty (i.e. None) filters

        kwargs = dict(filtus=self.filtus, model = model, closePairLimit = closePairLimit, columnfilters=cfs)

        if self.candidate_genes.on():     kwargs['restrict_genes_txt'] = self.candidate_genes.getvalue().strip()
        if self.exclusion_genes.on():     kwargs['exclude_genes_txt'] = self.exclusion_genes.getvalue().strip()
        if self.exclusion_variants.on():  kwargs['exclude_var_txt'] = self.exclusion_variants.getvalue().strip()
        if self.regions.on():             kwargs['regions_txt'] = self.regions.getvalue().strip()

        return Filter.Filter(**kwargs)

    def addColumnFilter(self):
        self.CFlist.append(ColumnFilter(parent=self.component('columnframe'), filtus = self.filtus, 
                           colnames = self.colnames, operatorNames = self.operatorNames))
        self.CFlist[-1].grid(sticky='new')
        self.nColumnFilters += 1
        self.filtus.scrollframe.reposition()
        
    def removeColumnFilter(self):
        if self.nColumnFilters <= 1:
            # filtus.mainMenu.item configure (Remove col filter state='disabled')
            return
        else:
            cf = self.CFlist.pop() #removes last item
            cf.destroy()
            self.nColumnFilters -= 1
            self.filtus.scrollframe.reposition()
        
    def unapplyFilters(self):
        self.candidate_genes.deselect()
        self.exclusion_genes.deselect()
        self.exclusion_variants.deselect()
        self.regions.deselect()
        self.filtus.fileSummary2.clearAll()
        self.filtus.filteredFiles = self.filtus.filteredFiles_initialcopy()

    def clearFilters(self):
        self.candidate_genes.deselect()
        self.exclusion_genes.deselect()
        self.exclusion_variants.deselect()
        self.regions.deselect()
        self._initColFilters(self.nColumnFilters)
        self.filtus.fileSummary2.clearAll()
        self.filtus.filteredFiles = self.filtus.filteredFiles_initialcopy()

    def loadFilterBrowse(self):
        fconfigfile = tkFileDialog.askopenfilename(initialdir=self.filtus.currentDir, title = "Select filter configuration file")#, filetypes = [('All files', '.*'), ("filter config files",".fconfig")])
        if not fconfigfile:
            return
        self.filtus.currentDir = os.path.dirname(fconfigfile)
        self.loadFilter(fconfigfile)

    def loadFilter(self, fconfigfile=None, candidate_genes=None, exclusion_genes=None, exclusion_variants=None, regions=None, columnfilters=None):
        self.clearFilters()
        cur = self.nColumnFilters
        
        if fconfigfile is not None:
            try:
                with open(fconfigfile, 'rU') as f:
                    lines = [line.strip() for line in f if line.strip()]
                
                if len(lines) < 6:
                    raise IndexError("Too few lines in filter configuration file")
                _candidate_genes = lines[0]
                _exclusion_genes = lines[1]
                _exclusion_variants = lines[2]
                
                new = len(lines) >= 5 and len(lines[4]) < 3
                _regions = lines[3] if new else '0'
                
                ncf = int(lines[4]) if new else int(lines[3])
                if 4 + new + ncf != len(lines):
                    raise ValueError("Wrong number of column filters:\nNumber indicated is %d, but the file contains %d" % (ncf, len(lines)-4-new))
                _columnfilters = lines[(4+new):(4+new+ncf)]
            except Exception as e:
                FiltusUtils.warningMessage("Cannot read filter configuration file.\n\n%s" % e)
                self.clearFilters()
                return   
        else:    
            _candidate_genes = '0' if candidate_genes is None else '1 %s' % candidate_genes
            _exclusion_genes = '0' if exclusion_genes is None else '1 %s' % exclusion_genes
            _exclusion_variants = '0' if exclusion_variants is None else '1 %s' % exclusion_variants
            _regions = '0' if regions is None else '1 %s' % regions
            _columnfilters = [] if columnfilters is None else columnfilters
            
        try:
            self.candidate_genes.setconfig(_candidate_genes)
            self.exclusion_genes.setconfig(_exclusion_genes)
            self.exclusion_variants.setconfig(_exclusion_variants)
            self.regions.setconfig(_regions)
            
            self._initColFilters(max(cur, len(_columnfilters)))
            try:
                for i, cf in enumerate(_columnfilters):
                    self.CFlist[i].setfilter(cf)
            except KeyError as rel:
                raise KeyError("Operator '%s' is unknown."%rel)
            
        except Exception as e:
            FiltusUtils.warningMessage("Cannot load filter.\n\n%s" % e)
            return
        self.applyButton.focus()

    def saveFilter(self):
        fconfigname = tkFileDialog.asksaveasfilename(initialdir=self.filtus.currentDir, title = "Save filter configuration as", defaultextension = 'fconfig')
        if not fconfigname:
            return
        self.filtus.currentDir = os.path.dirname(fconfigname)
        with open(fconfigname, 'w') as f:
            f.write('\n'.join(map(str, (self.candidate_genes, self.exclusion_genes, self.exclusion_variants, self.regions))) + '\n')
            CFs = [str(CF) for CF in self.CFlist if not CF.isEmpty()]
            f.write('%d\n' % len(CFs) + '\n'.join(CFs))

            
class MultiFilter(object):
    def __init__(self, filtus):
        self.filtus = filtus
        self.parent = filtus.parent
        self.dialog = Pmw.Dialog(filtus.parent, buttons=('Apply filters', 'Cancel'), 
            title='Individual filtering', command=self.execute, defaultbutton=0)
        self.dialog.withdraw()
        interior = Tkinter.Frame(self.dialog.interior())
        self.dialog.interior().columnconfigure(0, weight=1)
        interior.columnconfigure(0, weight=1)
        interior.grid(sticky='news', padx=10, pady=10)

        pmw_OPTIONS = dict(labelmargin=10, labelpos='w')
        grid_OPTIONS = dict(sticky='news', padx=15, pady=15)

        self.filter_entries, self.samples = [], []
        for i in range(4):
            filt = FileBrowser(interior, filtus=filtus,
                label="Filter:", checkbutton = False, labelpos='w', browsesticky='se',
                entryfield_command=None, entryfield_entry_width=32,
                browsetitle="Select filter configuration file")
            samp = Pmw.EntryField(interior, labelpos='w', labelmargin=5,
                label_text="Apply to samples:", entry_width=10)
            filt.grid(row=i, column=0, **grid_OPTIONS)
            samp.grid(row=i, column=1, **grid_OPTIONS)
            self.filter_entries.append(filt)
            self.samples.append(samp)
            
    def showDialog(self):
        FiltusUtils.activateInCenter(self.parent, self.dialog)

    def execute(self, button):
        filtus = self.filtus
        nFiles = len(filtus.files)
        if nFiles == 0 or button is None or button == 'Cancel':
            self.dialog.deactivate()
            return
        filter_txt = [entr.getvalue() for entr in self.filter_entries]
        sample_txt = [entr.getvalue() for entr in self.samples]
        try:
            filters, samples = self.validate(filter_txt, sample_txt)
        except Exception as e:
            FiltusUtils.warningMessage(e)
            return
            
        sample_indices = [next((j for j,s in enumerate(samples) if i in s), None) for i in range(nFiles)]
        filters_ordered = [filters[j] if j is not None else None for j in sample_indices]
        self.dialog.deactivate()
        summary2 = filtus.fileSummary2
        summary2.clearAll()
        
        filtus.filteredFiles = []
        for VF, filter in zip(filtus.files, filters_ordered):
            if filter:
                VFfilt = filter.apply(VF)
            else:
                VFfilt = VF.copyAttributes()
            filtus.filteredFiles.append(VFfilt)
            summary2.addVF(VFfilt)
        summary2.setTotal()
        
    
    def validate(self, filter_txt, sample_txt):
        filters, samples = [],[]
        filenames = self.filtus.longFileNameList
        for f, s in zip(filter_txt, sample_txt):
            if f:
                try:
                    filters.append(Filter.Filter(filterFile=f, filtus=self.filtus))
                except Exception as e:
                    raise IOError("Could not read filter file '%s'.\n\n%s" % (f, e))
                try:
                    samples.append(FiltusUtils.convert2indices(s, idlist=filenames))
                except ValueError as e:
                    raise ValueError("Invalid entry: %s\n\n%s" % (s, e))
            elif s:
                raise ValueError("No filter indicated for samples %s" % s)
                
        if all(len(s)==0 for s in samples):
            raise ValueError("No samples indicated")

        FiltusUtils.checkDuplicates(samples)
        return filters, samples
        

class ColumnFilter(Pmw.MegaWidget):
    def __init__(self, parent, filtus, **kw):
        # Define the megawidget options.
        optiondefs = (
            ('colnames', [], None),
            ('operatorNames', [], None),
        )
        self.defineoptions(kw, optiondefs)

        # Initialise base class (after defining options).
        Pmw.MegaWidget.__init__(self, parent)

        self.colnames = self['colnames']
        self.operatorNames = self['operatorNames']

        #create components
        interior = self.interior()
        width=12 # if filtus.windowingsystem != 'aqua' else 12
        button_OPTIONS = dict(menubutton_anchor = 'w', menubutton_padx=5, menubutton_pady=0, menubutton_width=width, menu_font = filtus.defaultfont)
        self.keepMissingVar = Tkinter.IntVar(self.winfo_toplevel(), value=0)

        self.columnOM = self.createcomponent('optionmenu1', (), None, OptionMenuExt, (interior,), **button_OPTIONS)
        self.relationOM = self.createcomponent('optionmenu2', (), None, Pmw.OptionMenu, (interior,), items = ['']+self.operatorNames, initialitem = '', **button_OPTIONS)
        self.relEntry = self.createcomponent('entryfield', (), None, Pmw.EntryField, (interior,), entry_width=15)
        self.checkbutton = self.createcomponent('checkbutton', (), None, Tkinter.Checkbutton, (interior,), variable = self.keepMissingVar)
        self.initialiseoptions()

        interior.columnconfigure(2, weight=1)
        self.columnOM.grid(row=0, sticky='w', padx=(5, 0))
        self.relationOM.grid(row=0, column=1, sticky='w', padx=(5, 0))
        self.relEntry.grid(row=0, column=2, sticky='ew', padx=(5, 0), pady=4)
        self.checkbutton.grid(row=0, column = 3, sticky='e', padx=(10, 0), pady=0)

        self.setcolnames(self.colnames)

    def __str__(self):
        return ' ::: '.join([self.getcol(), self.getrelation(), self.getvalue(), str(self.keepMissing())])

    def getcol(self):
        return self.columnOM.getvalue()

    def getrelation(self):
        return self.relationOM.getvalue()

    def getvalue(self):
        return self.relEntry.getvalue()

    def keepMissing(self):
        return self.keepMissingVar.get()

    def getfilter(self):
        if self.isEmpty():
            return None
        col, rel, val, keep = self.getcol(), self.getrelation(), self.getvalue(), self.keepMissingVar.get()
        if ('less' in rel or 'greater' in rel):
            try:
                val = float(val)
            except ValueError:
                FiltusUtils.warningMessage("Column filter ignored:\n\n'%s  %s  %s'\n\nNumerical value needed."%(col, rel, val))
                return None
        return (col, rel, val, keep)

    def setfilter(self, s):
        col, rel, entry, kim = s.strip().split(' ::: ')
        self.columnOM.setAndCheck(col.strip())
        self.relEntry.setvalue(entry)
        self.keepMissingVar.set(int(kim))
        rel = rel.strip()
        if rel in self.operatorNames:
            self.relationOM.setvalue(rel)
        else:
            raise KeyError(rel)

    def isEmpty(self):
        return bool(self.getcol() == '' or self.getrelation() == '')

    def setcolnames(self, items):
        self.colnames = items
        self.columnOM.setItems(items)


class OptionMenuExt(Pmw.OptionMenu):
    def __init__(self, parent, **kw):
        if not 'command' in kw:
            kw['command'] = self.check
        Pmw.OptionMenu.__init__(self, parent, **kw)

        self.items = kw.get('items', [])
        self.okColor = parent.cget('bg')
        self.errorColor = 'orange'
        self.inconsistent = False

    def setAndCheck(self, value):
        if value is None: value=''
        self.setvalue(value)
        self.check(value)

    def setItems(self, items):
        current = self.getvalue()
        self.items = items
        self.setitems(items)
        self.setAndCheck(current)
        if items:
            nCols = -(-len(items)/28)
            N = -(-len(items)/nCols)
            for i in range(1, nCols):
                self.component('menu').entryconfigure(i * N, columnbreak=1)

    def check(self, value):
        test = value in self.items
        self.setColor(test)

    def setColor(self, test):
        self.inconsistent = not test
        color = self.okColor if test else self.errorColor
        self.configure(menubutton_bg = color, menubutton_activebackground=color)

    def disable(self):
        self.component("menubutton").configure(state = 'disabled')

    def enable(self):
        self.component("menubutton").configure(state = 'normal')

        
class FileBrowser(Pmw.MegaWidget):
    '''Megawidget containing an entryfield, a 'browse' button, and possibly a checkbutton.
    '''
    def __init__(self, parent, filtus, **kw):

        # Define the megawidget options.
        optiondefs = (
            ('browsetitle', None, Pmw.INITOPT),
            ('filetypes', [('All files', '.*')], Pmw.INITOPT),
            ('checkbutton', True, Pmw.INITOPT),
            ('label', '', None),
            ('labelpos', 'w', Pmw.INITOPT),
            ('labelmargin', 5, Pmw.INITOPT),
            ('value', '', Pmw.INITOPT),
            ('state', 'normal', None),
            ('browsesticky', 'ew', Pmw.INITOPT)
        )
        self.defineoptions(kw, optiondefs)

        # Initialise base class (after defining options).
        Pmw.MegaWidget.__init__(self, parent)
        self.filtus = filtus
        self.content = set()
        self.modified = 0 #indicator variable, used in entryfield.modifiedcommand

        # Create the components.
        interior = self.interior()

        if self['checkbutton']:
            self.checkVar = Tkinter.IntVar(self.winfo_toplevel(), value=0)
            self.checkbutton = self.createcomponent('checkbutton', (), None, Tkinter.Checkbutton, (interior,),
                                        variable=self.checkVar, command=self._togglestatus, pady=0)
            self.checkbutton.grid(row=0, column=0, sticky='w')
            entrycol = 1
            browsecol = 2
        else:
            entrycol = 0
            browsecol = 1

        interior.columnconfigure(entrycol, weight=1)

        self.entryfield = self.createcomponent('entryfield', (), None, Pmw.EntryField, (interior,), entry_width=1,
                                    label_text = self['label'], labelpos=self['labelpos'], command=self.setcontent, labelmargin = self['labelmargin'],
                                    value = self['value'], entry_state=self['state'], modifiedcommand=self.modifiedIsTrue)
        buttwidth = 8 if filtus.windowingsystem == 'aqua' else 6
        self.browsebutton = self.createcomponent('button', (), None, Tkinter.Button, (interior,),
                                    command=self._browseCommand, text = 'Browse', width=buttwidth, pady=0, padx=1, highlightthickness = 0, font = filtus.smallfont)

        self.entryfield.grid(row=0, column = entrycol, sticky='ew', pady=0)
        self.browsebutton.grid(row=0, column = browsecol, sticky=self['browsesticky'], padx=(5, 0), pady=0)

        self.initialiseoptions()

    def __str__(self):
        return '%d %s' %(self.on(), self.getvalue())

    def modifiedIsTrue(self):
        self.modified = 1

    def setconfig(self, s): #s a string produced by __str__
        if s[0] == '0':
            self.deselect()
        elif s[0] == '1':
            self.select()
        else:
            raise ValueError("Expected string with initial character 0 or 1, but recieved %s" % s)
        path = os.path.normpath(s[2:]) if len(s)>2 else '' # because normpath('') gives '.'
        self.setvalue(path)

    def select(self):
        if self['checkbutton']:
            self.checkbutton.select()
            self._togglestatus()

    def deselect(self):
        if self['checkbutton']:
            self.checkbutton.deselect()
            self._togglestatus()

    def on(self):
        if self['checkbutton']:
            return self.checkVar.get()

    def setcontent(self):
        entrystring = self.getvalue().strip()
        try:
            self.content = self.read(entrystring)
            self.modified = 0
        except Exception as e:
            self.filtus.notbusy()
            typ = type(e).__name__
            FiltusUtils.warningMessage("There is a problem with the entry '%s':\n\n%s: %s" %(entrystring, typ, e))
            self.deselect()
            self.content = set()
            
    def getcontent(self):
        if self.modified: self.setcontent()
        return self.content

    def getvalue(self):
        return self.entryfield.getvalue()

    def setvalue(self, value):
        self.entryfield.setvalue(value)
        self.entryfield.component('entry').xview('end')

    def _togglestatus(self):
        self.entryfield.configure(entry_state='normal' if self.on() else 'disabled')

    def _browseCommand(self):
        fil = tkFileDialog.askopenfilename(initialdir=self.filtus.currentDir, title = self['browsetitle'])#, filetypes = self['filetypes'])
        if fil:
            self.filtus.currentDir = os.path.dirname(fil)
            self.setvalue(os.path.normpath(fil))
            self.select()
            self.entryfield.invoke()
            

class GeneFile(FileBrowser):
    def __init__(self, parent, filtus, **kw):
        FileBrowser.__init__(self, parent, filtus, **kw)

    @staticmethod
    def read(string):
        if ';' in string:
            return set.union(*[GeneFile.read(part.strip()) for part in string.split(';')])
        genes = set()
        if os.path.isfile(string):
            with open(string, 'rU') as genfil:
                genes = set(g.strip() for g in genfil)
        elif len(string)>0:
            genes = set(gen.strip() for gen in string.split(','))
        return genes
     
        
class RegionFile(FileBrowser):
    def __init__(self, parent, filtus, **kw):
        FileBrowser.__init__(self, parent, filtus, **kw)

    @staticmethod
    def read(string):

        def _unpack(reg):
            try:
                chr, startstop = reg.split(':')
                start, stop = startstop.split('-')
                return (chr, float(start), float(stop))
            except Exception as e:
                raise ValueError("Illegally specified region: %s"%str(reg))
        
        if ' AND ' in string:
            return FiltusUtils.region_intersection(*[RegionFile.read(part.strip()) for part in string.split(' AND ')])
        res = []
        if os.path.isfile(string):
            with open(string, 'rU') as regions:
                regs = [line.split() for line in regions if not line.startswith('##')]
            if all(x in regs[0] for x in ('CHR', 'POS1', 'POS2')): # reading plink .hom files!
                get3 = itemgetter(regs[0].index('CHR'), regs[0].index('POS1'), regs[0].index('POS2'))
                regs = [get3(r) for r in regs]
            else:
                regs = [r[:3] for r in regs]
            
            if all(len(r) == 1 for r in regs):
                res = [_unpack(r[0]) for r in regs]
            else:
                if 'chr' in regs[0][0].lower(): regs = regs[1:]
                res = [(chr, float(start), float(stop)) for chr, start, stop in regs]
        else:
            try:
                res = [_unpack(r.strip()) for r in string.split(';')]
            except ValueError as e:
                raise ValueError('%s\n\nRegions should be written as "chr:start-stop" and be separated by semicolon.\n\nExample:  "3:1-10000; 12:5e7-1e8"' %e)
        return res


class VariantFile(FileBrowser):
    def __init__(self, parent, filtus, **kw):
        kw['browsetitle'] = "Select file with exclusion variants"
        kw['state'] = 'disabled'
        FileBrowser.__init__(self, parent, filtus, **kw)

    def read(self, string):
        return self.read2(string, self.filtus)

    @staticmethod
    def read2(string, filtus):
        if ';' in string:
            res = set.union(*[VariantFile.read2(part.strip(), filtus=filtus) for part in string.split(';')])
            filtus.storage[string] = res
            return res
        if not os.path.isfile(string):
            raise IOError("File does not exist: %s."%string)
        res = set()
        dbReader = InputDialog.InputDialog(filtus=filtus)

        VF = dbReader.read(string, prompt=0, guess=1, sep='\t',  gtCol=None)
        if VF:
            res = VF.getUniqueVariants()
            filtus.storage[string] = res
        return res
    
    
class GeneLengthFile(FileBrowser):
    def __init__(self, parent, filtus, **kw):
        FileBrowser.__init__(self, parent, filtus, **kw)

    @staticmethod
    def read(file):
        file = os.path.normpath(file)
        res = {}
        if not os.path.isfile(file):
            raise IOError("Gene length file does not exist: %s."% file)
        try:
            with open(file, 'rU') as f:
                for line in f:
                    try:
                        gen, length = line.split('\t')[:2]
                        res[gen] = int(length)
                    except ValueError:
                        continue
            if not res:
                raise IOError("Empty file: %s."%file)
        except ValueError as e:
            raise ValueError("Gene length file does not have correct format (i.e. two tab-separated columns).\n\n%s" % e)
        return res
  
class GeneLookup(object):
    def __init__(self, parent, filtus):
        #self.genes = []
        self.filtus = filtus
        self.geneViewer = GeneView(filtus, title="Gene view of variants in all samples")
        self.prompt = Pmw.PromptDialog(parent, buttons=('OK', 'Cancel'), title='Gene lookup',
                            label_text='Gene name(s):', entryfield_labelpos='n', command=self.execute, defaultbutton=0)
                            #entryfield_value=', '.join(self.genes), 
        self.prompt.withdraw()

    def execute(self, button):
        self.prompt.deactivate()
        if button is None or button == 'Cancel': return

        genes = [g.strip() for g in self.prompt.get().split(',')]
        variants = FiltusAnalysis.geneLookup(genes, VFlist=self.filtus.filteredFiles)
        self.geneViewer.display(variants)

        
class ModelSelect(Pmw.RadioSelect):
    def __init__(self, parent, **kw):
        kw['buttontype'] = 'radiobutton'
        kw['orient'] = 'horizontal'
        kw['labelpos'] = 'w'
        kw['label_text'] = 'Model:'
        apply(Pmw.RadioSelect.__init__, (self, parent), kw)
        self.add('Dominant', text='Dominant')
        self.add('Recessive', text='Recessive c/h')
        self.add('Recessive homozygous', text='Recessive homoz')
        self.setvalue('Dominant')

    def getButtons(self):
        return [self.button(i) for i in range(3)]

        
def makeReadOnly(tkWdg):
    """Makes a Tk widget (typically an Entry or Text) read-only, in the
    sense that the user cannot modify the text (but it can still be set
    programmatically). The user can still select and copy text and key
    bindings for <<Copy>> and <<Select-All>> still work properly.
    Inputs:
    - tkWdg: a Tk widget
    """
    def killEvent(evt): return "break"
    def doCopy(evt):    tkWdg.event_generate("<<Copy>>")
    def selectAll(evt):
        try:
            tkWdg.tag_add("sel", '1.0', "end")
        except Exception as e:
            pass
        finally:
            return "break"


    tkWdg.bind("<<Cut>>", killEvent)
    tkWdg.bind("<<Paste>>", killEvent)
    tkWdg.bind("<<Paste-Selection>>", killEvent)
    tkWdg.bind("<<Clear>>", killEvent)
    tkWdg.bind("<Button-2>", killEvent)
    tkWdg.bind("<ButtonRelease-2>", killEvent)
    tkWdg.bind("<Key>", killEvent)

    tkWdg.bind('<Control-a>', selectAll)

    # restore copy and select all
    for evt in tkWdg.event_info("<<Copy>>"):  tkWdg.bind(evt, doCopy)


class AutEx_GUI(Pmw.Dialog):
    def __init__(self, filtus):
        
        Pmw.Dialog.__init__(self, filtus.parent, buttons = ('Compute','Close'), title = 'AutEx', defaultbutton=0,
            command=self.execute, activatecommand=self._prepare, buttonbox_pady=10, dialogchildsite_pady=0)
        self.withdraw()
        interior0 = self.interior()
        interior = Tkinter.Frame(interior0) #self.interior()
        interior.grid(row=0, column=0, pady=10, sticky='news')
        HelpButton(interior0, filtus=filtus, page="autex").grid(row=0, column=0, sticky="ne")
        
        self.filtus = filtus
        self.genmapfile = os.path.join(filtus.datadir, "DecodeMap_thin.txt")
        if not os.path.isfile(self.genmapfile): 
            print "File 'DecodeMap_thin.txt' not found"
            self.genmapfile = None
        self.autexComputer = FiltusAnalysis.AutExComputer(genmapfile = self.genmapfile)
        self.VF = None
        
        # Title label
        Tkinter.Label(interior, text='Autozygosity mapping', font=filtus.titlefont).grid(padx=20, pady=10)
        
        entry_opt = dict(entry_width=6, labelpos='w', labelmargin=5, entry_justify="center")
        grid_opt = dict(padx=10, pady=5, sticky='news')
        
        # Parameter input
        param_group = Pmw.Group(interior, tag_text = 'Model parameters')
        param_interior = param_group.interior()
        
        self.inbreedingOM = Pmw.OptionMenu(param_interior, labelpos = 'w', label_text = 'Parental relation:', items = ['Distant', '1st - 2nd cousins', 'Close'],
                initialitem='Distant', labelmargin=3, menubutton_pady=2, menubutton_anchor = 'w',)
        
        self.inbreedingOM.grid(**dict(grid_opt, sticky='w'))
        param_group.grid(ipady=2, **grid_opt)
        
        # Frequency input
        freq_group = Pmw.Group(interior, tag_text = 'Allele frequencies')
        freq_interior = freq_group.interior()
        
        self._altFreqMenu = OptionMenuExt(freq_interior, label_text = "ALT frequency column:", labelpos='nw', 
                                        menu_font=filtus.defaultfont, menubutton_anchor = 'w', labelmargin=2,
                                        menubutton_padx=5, menubutton_pady=1, menubutton_width=18)
        self._def_freq_entry = Pmw.EntryField(freq_interior, label_text = "Missing entry value:", entry_justify="center",
                                        value='0.1', entry_width=6, labelpos='nw', labelmargin=4)
        
        self._altFreqMenu.grid(row=0, column=0, sticky='news', padx=(10,25), pady=(5,10))
        self._def_freq_entry.grid(row=0, column=1, sticky='nes', padx=10, pady=(5,11))
        
        freq_group.grid(**grid_opt)
        
        # Output options
        out_group = Pmw.Group(interior, tag_text = 'Output')
        out_interior = out_group.interior()
        
        self._thresh_entry = Pmw.EntryField(out_interior, label_text = "Posterior threshold:", value='0.5', **entry_opt)
        self._minlength_entry = Pmw.EntryField(out_interior, label_text = "Minimum segment size:", value='0',**entry_opt)
        self._unitmenu = Pmw.OptionMenu(out_interior, items=['Mb*','Mb','cM*','cM'], initialitem=2, menubutton_anchor = 'e', menubutton_width = 4, menubutton_padx=0, menubutton_pady=1)
        self._mincount_entry = Pmw.EntryField(out_interior, label_text = " and", value='0', **entry_opt)
        
        self._thresh_entry.grid(padx=10, pady=5, sticky='w')
        self._minlength_entry.grid(row=1, column=0, padx=(10,0), pady=5, sticky='w')
        self._unitmenu.grid(row=1, column=1, padx=0, pady=5, sticky='w')
        self._mincount_entry.grid(row=1, column=2, padx=0, pady=5, sticky='w')
        Tkinter.Label(out_interior, text="variants ").grid(row=1, column=3, padx=0, pady=5, sticky='w')
        out_group.grid(ipady=2, **grid_opt)
        
        
        below = Tkinter.Frame(interior)
        below.columnconfigure(0, weight=1)
        below.columnconfigure(1, weight=1)
        summary_group = Pmw.Group(below, tag_text = 'Summary')
        summary_interior = summary_group.interior()
        self.summary = Pmw.ScrolledText(summary_interior, columnheader=1, columnheader_padx=2, columnheader_width = 42, text_wrap='none', text_padx=3, text_pady=3,
                        text_height=5, text_font=filtus.monofont, #labelpos='nw', label_text='Summary', 
                        text_width=42, borderframe=0, scrollmargin=0, vscrollmode='none', hscrollmode='none')
        self.summary.grid(**grid_opt)
        summary_group.grid(**grid_opt)
        
        plot_group = Pmw.Group(below, tag_text = 'Plot')
        plot_interior = plot_group.interior()
        self.plotCHR = OptionMenuExt(plot_interior, label_text = "Chrom:", menubutton_padx=1, menubutton_anchor="e", menubutton_pady=1, menubutton_width=3, 
                        labelmargin = 1, labelpos='nw')
        self.plotButton = Tkinter.Button(plot_interior, text="Plot now!", command=self.plot, width=5)
        if not PLOT_available: 
            self.plotCHR.configure(menubutton_text = " None")
            self.plotCHR.disable()
            self.plotButton.configure(state = "disabled")
        
        self.plotCHR.grid(row=0, column=0, sticky='ew', padx=12, pady=(0,5))
        self.plotButton.grid(row=1, sticky='ew', padx=12, pady=10)
        plot_group.grid(**dict(grid_opt, row=0, column=1))  
        
        starlabel = Tkinter.Label(below, text="* = Extending segments to nearest outside variants", font=filtus.smallfont)
        starlabel.grid(sticky='nw', padx=10, pady=0, columnspan=2)
        below.grid(sticky='news')
        
        
        Pmw.alignlabels([self.inbreedingOM, self._thresh_entry, self._minlength_entry]) #self.plotCHR, 
        
        for g in [param_group, freq_group, out_group, summary_group, plot_group]:
            g.interior().columnconfigure(0, weight=1)
            g.configure(tag_font = filtus.smallbold)
            
            
    def plot(self):
        filtus = self.filtus
        filtus.busy()
        autex = self.autexComputer
        try:
            VF = self.VF
            chrom = self.plotCHR.getvalue()
            if not chrom: raise RuntimeError("No chromosome indicated")
            params = self._readInput()
            pos_phys, obs, freqs, pos_cm, scores = autex.singleChrom_scores(VF, chrom=chrom, f=params['f'], a=params['a'], error=params['error'], 
                                                                            defaultFreq=params['defaultFreq'], altFreqCol=params['altFreqCol'])
            segs = autex.scores2segments(scores, pos_phys=pos_phys, pos_cM=pos_cm, threshold=params['threshold'], minlength=params['minlength'],
                                          mincount=params['mincount'], unit=params['unit'])
            title = "%s - chromosome %s" %(VF.shortName, chrom)
            FiltusQC.homPlot(pos_phys, obs, scores, freqs, title=title, segs=segs)
            filtus.notbusy()
        except Exception as e:
            filtus.notbusy()
            FiltusUtils.warningMessage("%s: %s"% (type(e).__name__, e))
            return
   
    def execute(self, button):
        if button is None or button=='Close':
            return self.deactivate()
        try:
            filtus = self.filtus
            filtus.busy()
            params = self._readInput()
            segments = self.autexComputer.autex_segments(self.VF, **params)
            filtus.text.prettyPrint(segments, label="Autozygous segments of %s"%self.VF.shortName)
            
            summ = self.autexComputer.summary(segments)
            summ = {k:round(v, 4) for k,v in summ.iteritems()}
            summ_header = '%-11s%-8s%-8s%-8s%-8s' %('', 'Mb*', 'Mb', 'cM*', 'cM')
            summary_txt = '\n'.join('%-11s%-8.3g%-8.3g%-8.3g%-8.3g' %(r+':', summ[r+' MB*'], summ[r+' MB'], summ[r+' CM*'], summ[r+' CM']) for r in ['Total', 'Longest'])
            summary_txt += '\n%-11s%-8s%-8s%-8.3g%-8.3g\n\n' %('Fraction:', '','', summ['Fraction*'], summ['Fraction'])
            summary_txt += '%d segments in total' % summ['Number of segments']
            self.summary.component('columnheader').delete('1.0', 'end')
            self.summary.component('columnheader').insert('1.0', summ_header)
            self.summary.settext(summary_txt)
            filtus.notbusy()
        except Exception as e:
            filtus.notbusy()
            FiltusUtils.warningMessage("%s: %s"% (type(e).__name__, e))
            return

    #### Private methods ###
    
    def _prepare(self):
        self.VF = VF = self.filtus.checkLoadedSamples(select="selection", VF=True, minimum=1, maximum=1)[0]
        self._altFreqMenu.setItems([''] + VF.columnNames)
        chr = VF.chromGetter
        chroms = [''] + sorted(set(chr(v) for v in VF.variants), key=FiltusUtils.chromInt)
        chroms = [x for x in chroms if not any(XYstring in x for XYstring in ("X","Y","x","y"))]
        self.plotCHR.setItems(chroms)
        
    def _readInput(self):
        inb = self.inbreedingOM.index(self.inbreedingOM.getvalue())
        if inb == 0: f,a = 0.01, 0.5  # mean IBD = 2 cM    
        elif inb == 1: f,a = 0.06, 0.06 # mean IBD = 17 cM
        elif inb == 2: f,a = 0.13, 0.04 # mean IBD = 28 cM
        
        error = 0.005 #float(self.err_entry.getvalue())
        threshold = float(self._thresh_entry.getvalue())
        minlength = float(self._minlength_entry.getvalue())
        unit = self._unitmenu.getvalue()
        mincount = int(self._mincount_entry.getvalue())
        altFreqCol = self._altFreqMenu.getvalue()
        try:
            def_freq = float(self._def_freq_entry.getvalue())
            if not 0 < def_freq < 1: raise ValueError
        except ValueError:
            raise ValueError("Please indicate a 'Missing entry value' strictly between 0 and 1. (This is used as the default allele frequency).")
        
        return dict(f=f, a=a, error=error, altFreqCol=altFreqCol, defaultFreq=def_freq, 
                    threshold=threshold, minlength=minlength, unit=unit, mincount=mincount)

                  
class PLINK_GUI(Pmw.Dialog):
    def __init__(self, filtus):
        Pmw.Dialog.__init__(self, filtus.parent, buttons = ('Compute','Close'), title = 'Plink homozygosity', defaultbutton=0,
            command=self.execute, activatecommand=self._prepare, buttonbox_pady=10, dialogchildsite_pady=10)
        self.withdraw()
        self.filtus = filtus
        self.plinkComputer = FiltusAnalysis.PlinkComputer()
        
        self.VF = None
        self.plinkDefaults = {'homozyg-window-kb' : 5000, 'homozyg-window-snp' : 50, 'homozyg-window-het' : 1,
                                'homozyg-window-missing' : 5, 'homozyg-window-threshold' : 0.05, 'homozyg-snp' : 100,
                                'homozyg-kb' : 1000, 'homozyg-density' : 50, 'homozyg-gap' : 1000}
        self.strictDefaults = {'homozyg-window-kb' : 0, 'homozyg-window-snp' : 1, 'homozyg-window-het' : 0,
                                'homozyg-window-missing' : 0, 'homozyg-window-threshold' : 1, 'homozyg-snp' : 10,
                                'homozyg-kb' : 0, 'homozyg-density' : 10000, 'homozyg-gap' : 1000000}
        self.paramList = ['homozyg-window-kb','homozyg-window-snp', 'homozyg-window-het','homozyg-window-missing',
                                'homozyg-window-threshold', 'homozyg-snp', 'homozyg-kb','homozyg-density', 'homozyg-gap']
        interior = self.interior()
        
        # Title label
        Tkinter.Label(interior, text='Homozygosity mapping (PLINK)', font=filtus.titlefont).grid(padx=20, pady=10)
        
        entry_opt = dict(entry_width=6, labelpos='w', labelmargin=5, entry_justify="center")
        grid_opt = dict(padx=10, pady=5, sticky='news')
        
        # Parameter input
        param_group = Pmw.Group(interior, tag_text = 'Model parameters')
        param_interior = param_group.interior()
        
        self.entryFields = []
        for param in self.paramList:
            entr = Pmw.EntryField(param_interior, label_text = param+':', labelpos='w', entry_width=9,  entry_state='normal', 
                                  validate = 'real' if 'threshold' in param else 'numeric', value = self.plinkDefaults[param])
            entr.grid(padx=20, pady=5)
            self.entryFields.append(entr)
        
        reset = Pmw.RadioSelect(param_interior, label_text = "Quick select:", labelpos='w', buttontype = 'radiobutton',
                                        orient = 'horizontal', command=self.resetValues)
        reset.add('PLINK defaults')
        reset.add('Strict ROHs')
        reset.invoke('Strict ROHs')
        reset.grid(pady=5)
        self.align()
        param_group.grid(ipady=2, **grid_opt)
        
        out = Tkinter.Frame(interior)
        out.columnconfigure(0, weight=1)
        out.columnconfigure(1, weight=1)
        
        summary_group = Pmw.Group(out, tag_text = 'Summary')
        summary_interior = summary_group.interior()
        self.summary = Pmw.ScrolledText(summary_interior, columnheader=0, #columnheader_padx=2, columnheader_width = 42,
                        text_height=5, text_font=filtus.monofont, text_wrap='none', text_padx=3, text_pady=3,
                        text_width=42, borderframe=0, scrollmargin=0, vscrollmode='none', hscrollmode='none')
        self.summary.grid(**grid_opt)
        summary_group.grid(**grid_opt)
        
        plot_group = Pmw.Group(out, tag_text = 'Plot')
        plot_interior = plot_group.interior()
        self.plotCHR = OptionMenuExt(plot_interior, label_text = "Chrom:", menubutton_padx=1, menubutton_anchor="e", menubutton_pady=1, menubutton_width=3, 
                        labelmargin = 1, labelpos='nw')
        self.plotButton = Tkinter.Button(plot_interior, text="Plot now!", command=self.plot, width=5)
        if not PLOT_available: 
            self.plotCHR.configure(menubutton_text = " None")
            self.plotCHR.disable()
            self.plotButton.configure(state = "disabled")
        
        self.plotCHR.grid(row=0, column=0, sticky='ew', padx=12, pady=(0,5))
        self.plotButton.grid(row=1, sticky='ew', padx=12, pady=10)
        plot_group.grid(**dict(grid_opt, row=0, column=1))  
        
        out.grid(sticky='news')
        
        for g in [param_group, summary_group, plot_group]:
            g.interior().columnconfigure(0, weight=1)
            g.configure(tag_font = filtus.smallbold)
            
    def align(self):
        Pmw.alignlabels(self.entryFields)

    def resetValues(self, button):
        newvals = self.plinkDefaults if button == 'PLINK defaults' else self.strictDefaults
        for param, entry in zip(self.paramList, self.entryFields):
                entry.setvalue(newvals[param])

            
    def plot(self):
        filtus = self.filtus
        filtus.busy()
        try:
            inputchrom = self.plotCHR.getvalue()
            if not inputchrom: raise RuntimeError("No chromosome indicated")
            VF = Filter.Filter(columnfilters=[(self.VF.chromCol, 'equal to', inputchrom, 0)]).apply(self.VF)
            
            cleanup, dir = True, filtus.currentDir
            pedprefix = outprefix = 'filtus__plink__homozyg'
            segments = self.plinkComputer.runPlink(VF, pedprefix=pedprefix, outprefix=outprefix, dir=dir, verbose=False, **self._readInput())
            chr, pos, gtnum = VF.chromGetter, VF.posGetter, VF.GTnum()
            chromDataSorted = sorted((int(pos(v)), gtnum(v)) for v in VF.variants)
            pos_phys, obs = zip(*chromDataSorted)
            segs = [(float(s[2]), float(s[3]), float(s[4])/1000) for s in segments.variants]
            title = "PLINK homozygosity: %s - chromosome %s" %(VF.shortName, inputchrom)
            FiltusQC.homPlotSimple(pos_phys, obs, title=title, segs=segs)
            filtus.notbusy()
        except Exception as e:
            filtus.notbusy()
            FiltusUtils.warningMessage("%s: %s"% (type(e).__name__, e))
            return
   
    def execute(self, button):
        if button is None or button=='Close':
            return self.deactivate()
        try:
            filtus = self.filtus
            filtus.busy()
            params = self._readInput()
            
            cleanup = True
            pedprefix = outprefix = 'filtus__plink__homozyg'
            dir = filtus.currentDir
            segments = self.plinkComputer.runPlink(self.VF, pedprefix=pedprefix, outprefix=outprefix, dir=dir, verbose=False, **params)
            filtus.text.prettyPrint(segments, label="PLINK homozygous segments for %s" %self.VF.shortName)
    
            summ = self.plinkComputer.summary(segments)
            summ = {k:round(v, 4) for k,v in summ.iteritems()}
            summary_txt = '\n'.join('%-15s%.3g' %(k+':', summ[k]) for k in ['Total (MB)', 'Longest (MB)', 'Segments', 'Fraction'])
            self.summary.settext(summary_txt)
            
            filtus.notbusy()
        except Exception as e:
            filtus.notbusy()
            FiltusUtils.warningMessage("%s: %s"% (type(e).__name__, e))
            return
        finally:
            if cleanup: self._cleanup(dir, pedprefix, outprefix)
                
                
    #### Private methods ###
    def _cleanup(self, dir, pedprefix, outprefix):
        if not dir: dir = os.getcwd()
        pedfiles = [pedprefix + '.ped', pedprefix + '.map']
        for f in os.listdir(dir):
            if f.startswith(outprefix) or f in pedfiles:
                os.remove(f)
    
    def _prepare(self):
        self.VF = VF = self.filtus.checkLoadedSamples(select="selection", VF=True, minimum=1, maximum=1)[0]
        chr = VF.chromGetter
        chroms = [''] + sorted(set(chr(v) for v in VF.variants), key=FiltusUtils.chromInt)
        chroms = [x for x in chroms if not any(XYstring in x for XYstring in ("X","Y","x","y"))]
        self.plotCHR.setItems(chroms)
        
    def _readInput(self):
        return {param:entry.getvalue().strip() for param, entry in zip(self.paramList, self.entryFields)}
        

class DeNovo_GUI(Pmw.Dialog):
    def __init__(self, filtus):
        
        Pmw.Dialog.__init__(self, filtus.parent, buttons = ('Compute','Close'), title = 'DeNovo', defaultbutton=0,
            command=self.execute, activatecommand=self._prepare, buttonbox_pady=10, dialogchildsite_pady=0)
        self.withdraw()
        interior0 = self.interior()
        interior = Tkinter.Frame(interior0) #self.dialog.interior()
        interior.columnconfigure(0, weight=1)
        interior.grid(row=0, column=0, pady=10, sticky='news')
        HelpButton(interior0, filtus=filtus, page="denovo").grid(row=0, column=0, sticky="ne")
        
        self.filtus = filtus
        self.denovoComputer = FiltusAnalysis.DeNovoComputer()
        
        # Title label
        Tkinter.Label(interior, text='Detection of de novo variants', font=filtus.titlefont).grid(padx=20, pady=10)
        
        entry_opt = dict(entry_width=6, labelpos='w', labelmargin=5, label_anchor='w')
        grid_opt = dict(padx=10, pady=5, sticky='news')
        
        # Parameter input
        trio_group = Pmw.Group(interior, tag_text = 'Trio samples')
        trio_interior = trio_group.interior()
        trio_interior.columnconfigure(0, weight=1)
        trio_interior.columnconfigure(1, weight=1)
        trio_interior.columnconfigure(2, weight=1)
            
        self.child = Pmw.EntryField(trio_interior, label_text = "Child:", entry_justify="left", **entry_opt)#, modifiedcommand=self._setFreqCols)
        self.father = Pmw.EntryField(trio_interior, label_text = "Father:", entry_justify="left", **entry_opt)
        self.mother = Pmw.EntryField(trio_interior, label_text = "Mother:", entry_justify="left", **entry_opt)
        self.boygirl = Pmw.RadioSelect(trio_interior, label_text = "Child gender:", labelpos='w', buttontype = 'radiobutton', orient = 'horizontal')
        self.boygirl.add("Boy")
        self.boygirl.add("Girl")
        self.boygirl.invoke("Girl")
        
        self.child.grid(row=0, column=0, **grid_opt)
        self.father.grid(row=0, column=1, **grid_opt)
        self.mother.grid(row=0, column=2, **grid_opt)
        self.boygirl.grid(row=1, column=0, columnspan=3, **grid_opt)
        
        mut_group = Pmw.Group(interior, tag_text = 'Mutation rate')
        mut_interior = mut_group.interior()
        
        self.mutEntry = Pmw.EntryField(mut_interior, label_text = "Prior mutation rate:", value="1e-8", entry_justify="center", **entry_opt)
        self.mutEntry.grid(padx=10, pady=5, sticky='nw')
        
        # Frequency input
        freq_group = Pmw.Group(interior, tag_text = 'Allele frequencies')
        freq_interior = freq_group.interior()
        self._altFreqMenu = OptionMenuExt(freq_interior, label_text = "ALT frequency column:", labelpos='nw', 
                                        menu_font=filtus.defaultfont, menubutton_anchor = 'w', labelmargin=2,
                                        menubutton_padx=5, menubutton_pady=1, menubutton_width=18)
        self._def_freq_entry = Pmw.EntryField(freq_interior, label_text = "Missing entry value:", entry_justify="center",
                                        value='0.1', entry_width=6, labelpos='nw', labelmargin=4)
        
        self._altFreqMenu.grid(row=1, column=0, sticky='news', padx=(10,25), pady=(5,10))
        self._def_freq_entry.grid(row=1, column=1, sticky='nes', padx=10, pady=(5,11))
        
        # Output options
        out_group = Pmw.Group(interior, tag_text = 'Output filters')
        out_interior = out_group.interior()
        
        ALTframe = Tkinter.Frame(out_interior)
        ALTframe.columnconfigure(0, weight=1)
        
        Tkinter.Label(ALTframe, text="Percentage of reads with ALT allele:").grid(rowspan=2, sticky='nsw', padx=(10,5), pady=5)
        self._minALTchild_entry = Pmw.EntryField(ALTframe, label_text = "Child >", value='', **entry_opt)
        self._maxALTparent_entry = Pmw.EntryField(ALTframe, label_text = "Parents <", value='', **entry_opt)
        
        self._minALTchild_entry.grid(row=0, column=1, padx=(0,10), pady=2, sticky='w')
        self._maxALTparent_entry.grid(row=1, column=1, padx=(0,10), pady=2, sticky='w')
        Pmw.alignlabels([self._minALTchild_entry, self._maxALTparent_entry], 'e') 
        ALTframe.grid(columnspan=2, sticky='news')
        
        self._thresh_entry = Pmw.EntryField(out_interior, label_text = "Posterior probability >", value='', **entry_opt)#labelpos='nw', labelmargin=4, entry_width=6)
        self._thresh_entry.grid(padx=10, pady=5, sticky='w')
        
        summary_group = Pmw.Group(interior, tag_text = 'Summary')
        summary_interior = summary_group.interior()
        self.summary = Pmw.ScrolledText(summary_interior, columnheader=0, text_height=2, text_font=filtus.monofont, text_wrap='none', text_padx=3, text_pady=3,
                        text_width=50, borderframe=0, scrollmargin=0, vscrollmode='none', hscrollmode='none')
        makeReadOnly(self.summary.component('text'))
        self.summary.grid(**grid_opt)
        
        #Pmw.alignlabels([self._thresh_entry, self.mutEntry, self._minALTchild_entry, self._maxALTparent_entry]) 
        Pmw.alignlabels([self.mutEntry, self._thresh_entry]) 
        
        for g in [trio_group, mut_group, freq_group, out_group, summary_group]:
            g.interior().columnconfigure(0, weight=1)
            g.configure(tag_font = filtus.smallbold)
            g.grid(ipady=2, **grid_opt)
        
    def execute(self, button):
        if button is None or button=='Close':
            return self.deactivate()
        try:
            filtus = self.filtus
            filtus.busy()
            input = self._readInput()
            result = self.denovoComputer.analyze(**input) 
            self.summary.settext("%d variants found!\nVariants are shown in the main window." % result.length)
            filtus.text.prettyPrint(result, rightClick="variantMenu", label="De novo variants in %s"%input['VFch'].shortName)
            filtus.notbusy()
        except Exception as e:
            filtus.notbusy()
            FiltusUtils.warningMessage("%s: %s"% (type(e).__name__, e))
            return

    #### Private methods ###
    
    def _prepare(self):
        VF = self.filtus.checkLoadedSamples(select="all", VF=True, minimum=3)[0]
        self._altFreqMenu.setItems([''] + VF.columnNames)
        
    def _readInput(self):
        loadedVFs = self.filtus.filteredFiles
        loadedFilenames = [VF.longName for VF in loadedVFs]
        ch_fa_mo = []
        for name, entryfield in zip(['Child', 'Father', 'Mother'], [self.child, self.father, self.mother]):
            try:
                entry = entryfield.getvalue().strip()
                indices = FiltusUtils.convert2indices(entry, idlist=loadedFilenames)
            except Exception as e:
                raise ValueError("Invalid %s: %s\n\nDetails:\n%s"%(name.lower(), entry, e))
            if not indices: # indices is a list, so [0] will pass this.
                raise ValueError("Invalid %s: Field cannot be empty." % name.lower())
            if len(indices) > 1:
                raise ValueError("Invalid %s: Please enter the ID of exactly one individual." % name.lower())    
            ch_fa_mo.append(indices)
            
        FiltusUtils.checkDuplicates(ch_fa_mo)
        ch_fa_mo_flatten = [a[0] for a in ch_fa_mo]
        VFch, VFfa, VFmo = [loadedVFs[i] for i in ch_fa_mo_flatten]
        
        boygirl = self.boygirl.getvalue()
        threshold = float(self._thresh_entry.getvalue()) if self._thresh_entry.getvalue() else None 
        minALTchild = float(self._minALTchild_entry.getvalue()) if self._minALTchild_entry.getvalue() else None 
        maxALTparent = float(self._maxALTparent_entry.getvalue()) if self._maxALTparent_entry.getvalue() else None 
        
        altFreqCol = self._altFreqMenu.getvalue()
        def_freq = self._def_freq_entry.getvalue()
        if not def_freq or not 0 < float(def_freq) < 1: 
            raise ValueError("Please indicate a 'Missing entry value' strictly between 0 and 1. (This is used as the default allele frequency).")
        def_freq = float(def_freq)
        mut = self.mutEntry.getvalue()
        if not mut or not 0 < float(mut) < 1: 
            raise ValueError("Please indicate a prior mutation rate strictly between 0 and 1.")
        mut = float(mut)
        
        return dict(VFch=VFch, VFfa=VFfa, VFmo=VFmo, boygirl=boygirl, trioID=ch_fa_mo_flatten, mut=mut, altFreqCol=altFreqCol, defaultFreq=def_freq, 
                    threshold=threshold, minALTchild=minALTchild, maxALTparent=maxALTparent)
