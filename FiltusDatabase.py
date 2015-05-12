import os.path
import time
import tkFileDialog
import collections
import itertools
from operator import itemgetter

import Tkinter
import Pmw
import FiltusUtils
import FiltusWidgets
import Filter
import DataContainer

def _readTop(file):
    '''
    Input: An open database file object
    Output: List of meta lines; list of column names; remaining file object
    '''
    meta = []
    line = file.next()
    while line.startswith('##'):
        meta.append(line.strip())
        line = file.next()
    colNames = line.strip().split('\t')
    return meta, colNames, file

    
class TwoListWidget(Pmw.MegaWidget):
    def __init__(self, parent, filtus, **kw):
        # Define the megawidget options.
        optiondefs = (
            ('height', 8, Pmw.INITOPT),
            ('width', 32, Pmw.INITOPT),
            ('items', (), None),
            ('lefttoptext', '', None),
        )
        self.defineoptions(kw, optiondefs)

        # Initialise base class (after defining options).
        Pmw.MegaWidget.__init__(self, parent)
        self.filtus = filtus
        self.selection = []
        
        #create components
        interior = self.interior()

        self._editdialog = Pmw.PromptDialog(filtus.parent, title = 'Edit', label_text = 'Edit sample name:', entryfield_labelpos = 'nw',
            entryfield_entry_width=50, defaultbutton = 0, buttons = ('OK', 'Cancel'), activatecommand=self.insertselected, command=self.doEdit)
        self._editdialog.withdraw()
        
        self._leftlist = self.createcomponent('leftlist', (), None, FiltusWidgets.LabeledListBox, (interior,), filtus=filtus,
                toptext=self['lefttoptext'], width=self['width'], height=self['height'], 
                scrolledlistbox_selectioncommand=self.setButtonStates, scrolledlistbox_items=self['items'], 
                scrolledlistbox_dblclickcommand=self.select)
        
        self._leftlist.bind("<Return>", self.select)
        
        self._rightlist = self.createcomponent('rightlist', (), None, FiltusWidgets.LabeledListBox, (interior,), filtus=filtus,
                                toptext="Selected: 0", bottomtext="Double click entry to edit names", width=self['width'], 
                                height=self['height'], scrolledlistbox_selectioncommand=self.setButtonStates,
                                scrolledlistbox_dblclickcommand=self._editdialog.activate)
                                
        self._buttonColumn =  self.createcomponent('buttoncolumn', (), None, Pmw.ButtonBox, (interior,), orient="vertical")
        
        self._buttonColumn.add("selectButton", text="Select", command=self.select)
        self._buttonColumn.add("deselectButton", text="Deselect", command=self.deselect)
        self._buttonColumn.add("selectAllButton", text="Select all", command=self.selectAll)
        self.initialiseoptions()

        interior.columnconfigure(0, weight=1)
        interior.columnconfigure(2, weight=1)
        interior.rowconfigure(0, weight=1)
        
        self._leftlist.grid(sticky='news', row=0, column=0)
        self._buttonColumn.grid(sticky='news', row=0, column=1, pady=5, padx=30)
        self._rightlist.grid(sticky='news', row=0, column=2)
    
    
    def setItems(self, items):
        self._leftlist.setlist(items)
        
    def setButtonStates(self, state=None):
        selButt = self._buttonColumn.button("selectButton")
        deselButt = self._buttonColumn.button("deselectButton")
        allButt = self._buttonColumn.button("selectAllButton")
        if state is not None: 
            # if called explicitly - set all buttons
            for b in (selButt, deselButt, allButt): 
                b.configure(state=state)
                self._leftlist.component('scrolledlistbox_label').configure(state=state)
                self._rightlist.component('scrolledlistbox_label').configure(state=state)
                self._rightlist.component('bottomlabel').configure(state=state)
        else: 
            # callback from list widgets
            if self._leftlist.curselection():
                selButt.configure(state="normal")
                deselButt.configure(state="disabled")
            elif self._rightlist.curselection():
                selButt.configure(state="disabled")
                deselButt.configure(state="normal")
            
    def select(self):
        box1, box2 = self._leftlist, self._rightlist
        filelist = box1.get()
        sel = box1.getcurselection()
        sel_ind = [filelist.index(s) for s in sel]
        taken = [i for i in sel_ind if i in self.selection]
        
        if len(taken)>0: 
            FiltusUtils.warningMessage("Samples already selected:\n\n%s" %'\n'.join(filelist[i] for i in taken))
            sel = [s for s,i in zip(sel, sel_ind) if not i in taken]
            sel_ind = [i for i in sel_ind if not i in taken]
            
        box2.insert('end', *sel)
        self.selection.extend(sel_ind)
        box2.settoptext("Selected: %d" % box2.size())
        
    def selectAll(self):
        box1, box2 = self._leftlist, self._rightlist
        box2.setlist(box1.get())
        box2.settoptext("Selected: %d" % box2.size())
        self.selection = range(box2.size())
        
    def deselectAll(self):
        self._rightlist.setlist([])
        self._rightlist.settoptext("Selected: 0")
        self.selection = []
        
    def deselect(self):
        box2 = self._rightlist
        for i in reversed([int(i) for i in box2.curselection()]):
             box2.delete(i)
             del self.selection[i]
        box2.settoptext("Selected: %d" % box2.size())
        
    def insertselected(self):
        box2 = self._rightlist
        self._current = box2.index("active")
        oldname = box2.get(self._current)
        self._editdialog.setentry(oldname)
        #self._editdialog.component('entry').icursor('end')
        
    def doEdit(self, button):
        self._editdialog.deactivate()
        if button and button=="OK":
            box2 = self._rightlist
            newname = self._editdialog.get().strip()
            box2.delete(self._current)
            box2.insert(self._current, newname)
    
    def getright(self):
        return list(self._rightlist.get())
        
class DatabaseWidget(Pmw.MegaToplevel):
    def __init__(self, filtus):
        self.filtus = filtus
        Pmw.MegaToplevel.__init__(self, filtus.parent, master=filtus.parent, title="Variant databases", activatecommand=self._prepare)
        self.withdraw()
        self.interior().columnconfigure(0, weight=1)
        self.interior().rowconfigure(0, weight=1)
        self.notebook = Pmw.NoteBook(self.interior())
        self.notebook.add("New database", page_pyclass=DBpage_new, page_filtus=filtus, page_notebook=self)
        self.notebook.add("Add samples", page_pyclass=DBpage_add, page_filtus=filtus, page_notebook=self)
        self.notebook.add("Extract subset", page_pyclass=DBpage_extract, page_filtus=filtus, page_notebook=self)
        self.notebook.add("Search", page_pyclass=DBpage_search, page_filtus=filtus, page_notebook=self)
        self.notebook.setnaturalsize()
        self.notebook.grid(sticky='news')
    
    def _prepare(self):
        filenames = self.filtus.shortFileNameList
        for p in [self.notebook.page(y) for y in ("New database", "Add samples")]:
            p.lists.deselectAll()
            p.lists.setItems(filenames)
            p.lists.component('leftlist').settoptext("Available samples: %d" %len(filenames))
            
   
class DBpage_new(Tkinter.Frame):
    def __init__(self, parent, filtus, notebook):
        Tkinter.Frame.__init__(self, parent)
        self.filtus = filtus
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        
        Tkinter.Label(self, text="New variant database", font=filtus.titlefont, anchor='c').grid(sticky='news', padx=10, pady=15)
        
        filenames = [VF.shortName for VF in filtus.files]
        self.lists = TwoListWidget(self, filtus=filtus, lefttoptext="Available samples: %d"%len(filenames), items=filenames)
        self.lists.grid(sticky='news', padx=10, pady=0)
        
        ####### OUTPUT
        Tkinter.Frame(self, height=2, borderwidth=2, relief="sunken").grid(padx=10, pady=10, sticky="new")
        save_interior = Tkinter.Frame(self)
        save_interior.columnconfigure(0, weight=1)
        self.save_browser = FiltusWidgets.FileBrowser(save_interior, filtus=filtus, label="Output database file name:", 
                checkbutton = False, labelpos='nw', browsesticky='se', entryfield_entry_width=20, 
                browsetitle="")
        self.save_browser.browsebutton.configure(command = self._browseSave)
        self.save_browser.entryfield.configure(command = None)
        
        self.formatSelect = Pmw.RadioSelect(save_interior, buttontype="radiobutton", pady=0, 
                        labelpos="nw", labelmargin=0, label_text = "Format:")#, hull_borderwidth = 2, hull_relief = 'ridge')
        self.formatSelect.add("Simple")
        self.formatSelect.add("Extended")
        self.formatSelect.invoke("Extended")
        
        self.save_browser.grid(sticky='news', padx=(10,25), pady=(5,10))
        self.formatSelect.grid(row=0, column=1, sticky='e', padx=10, pady=0)
        save_interior.grid(sticky='news', padx=10, pady=0)
        
        #######
        Tkinter.Frame(self, height=2, borderwidth=3, relief="raised").grid(sticky='ew', pady=15)
        self.lowerButtons = Pmw.ButtonBox(self, padx=50)
        self.lowerButtons.add('createButt', text="Create database", command=self.createdb)
        self.lowerButtons.add('cancelButt', text="Cancel", command=notebook.deactivate)
        self.lowerButtons.grid(sticky='news', pady=(0,10))
        
    def _browseSave(self):
        fil = tkFileDialog.asksaveasfilename(initialdir=self.filtus.currentDir, title = "Save database as")
        if fil:
            self.filtus.currentDir = os.path.dirname(fil)
            self.save_browser.setvalue(os.path.normpath(fil))
        
    def createdb(self):
        try:
            self._doCreate()
        except Exception as e: 
            FiltusUtils.warningMessage(e)
           
    def _doCreate(self):
        st = time.time()
        filtus = self.filtus
        selection = self.lists.selection
        VFlist = [filtus.filteredFiles[i] for i in selection]
        sampleNames = self.lists.getright()
        outFormat = self.formatSelect.getvalue()
        outFilename = self.save_browser.getvalue()
        if not selection:
            raise RuntimeError("No samples selected")
        if not outFilename:
            raise RuntimeError("Please specify output file name and format")
        db = VariantDatabase.buildFromSamples(filtus, VFlist, outFormat, sampleNames)
        db.save(outFilename)
        message = "Variant database written to %s.\n\n" % outFilename \
                + "\n".join(db.infoSummary()) \
                + '\n\nTime used: %.2f seconds' %(time.time()-st,) 
        FiltusUtils.infoMessage(message)
        
class DBbrowser(Tkinter.Frame):
    def __init__(self, parent, filtus, label, updates=None, **kw):
        Tkinter.Frame.__init__(self, parent)
        self.columnconfigure(0, weight=1)
        self.filtus = filtus
        self.filename = None
        self.updates = updates
        self.nSamples, self.nVariants, self.format, self.colNames = None, None, '', []
        self.browser = FiltusWidgets.FileBrowser(self, filtus=filtus, label=label, 
                checkbutton=kw.get('checkbutton', False), labelpos=kw.get('labelpos', 'nw'), 
                browsesticky=kw.get('browsesticky', 'se'),
                entryfield_command=self.loadMeta_and_update, entryfield_entry_width=kw.get('width', 20), 
                browsetitle="Select variant database file")
        
        self.browser.grid(sticky='news', padx=0)
        self.summaryLabel = Tkinter.Label(self, text = "Summary:")
        self.summaryLabel.grid(sticky='w', padx=0, pady=(5,0))
        
    
    def getInfo(self):
        return (self.filename, self.nSamples, self.nVariants, self.format)
        
    def getColnames(self):
        return self.colNames
        
    def updateSummary(self, nSamples, nVariants, format):
        self.summary = "Summary:     Format = %s.     Samples = %s.     Variants = %s." %(format, str(nSamples), str(nVariants))
        self.summaryLabel.configure(text=self.summary)
        
    def loadMeta_and_update(self):
        self.filename = self.browser.getvalue()
        updates = self.updates
        try:
            meta, self.nSamples, self.nVariants, self.format, self.colNames = VariantDatabase.readMeta(self.filename)
            self.updateSummary(self.nSamples, self.nVariants, self.format)
            if updates: updates(inFormat=self.format, colNames=self.colNames)
        except Exception as e: 
            FiltusUtils.warningMessage("Could not load database.\n\n %s" %e)
            self.summaryLabel.configure(text="Summary: Error")
            self.filename = None
            if updates: updates(reset=True)
    
    def on(self):
        return self.browser.on()
    
    def modified(self):
        return self.browser.modified
    
    def setUnmodified(self):
        self.browser.modified = 0
        
class DBpage_add(Tkinter.Frame):
    def __init__(self, parent, filtus, notebook):
        Tkinter.Frame.__init__(self, parent)
        self.filtus = filtus
        self.db = None
        self.columnconfigure(0, weight=1)
        self.rowconfigure(3, weight=1)
        
        Tkinter.Label(self, text="Add samples to database", font=filtus.titlefont, anchor='c').grid(sticky='news', padx=10, pady=(15,10))
        
        self.browser = DBbrowser(self, filtus, label="Load database:", updates=self.updateStuff)
        self.browser.grid(sticky='news', padx=10)
        Tkinter.Frame(self, height=2, borderwidth=2, relief="sunken").grid(padx=10, pady=10, sticky="new")
        filenames = [VF.shortName for VF in filtus.files]
        self.lists = TwoListWidget(self, filtus=filtus, lefttoptext="Available samples: %d"%len(filenames), items=filenames)
        self.lists.grid(sticky='news', padx=10, pady=0)
        
        ####### OUTPUT
        Tkinter.Frame(self, height=2, borderwidth=2, relief="sunken").grid(padx=10, pady=10, sticky="new")
        save_interior = Tkinter.Frame(self)
        save_interior.columnconfigure(0, weight=1)
        self.save_browser = FiltusWidgets.FileBrowser(save_interior, filtus=filtus, label="Output database file name:", 
                checkbutton = False, labelpos='nw', browsesticky='se', entryfield_entry_width=20, 
                browsetitle="")
        self.save_browser.browsebutton.configure(command = self._browseSave)
        self.save_browser.entryfield.configure(command = None)
        
        self.formatSelect = Pmw.RadioSelect(save_interior, buttontype="radiobutton", pady=0, 
                        labelpos="nw", labelmargin=0, label_text = "Format:")#, hull_borderwidth = 2, hull_relief = 'ridge')
        self.formatSelect.add("Simple")
        self.formatSelect.add("Extended")
        self.formatSelect.invoke("Simple")
        
        self.save_browser.grid(sticky='news', padx=(10,25), pady=(5,10))
        self.formatSelect.grid(row=0, column=1, sticky='e', padx=10, pady=0)
        save_interior.grid(sticky='news', padx=10, pady=0)
        #######
        
        Tkinter.Frame(self, height=2, borderwidth=3, relief="raised").grid(sticky='ew', pady=15)
        self.lowerButtons = Pmw.ButtonBox(self, padx=50)
        self.lowerButtons.add('createButt', text="Add to database", command=self.addSamples)
        self.lowerButtons.add('cancelButt', text="Cancel", command=notebook.deactivate)
        self.lowerButtons.grid(sticky='news', pady=(0,10))
    
    def updateStuff(self, inFormat=None, colNames=None, reset=False):
        extButt = self.formatSelect.button("Extended")
        if not reset:
            if inFormat=="Simple": 
                self.formatSelect.invoke("Simple")
                extButt.configure(state="disabled")
            else:
                extButt.configure(state="normal")
        else: 
            extButt.configure(state="normal")
            
    def _browseSave(self):
        fil = tkFileDialog.asksaveasfilename(initialdir=self.filtus.currentDir, title = "Save database as")
        if fil:
            self.filtus.currentDir = os.path.dirname(fil)
            self.save_browser.setvalue(os.path.normpath(fil))
     
    def addSamples(self):
        try:
            self.doAddSamples()
        except Exception as e: 
            FiltusUtils.warningMessage(e)
        
    def doAddSamples(self):
        st = time.time()
        filtus = self.filtus
        inFilename, inNS, inNV, inFormat = self.browser.getInfo()
        outFilename = self.save_browser.getvalue()
        outFormat = self.formatSelect.getvalue()
        selection = self.lists.selection
        VFlist = [filtus.filteredFiles[i] for i in selection]
        sampleNames = self.lists.getright()
        if not inFilename:
            raise RuntimeError("Please specify existing database")
        if not selection:
            raise RuntimeError("No samples selected")
        if not outFilename:
            raise RuntimeError("Please specify output file name and format")
        newmeta = ''
        db = VariantDatabase.readFileAndAdd(filtus, inFilename, inFormat=inFormat, 
            inNS=inNS, outFormat=outFormat, VFlist=VFlist, sampleNames=sampleNames)
        if db.nSamples == inNS: raise IndexError("No samples to add")
        db.save(outFilename)
        message = "Variant database written to %s.\n\n" % outFilename \
                + "\n".join(db.infoSummary()) \
                + '\n\nTime used: %.2f seconds' %(time.time()-st,)
        FiltusUtils.infoMessage(message)
        
    
class DBpage_extract(Tkinter.Frame):
    def __init__(self, parent, filtus, notebook):
        Tkinter.Frame.__init__(self, parent)
        self.filtus = filtus
        self.db = None
        self.columnconfigure(0, weight=1)
        self.rowconfigure(3, weight=1)
        
        Tkinter.Label(self, text="Extract subset database", font=filtus.titlefont, anchor='c').grid(sticky='news', padx=10, pady=(15,10))
        
        self.lists = TwoListWidget(self, filtus=filtus, lefttoptext="Samples:")
        self.browser = DBbrowser(self, filtus, label="Load database:", updates=self.updateStuff)
        
        self.browser.grid(sticky='news', padx=10, pady=0)
        Tkinter.Frame(self, height=2, borderwidth=2, relief="sunken").grid(padx=10, pady=10, sticky="new")
        self.lists.grid(sticky='news', padx=10, pady=0)
        
        ### Column filter
        colnames = ('','CHROM','POS','OBS','HET','HOM','AFREQ')
        operatorNames = ["equal to", "not equal to", "greater than", "less than"]
        Tkinter.Label(self, text = "Variant filter:").grid(sticky='w', pady=(5, 2), padx=10)
        
        self.columnFilter = FiltusWidgets.ColumnFilter(self, filtus=self.filtus, colnames=colnames, operatorNames=operatorNames)
        self.columnFilter.component('checkbutton').grid_remove()
        self.columnFilter.grid(sticky='nw', padx=10)
        ####### OUTPUT
        Tkinter.Frame(self, height=2, borderwidth=2, relief="sunken").grid(padx=10, pady=10, sticky="new")
        save_interior = Tkinter.Frame(self)
        save_interior.columnconfigure(0, weight=1)
        self.save_browser = FiltusWidgets.FileBrowser(save_interior, filtus=filtus, label="Output database file name:", 
                checkbutton = False, labelpos='nw', browsesticky='se', entryfield_entry_width=20, 
                browsetitle="")
        self.save_browser.browsebutton.configure(command = self._browseSave)
        self.save_browser.entryfield.configure(command = None)
        
        self.formatSelect = Pmw.RadioSelect(save_interior, buttontype="radiobutton", pady=0, 
                        labelpos="nw", labelmargin=0, label_text = "Format:")#, hull_borderwidth = 2, hull_relief = 'ridge')
        self.formatSelect.add("Simple")
        self.formatSelect.add("Extended")
        self.formatSelect.invoke("Simple")
        
        self.save_browser.grid(sticky='news', padx=(10,25), pady=(5,10))
        self.formatSelect.grid(row=0, column=1, sticky='e', padx=10, pady=0)
        save_interior.grid(sticky='news', padx=10, pady=0)
        
        #######
        
        Tkinter.Frame(self, height=2, borderwidth=3, relief="raised").grid(sticky='ew', pady=15)
        self.lowerButtons = Pmw.ButtonBox(self, padx=50)
        self.lowerButtons.add('createButt', text="Extract subset", command=self.extractdb)
        self.lowerButtons.add('cancelButt', text="Cancel", command=notebook.deactivate)
        self.lowerButtons.grid(sticky='news', pady=(0,10))
    
    def updateStuff(self, inFormat=None, colNames=None, reset=False):
        lists = self.lists
        extButt = self.formatSelect.button("Extended")
        if not reset:
            if inFormat=="Simple": 
                lists.setItems([])
                lists.deselectAll()
                lists.setButtonStates(state="disabled")
                self.formatSelect.invoke("Simple")
                extButt.configure(state="disabled")
            else:
                lists.setButtonStates(state="normal")
                lists.setItems(colNames[(colNames.index('AFREQ') + 1):])
                lists.deselectAll()
                extButt.configure(state="normal")
        else: 
            lists.setButtonStates(state="normal")
            lists.setItems([])
            lists.deselectAll()
            extButt.configure(state="normal")
        
    def _browseSave(self):
        fil = tkFileDialog.asksaveasfilename(initialdir=self.filtus.currentDir, title = "Save database as")
        if fil:
            self.filtus.currentDir = os.path.dirname(fil)
            self.save_browser.setvalue(os.path.normpath(fil))
        
    # def extractdb(self):
        # try:
            # self.doExtractAndSave()
        # except Exception as e: 
            # FiltusUtils.warningMessage(e)
        
    def extractdb(self):
        try:
            st = time.time()
            inFilename, inNS, inNV, inFormat = self.browser.getInfo()
            subset = self.lists.selection
            sampleNames = self.lists.getright()
            outFilename = self.save_browser.getvalue()
            outFormat = self.formatSelect.getvalue()
            cfilter = self.columnFilter.getfilter()
            if cfilter:
                filter = Filter.Filter(columnfilters=[cfilter])
            else:
                filter = None
            if not inFilename:
                raise RuntimeError("Please specify existing database")
            if inFormat == "Extended" and not subset:
                raise RuntimeError("No samples selected")
            if not outFilename:
                raise RuntimeError("Please specify output file name and format")
            db = VariantDatabase.readFileAndExtract(self.filtus, inFilename, inFormat, inNS, 
                    subset, sampleNames, outFormat, filter=filter)
            
            db.save(outFilename)
            message = "Variant database written to %s.\n\n" % outFilename \
                    + "\n".join(db.infoSummary()) \
                    + '\n\nTime used: %.2f seconds' %(time.time()-st,)
            FiltusUtils.infoMessage(message)
        except Exception as e: 
            FiltusUtils.warningMessage(e)
            
class DBpage_search(Tkinter.Frame):
    def __init__(self, parent, filtus, notebook):
        Tkinter.Frame.__init__(self, parent)
        self.filtus = filtus
        self.db = None
        self.results = None
        self.columnconfigure(0, weight=1)
        self.rowconfigure(3, weight=1)
        
        Tkinter.Label(self, text="Search database", font=filtus.titlefont, anchor='c').grid(sticky='news', padx=10, pady=(15,10))
        
        self.browser = DBbrowser(self, filtus, label="Enter database file name")
        self.browser.grid(sticky='news', padx=10, pady=0)
        Tkinter.Frame(self, height=2, borderwidth=2, relief="sunken").grid(padx=10, pady=10, sticky="new")
        
        search_frame = Tkinter.Frame(self)
        search_frame.columnconfigure(1, weight=1)
        search_frame.rowconfigure(1, weight=1)
        sgroup = Pmw.Group(search_frame, tag_text="Query")
        search_int = sgroup.interior()
        self.chrom = Pmw.EntryField(search_int, label_text = "Chromosome:", labelpos='w', entry_width=9, labelmargin=5)
        self.pos = Pmw.EntryField(search_int, label_text = "Position:", labelpos='w', entry_width=9, labelmargin=5)
        Pmw.alignlabels([self.chrom, self.pos])
        self.chrom.grid(padx=(5,10), pady=(10,5))
        self.pos.grid(padx=(5,10), pady=5)
        sgroup.grid(sticky='new', padx=(0,10), pady=12)
        
        self.searchButt = Tkinter.Button(search_frame, text="Search", command=self.doSearch)
        self.searchButt.grid(row=1, column=0, padx=5, pady=20, sticky='n')
        
        self.resultWindow = Pmw.ScrolledText(search_frame, borderframe=1, text_padx=2,
                scrollmargin=2, hscrollmode='dynamic', label_text="Results", labelpos="nw",
                text_width=10, text_height=1, text_wrap='none', text_font=filtus.monofont)
        FiltusWidgets.makeReadOnly(self.resultWindow.component('text'))
        
        self.resultWindow.grid(row=0, column=1, sticky='news', padx=10, pady=0, rowspan=2)
        search_frame.grid(sticky='news', padx=10, pady=(0,10))
        
        Tkinter.Frame(self, height=2, borderwidth=3, relief="raised").grid(sticky='ew', pady=15)
        self.lowerButtons = Pmw.ButtonBox(self, padx=50)
        self.lowerButtons.add('createButt', text="Save result", command=self.save)
        self.lowerButtons.add('cancelButt', text="Cancel", command=notebook.deactivate)
        self.lowerButtons.grid(sticky='news', pady=(0,10))
        
    def save(self):
        if self.results is None: return
        filtus = self.filtus
        db_summary = '## Database file: %s\n## Format: %s\n## Number of samples: %d\n## Number of variants: %d\n##\n' %(self.filename, self.format, self.nSamples, self.nVariants)
        query = '## Query: Chromosome %s, position: %s' % tuple(self.query)
        meta = FiltusUtils.preambleNY(VFlist=None, analysis="VARIANT DATABASE - SEARCH\n##\n" + db_summary + query)
        
        fname = tkFileDialog.asksaveasfilename(initialdir=filtus.currentDir, title = "Save search results as")
        if not fname:
            return
        filtus.currentDir = os.path.dirname(fname)
        
        includePre = self.filtus.includePreamble
        try:
            with open(fname, 'w') as utfil:
                if 'Top' in includePre: utfil.write(meta)
                utfil.write(self.results)
                if 'Bottom' in includePre: utfil.write('\n' + meta)
        except Exception as e:
            FiltusUtils.warningMessage('%s\n\nFile not saved.' % e)
            return
    
    def _chromStartIndex(self, inFilename):
        with open(inFilename, 'rU') as database:
            first_index = {}
            last_index = {}
            curr_2char = '##'
            curr_chr = None
            for i,v in enumerate(database):
                if v[:2] != curr_2char: 
                    if curr_chr: last_index[curr_chr] = i-1
                    curr_chr = v[:2].strip()
                    first_index[curr_chr] = i
                    curr_2char = v[:2]
            last_index[curr_chr] = i
        return first_index, last_index
        
    def doSearch(self):
        self.results = None
        self.resultWindow.clear()
        chrom = self.chrom.getvalue().strip()
        pos = self.pos.getvalue().strip()
        inFilename, inNS, inNV, inFormat = self.browser.getInfo()
        if not all(x for x in (inFilename, chrom, pos)):
            if not inFilename: FiltusUtils.warningMessage("No database loaded")
            elif not chrom: FiltusUtils.warningMessage("Please indicate chromosome")
            elif not pos: FiltusUtils.warningMessage("Please indicate position")
            return
        query = self.query = [chrom, pos]
        self.filename = inFilename
        self.nSamples = inNS
        self.nVariants = inNV
        self.format = inFormat
        self.colNames = self.browser.getColnames()
        
        try:
            if self.browser.modified():
                self.firstIndex, self.lastIndex = self._chromStartIndex(inFilename)
                self.browser.setUnmodified()
            with open(self.filename, 'rU') as database:
                slice_database = itertools.islice(database, self.firstIndex[chrom], self.lastIndex[chrom])
                data = next((v for v in slice_database if v.split('\t')[:2] == query), None)
            if data is None:
                self.results = "Not found in the database"
                self.resultWindow.settext(self.results)
                return
            
            data = data.strip().split('\t')
            fields = ['Total observations', 'Heterozygous', 'Homozygous', 'Allele frequency in database']
            results = ['%s: %s' %x for x in zip(fields, data[2:6])]
            results[0] += ' (out of %d)' % inNS
            
            if self.format=="Extended":
                allSamples = self.colNames[-inNS:]
                allObs = map(int, data[-inNS:])
                observations = [(s, obs) for s, obs in zip(allSamples, allObs) if obs != 0]
                
                samples, gtCode = zip(*observations)
                gt = [('heterozygous', 'homozygous')[x-1] for x in gtCode]
                width = max(map(len, samples))
                
                results.extend(['','Samples:'] + [s.ljust(width) + ' - ' + g for s,g in zip(samples, gt)])
            
            self.results = '\n'.join(results)
            self.resultWindow.settext(self.results)
        
        except Exception as e:
            FiltusUtils.warningMessage(e)
            return
        

class VariantDatabase(DataContainer.ColumnData):
    def __init__(self, format, nSamples, columnNames, variantMatrix=None, extendedDict=None, filter=None, meta=''):
        '''Either variantMatrix or extendedDict should be None'''
        format = formatInit(format)
        if format=="S": columnNames = columnNames[:6]
        DataContainer.ColumnData.__init__(self, columnNames=columnNames, variants=variantMatrix, columnDescriptions=None)
        
        if variantMatrix is None:
            vars = self.convertFromDict(extendedDict, format, 2.0*nSamples)
            self.setVariants(vars)
        
        if filter: 
            filter.apply(self, checks = False, inplace=True)
        
        self.sort(column=self.columnNames[1], descending=False)
        self.sort(column=self.columnNames[0], descending=False)
        
        self.format = format
        self.nSamples = nSamples
        self.sampleNames = columnNames[6:]
        if format=="E" and len(self.sampleNames) != nSamples:
            raise IndexError("Number of samples (%d) not consistent with the given sample names:\n" % nSamples +'\n'.join(self.sampleNames))
        self.meta = meta
    
        
    @classmethod
    def loadSimple(cls, filename, nSamples, filter=None, keepMeta=True, newMeta=''):
        with open(filename, 'rU') as dbfil:
            oldmeta, colNames, dbfil = _readTop(dbfil)
            if len(colNames)==6: 
                # if the database base actually has Simple format.
                variants = [line.strip().split('\t') for line in dbfil]
            else:
                # if the database is extended, extract first 6 columns to make it simple.
                variants = [line.strip().split('\t')[:6] for line in dbfil]
        
        meta = '\n'.join(oldmeta) if keepMeta else ''
        meta += newMeta
        
        return cls("Simple", nSamples, colNames, variantMatrix=variants, filter=filter, meta=meta)
        
    @classmethod
    def buildFromSamples(cls, VFlist, outFormat, sampleNames=None):
        r"""Build database.

        If format is "Simple", the database has 6 columns by default: 
        CHROM 
        POS 
        OBS (# samples with a variant at this position)
        HET (# samples with het variant),
        HOM (# samples with hom variant)
        AFREQ (allele frequency)
        
        In extended format, the first 6 columns are as above, followed by one column per sample. 
        Entries in these columns are 0, 1 or 2 (= not present, het, hom)
        """
        
        N = len(VFlist)
        # First creating the extended part, i.e. a matrix (list of lists) with 0,1,2. 
        # (This turned out to be quicker than getUniqueVariants(), also for simple format.)
        extended = collections.defaultdict(lambda : [0]*N)
        for i,VF in enumerate(VFlist):
            gt = VF.GTnum()
            vDef = VF.varDefGetter
            for v in VF.variants:
                extended[vDef(v)][i] = gt(v)
        
        colNames = ['CHROM','POS', 'OBS', 'HET', 'HOM', 'AFREQ']
        
        outFormat = formatInit(outFormat)
        if outFormat=="E":
            if sampleNames is None: sampleNames = [VF.shortName for VF in VFlist]
            colNames += sampleNames
        meta = FiltusUtils.preambleNY(VFlist=VFlist, analysis = "NEW VARIANT DATABASE", sort=False)
        
        return cls(outFormat, nSamples=N, columnNames=colNames, extendedDict=extended, meta=meta)
    
    @classmethod
    def readFileAndAdd(cls, filename, inFormat, inNS, outFormat, VFlist, sampleNames=None):
        inFormat = formatInit(inFormat)
        outFormat = formatInit(outFormat)
        
        if sampleNames is None: sampleNames = [VF.shortName for VF in VFlist]
        if inFormat == "S":
            old = VariantDatabase.loadSimple(filename, inNS)
            new = VariantDatabase.buildFromSamples(VFlist, outFormat, sampleNames)
            meta = FiltusUtils.preambleNY(VFlist=VFlist, analysis="ADDED TO DATABASE", sort=False, appendTo=old.meta)
            return old.addSimple(new, meta=meta)
        
        N_add = len(VFlist)
        with open(filename, 'rU') as dbfil:
            meta, colNames, dbfil = _readTop(dbfil)
            
            OBS_ind = colNames.index('OBS')
            v_extractor = itemgetter(*range(OBS_ind))
            
            sample_ind = OBS_ind + 4
            N = inNS + N_add
            extended = collections.defaultdict(lambda : [0]*N)
            for line in dbfil:
                dat = line.strip().split('\t')
                vdef = v_extractor(dat)
                extended[vdef][:inNS] = map(int, dat[sample_ind:])
        
        simpleColnames = colNames[:sample_ind]
        sampleNames = colNames[sample_ind:] + sampleNames
        k = inNS
        skip = []
        for VF in VFlist: # not enumerate() since the sample might be skipped
            gt = VF.GTnum()
            vDef = VF.varDefGetter
            for v in VF.variants:
                extended[vDef(v)][k] = gt(v)
            
            # check if already included
            copy = next((j for j in range(inNS) if all(x[j]==x[k] for x in extended.itervalues())), None)
            if copy is None:
                k += 1
                continue
            else:    
                question = "The new sample '%s' is exactly equal to '%s' in the database." % (sampleNames[k], sampleNames[copy]) \
                        + "\n\nSkip this sample?"
                action = FiltusUtils.yesnoMessage(question)
                if action:
                    skip.append(VF)
                    for v in extended:
                        del extended[v][k]
                    del sampleNames[k]
                    extended.default_factory = lambda: [0]*len(sampleNames)
        
        new_meta = FiltusUtils.preambleNY(VFlist = [VF for VF in VFlist if not VF in skip], sort=False, 
                                          analysis = "ADDED TO DATABASE", appendTo='\n'.join(meta))
        colNames = simpleColnames + sampleNames
        return cls(outFormat, nSamples=len(sampleNames), columnNames=colNames, extendedDict=extended, meta=new_meta)
       
    @classmethod
    def readFileAndExtract(cls, filename, inFormat, inNS, subset, sampleNames, outFormat, filter):
        inFormat = formatInit(inFormat)
        outFormat = formatInit(outFormat)
        if inFormat == "S" or (outFormat=="S" and inNS == len(subset)):
            return VariantDatabase.loadSimple(filename, nSamples=inNS, filter=filter, keepMeta=False)
        
        extended = {}
        with open(filename, 'rU') as dbfil:
            meta, colNames, dbfil = _readTop(dbfil)
            
            OBS_ind = colNames.index('OBS')
            v_extractor = itemgetter(*range(OBS_ind))
            gt_extractor = itemgetter(*[s + OBS_ind + 4 for s in subset])
            
            for line in dbfil:
                dat = line.strip().split('\t')
                gt = map(int, gt_extractor(dat)) # map() makes a list, also of singletons!
                if any(gt):
                    extended[v_extractor(dat)] = gt
        
            colNames = colNames[:OBS_ind + 4] + sampleNames
        
        return cls(outFormat, len(subset), colNames, extendedDict=extended, filter=filter)
        
    def addSimple(self, other, meta):
        '''Adds a simple database to this (simple) database. Returns a new database.'''
        vdef = itemgetter(0,1)
        dic = {vdef(v) : map(int, v[3:5]) for v in self.variants}
        for v in other.variants:    
            vv = vdef(v)
            if vv in dic:
                dic[vv][0] += v[3]
                dic[vv][1] += v[4]
            else:
                dic[vv] = v[3:5]
        
        def stats(hethom, Ndip):
            '''gt: vector with genotype codes 0,1,2'''
            het = hethom[0]
            hom = hethom[1]
            freq = round((het+2*hom)/Ndip, 4)
            return [het+hom, het, hom, freq]
        
        N = self.nSamples + other.nSamples
        Ndip = N*2.0
        newvars = [list(v) + stats(hethom, Ndip) for v, hethom in dic.iteritems()]

        return VariantDatabase("Simple", nSamples=N, columnNames=self.columnNames[:], variantMatrix=newvars, meta=meta)
        
        
    def convertFromDict(self, dic, format, Ndip):
        '''takes dict in extended format into matrix, adding stats columns'''
        
        format = formatInit(format)
        
        def stats(gtvec, Ndip):
            '''gt: vector with genotype codes 0,1,2'''
            het = sum(d==1 for d in gtvec)
            hom = sum(d==2 for d in gtvec)
            freq = round((het+2*hom)/Ndip, 4)
            return [het+hom, het, hom, freq]
        
        if format=="S":
            return [list(v) + stats(gtvec, Ndip) for v, gtvec in dic.iteritems()]
        else:    
            return [list(v) + stats(gtvec, Ndip) + gtvec for v, gtvec in dic.iteritems()]
 
    @staticmethod
    def readMeta(filename):
        allvars = []
        with open(filename, 'rU') as dbfil:
            meta, colNames, dbfil = _readTop(dbfil)
            
            if not colNames or len(colNames) < 5:
                raise RuntimeError("I don't think this is a database file.")
            if any(col not in colNames for col in ('OBS', 'HET', 'HOM', 'AFREQ')):
                raise RuntimeError("Database error: Column names must include OBS, HET, HOM and AFREQ.")
            if len(meta) < 4:
                allvars = [line for line in dbfil]
                
        afreq_ind = colNames.index('AFREQ')
        sampleNames = colNames[afreq_ind + 1:]
        
        if len(meta) >= 4:
            general = '\n'.join(meta[:-4]) + '\n'
            nSamples = int(meta[-4].split(':')[1])
            nVariants = int(meta[-3].split(':')[1])
            format = meta[-2].split(':')[1].strip()
   
            if format == "Simple" and len(sampleNames) > 0:
                raise RuntimeError("Database error: Format is inconsistent with column names.")
            if format == "Extended" and nSamples != len(sampleNames):
                raise RuntimeError("Number of samples is inconsistent with column names")
        else:
            general = ''
            if len(sampleNames)==0:
                format = "Simple"
                v = allvars[0].strip().split('\t')[:afreq_ind + 1]
                nSamples = round((int(v[-3]) + 2*int(v[-2]))/(2*float(v[-1])))
            else:
                fornat = "Extended"
                nSamples = len(sampleNames)
            nVariants = len(allvars)
                
        return general, nSamples, nVariants, format, colNames
     
    
    def infoSummary(self):
        return ['Samples: %d' % self.nSamples, 'Variants: %d'% self.length, 'Format: %s' % self.format]
        
    def save(self, outFilename, sep='\t', colnames=True, preamblePos=("Top",)): 
        info = self.infoSummary()
        meta = self.meta + '##  ' + '\n##  '.join(info) + '\n##\n'
        DataContainer.ColumnData.save(self, outFilename, sep=sep, colnames=colnames, preamblePos=preamblePos)
        
def formatInit(format):
    '''return the uppercase of the first letter of 'format'. Raises error if not 'S' (simple) or 'E' (extended)'''
    format0 = format[0].upper()
    if not format0 in ['S', 'E']: 
        raise ValueError("Wrong database format indicated: %s\nShould be either 'Simple' or 'Extended' (possibly abbreviated)." %format)
    return format0    
    
def okForDB(VFlist, ndef=None):
    message = ""
    if ndef is not None and any(len(VF.varDefColNames) != ndef for VF in VFlist):
        message = 'The existing database has %d variant-defining columns, but at least one of the selected samples does not match this. \
To extend this database, make sure to indicate matching columns in the "Columns uniquely defining a variant" entry when loading new files.' % ndef
    elif len(set(len(VF.varDefColNames) for VF in VFlist)) > 1:
        message = 'The selected files do not have the same number of variant-defining columns. To create the database, please load the files again, making sure the "Columns uniquely defining a variant" entries match.'
    elif any(VF.varDefGetter is None for VF in VFlist):
        message = 'There is a problem with the variant-defining columns of (at least one of) the selected samples.'
    if message:
        FiltusUtils.warningMessage(message)
        return False
    else:
        return True

    def readDB(self, filename):
        with open(filename, 'rU') as dbfil:
            db = [line.strip().split('\t') for line in dbfil]
        
        m = next(i for i in xrange(len(old)) if len(old[i]) > 1)
        meta = '\n'.join(old[:m]) + '\n'
        old[:] = old[m:]
    
        len_old = len(old) - 1
        old_heads = old[0]
        type = 'simple' if old_heads[-3:] == ['Total','Heterozygous','Homozygous'] else 'extended'

        
        
if __name__=='__main__':
    import VariantFileReader
    reader = VariantFileReader.VariantFileReader()

    vcftest = "C:\\Projects\\FILTUS\\Testfiles\\vcf_example2.vcf"
    VFlist = reader.readVCFlike(vcftest, sep="\t", chromCol="CHROM", posCol="POS", geneCol="Gene_INFO", splitAsInfo="INFO", keep00=0)
    dbe = VariantDatabase.buildFromSamples(VFlist, "Ex")
    dbe.save("kast.txt")
    dbs1 = VariantDatabase.readFileAndAdd("kast.txt", inFormat="E", inNS=3, outFormat="S", VFlist=VFlist)
    dbs1.save("kast2.txt")
    dbs2 = VariantDatabase.buildFromSamples(VFlist, "S")
    