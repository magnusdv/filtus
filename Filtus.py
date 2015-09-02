# -*- coding: utf-8 -*-
from Tkinter import *
import Pmw
import tkFont
import tkFileDialog

#import cProfile
#import pstats

import gc
import sys
import os.path
from math import copysign

import FiltusWidgets
import FiltusDatabase
import FiltusAnalysis
import FiltusUtils
import InputDialog

# module FiltusQC is imported in a special method - to check availability of numpy and matplotlib.
try:
    import FiltusQC
    PLOT_error = None
except ImportError as e:
    print 'plot import error'
    PLOT_error = e
 
class Filtus(object):
    def __init__(self, parent, version):
        self.parent = parent
        self.version = version
        self.busyManager = BusyManager(parent)
        self.windowingsystem = parent.tk.call('tk', 'windowingsystem')
        self.rightclickevents = ['<2>', '<Control-1>'] if self.windowingsystem == 'aqua' else ['<3>']
        parent.rowconfigure(1, weight=1)
        parent.columnconfigure(0, weight=1)
        
        self.scrollframe = Pmw.ScrolledFrame(parent, borderframe=0, clipper_borderwidth=0, vertflex='expand', horizflex='expand')
        frame = self.scrollframe.interior()
        frame.rowconfigure(2, weight=1)
        frame.columnconfigure(1, weight=1)
        self.frame = frame
    
        ### fonts
        self.defaultfont = tkFont.nametofont("TkDefaultFont")
        self.smallfont, self.smallbold, self.tinyfont, self.titlefont = self.defaultfont.copy(), self.defaultfont.copy(), self.defaultfont.copy(), self.defaultfont.copy()
        self.smallbold['weight'] = 'bold'
        self.monofont = tkFont.nametofont("TkFixedFont")
        self.monobold = self.monofont.copy()
        self.monobold['weight'] = 'bold'
        
        self.textfont = tkFont.nametofont("TkTextFont")
        #self.menufont = tkFont.nametofont("TkMenuFont") # this didn't respond to change... Used workaround by setting menu label fonts manually to defaultfont.
        for opt in ['family','weight','slant','underline','overstrike']:
            self.textfont[opt] = self.monofont[opt]
        self.setFontSizes(self.defaultfont['size'], self.textfont['size'], init=True)

        self.files = []
        self.filteredFiles = []
        self.shortFilenames = False
        self.longFileNameList = []
        self.shortFileNameList = []
        self.currentFileNameList = []
        
        self.currentDir = ""
        self.currentFileDir = ""
        self.storage = {} # storage for  variant databases (to avoid reloading when filtering)
        
        ############## The file group
        self.fileGroup = Frame(frame)
        self.fileGroup.columnconfigure(0, weight=1)

        self.fileListbox = FiltusWidgets.LabeledListBox(self.fileGroup, filtus=self, toptext="Loaded files: 0", width=50)
        self.fileSummary1 = FiltusWidgets.SummaryBox(self.fileGroup, filtus=self, toptext="Unfiltered summaries", width=36)
        self.fileSummary2 = FiltusWidgets.SummaryBox(self.fileGroup, filtus=self, toptext="Filtered summaries", width=36)

        self.fileListbox.grid(sticky='new')
        self.fileSummary1.grid(row=0, column=1, sticky='nw', padx=(10, 0))
        self.fileSummary2.grid(row=0, column=2, sticky='nw', padx=(10, 0))
        self.fileGroup.grid(sticky='news', columnspan=2, pady=(0, 10))

        ############## The filter group
        self.FM = FiltusWidgets.FilterMachine(frame, filtus=self)
        self.FM.grid(sticky='new')

        ############ Sharing notebook
        self.sharingNotebook = Pmw.NoteBook(frame, arrownavigation=False, pagemargin=0)

        self.gs = FiltusWidgets.GeneSharingPage(self.sharingNotebook, self, 'Gene sharing')
        self.fs = FiltusWidgets.GeneSharingPage(self.sharingNotebook, self, 'Gene sharing fam', family=True)
        self.vs = FiltusWidgets.VariantSharingPage(self.sharingNotebook, self, 'Variant sharing')
        self.dn = FiltusWidgets.DeNovoPage(self.sharingNotebook, self, 'De novo')
        
        self.sharingNotebook.setnaturalsize()
        self.sharingNotebook.grid(row=2, sticky='new', pady=(10, 0))

        ############## The big text field
        self.text = FiltusWidgets.FiltusText(frame, filtus=self)
        self.text.grid(row=1, column=1, rowspan=2, sticky='news', padx=(20, 0), pady=(10, 0))

        ###### Settings
        self.fileListbox.fixselectall()
        self.sepOutput = '\t'
        self.truncate = 50
        self.makeSettingsDialog()
        self.settingsDialog.invoke()

        self.makeMainMenu()
        self.menuBar.bind('<Triple-1>', self._run)
        self.parent.bind('<Shift-Return>', self._run)

        self.scrollframe.component('clipper').configure(height=min(frame.winfo_height(), frame.winfo_screenheight()-50),
                                                    width=min(frame.winfo_width(), frame.winfo_screenwidth()-50))
        self.scrollframe.grid(row=1, sticky='news', padx=20, pady=(10, 20))
        if PLOT_error:
            FiltusUtils.warningMessage("Plotting functionality is disabled. Error message:\n\n%s\n\nNote: On MAC and Linux the modules 'numpy' and 'matplotlib' must be installed separately to make the plots work.\nSee also Filtus homepage: folk.uio.no/magnusv/filtus.html"%PLOT_error)
        
             
    def _run(self, event):
        import autorun
        reload(autorun)

        def execute(button):
            if button == 'Cancel': prompt.deactivate()
            else: prompt.deactivate(prompt.get().strip())

        if event.keysym == 'Return':
            f = 'x'
        else:
            prompt = Pmw.PromptDialog(self.parent, buttons=('OK', 'Cancel'), title='Run command',
                label_text='Function to run:', entryfield_labelpos='n', command=execute, defaultbutton=0)
            f = FiltusUtils.activateInCenter(self.parent, prompt)
        if f: getattr(autorun, f, FiltusUtils.ignore_break)(self)


    def makeSettingsDialog(self):
        self.settingsDialog = Pmw.Dialog(self.parent, buttons=('OK', 'Cancel'), defaultbutton='OK', title='Preferences',
                                                    command=self.setSettings, dialogchildsite_pady=10, buttonbox_pady=10)
        self.settingsDialog.withdraw()
        self.settingsDialog.interior().columnconfigure(0, weight=1)
        pmw_OPTIONS = dict(labelmargin=10, labelpos='w')
        grid_OPTIONS = dict(sticky='news', padx=10, pady=10)

        ## Output file settings
        self.outputGroup = Pmw.Group(self.settingsDialog.interior(), tag_text='Output files')
        self.sepOutputOM = Pmw.OptionMenu(self.outputGroup.interior(), label_text="Column separator:",
                                    menubutton_padx=0, menubutton_pady=1, menubutton_width=10, items=["comma", "tab", "semicolon", "space"],
                                    initialitem='tab', **pmw_OPTIONS)
        self.preambleCheck = Pmw.RadioSelect(self.outputGroup.interior(), labelpos='w', label_text="Attach meta information:", buttontype='checkbutton', pady=0)
        self.preambleCheck.add("Top"); self.preambleCheck.add("Bottom")
        self.preambleCheck.setvalue(["Top"])
        self.sepOutputOM.grid(**grid_OPTIONS)
        self.preambleCheck.grid(**grid_OPTIONS)

        ## Font settings
        self.fontGroup = Pmw.Group(self.settingsDialog.interior(), tag_text='Font sizes')
        counter_OPTIONS = dict(entry_width=3, entry_justify='center', datatype='integer', entryfield_value=abs(self.textfont['size']),
                               entryfield_validate=dict(validator='integer', min='1', max='16'), **pmw_OPTIONS)
        self.generalfontSizeEntry = Pmw.Counter(self.fontGroup.interior(), label_text="General font size:", **counter_OPTIONS)
        self.textfontSizeEntry = Pmw.Counter(self.fontGroup.interior(), label_text="Font size for text fields:", **counter_OPTIONS)
        self.generalfontSizeEntry.grid(**grid_OPTIONS)
        self.textfontSizeEntry.grid(**grid_OPTIONS)

        self.otherGroup = Pmw.Group(self.settingsDialog.interior(), tag_text='Other settings')
        counter_OPTIONS.update(dict(entryfield_value=self.truncate, entryfield_validate=dict(validator='integer', min='0', max='1000')))
        self.truncateEntry = Pmw.Counter(self.otherGroup.interior(), label_text="Truncate columns wider than:", **counter_OPTIONS)
        self.truncateEntry.grid(**grid_OPTIONS)
        for g in [self.outputGroup, self.fontGroup, self.otherGroup]:
            g.interior().columnconfigure(0, weight=1)
            g.configure(ring_borderwidth=1)
            g.configure(tag_font = self.smallbold)
            g.grid(**grid_OPTIONS)

        Pmw.alignlabels([self.sepOutputOM, self.preambleCheck, self.generalfontSizeEntry, self.textfontSizeEntry, self.truncateEntry], sticky='w')

    def settingsPrompt(self):
        FiltusUtils.activateInCenter(self.parent, self.settingsDialog)

    ### Menus

    def makeMainMenu(self):
        menuBar = Pmw.MenuBar(self.parent, hull_relief='raised', hull_borderwidth=2)
        menuBar.grid(row=0, column=0, sticky='ew')
        menuBar.addmenu('File', None)
        menuBar.addmenuitem('File', 'command', None, command=self.browseFiles, label='Load variant files (simple)', font=self.defaultfont)
        menuBar.addmenuitem('File', 'command', None, command=self.advLoad_prompt, label='Load variant files (advanced)', font=self.defaultfont)
        menuBar.addmenuitem('File', 'separator')
        menuBar.addmenuitem('File', 'command', None, command=self.clearAll, label='Unload all samples', font=self.defaultfont)
        menuBar.addmenuitem('File', 'command', None, command=self.unloadFiles, label='Unload selected samples', font=self.defaultfont)
        menuBar.addmenuitem('File', 'separator')
        menuBar.addmenuitem('File', 'command', None, command=self.pedwriter_prompt, label='Export to MERLIN format (ped/map/dat/freq)', font=self.defaultfont)
        menuBar.addmenuitem('File', 'separator')
        menuBar.addmenuitem('File', 'command', None, command=self.text.save, label='Save main window content', font=self.defaultfont)
        menuBar.addmenuitem('File', 'separator')
        menuBar.addmenuitem('File', 'command', None, command=self.parent.destroy, label="Quit", font=self.defaultfont)

        menuBar.addmenu('View', None)
        menuBar.addcascademenu('View', 'columnsum', label="Summarize column", font=self.defaultfont)

        menuBar.addmenuitem('View', 'command', None, command=self.geneLookup_prompt, label="Gene lookup", font=self.defaultfont)
        menuBar.addmenuitem('View', 'separator')
        menuBar.addmenuitem('View', 'command', None, command=self.mergeCommand(collapse=False), label='Joint view of selected samples', font=self.defaultfont)
        menuBar.addmenuitem('View', 'command', None, command=self.mergeCommand(collapse=True), label='Joint view - unique variants only', font=self.defaultfont)
        menuBar.addmenuitem('View', 'separator')
        menuBar.addmenuitem('View', 'command', None, command=self.toggleShortNames, label='Toggle long/short sample names', font=self.defaultfont)
        
        menuBar.addmenu('Filters', None)
        menuBar.addmenuitem('Filters', 'command', None, command=self.FM.addColumnFilter, label="Add column filter", font=self.defaultfont)
        menuBar.addmenuitem('Filters', 'command', None, command=self.FM.removeColumnFilter, label="Remove column filter", font=self.defaultfont)
        menuBar.addmenuitem('Filters', 'separator')
        menuBar.addmenuitem('Filters', 'command', None, command=self.FM.unapplyFilters, label='Unapply (but keep) filters', font=self.defaultfont)
        menuBar.addmenuitem('Filters', 'command', None, command=self.FM.clearFilters, label='Clear all filters', font=self.defaultfont)
        menuBar.addmenuitem('Filters', 'command', None, command=self.FM.loadFilterBrowse, label='Load filter configuration', font=self.defaultfont)
        menuBar.addmenuitem('Filters', 'command', None, command=self.FM.saveFilter, label='Save current filter configuration', font=self.defaultfont)
        menuBar.addmenuitem('Filters', 'separator')
        menuBar.addmenuitem('Filters', 'command', None, command=self.multiFilter_prompt, label='Individual filtering', font=self.defaultfont)
        
        menuBar.addmenu('Analysis', None)
        menuBar.addmenuitem('Analysis', 'command', None, command=self.plink_prompt, label="Runs of homozygosity - PLINK", font=self.defaultfont)
        menuBar.addmenuitem('Analysis', 'command', None, command=self.autozyg_prompt, label="Regions of autozygosity - AutEx", font=self.defaultfont)
        menuBar.addmenuitem('Analysis', 'separator')
        menuBar.addmenuitem('Analysis', 'command', None, command=self.pairwiseSharing, label='Pairwise variant sharing', font=self.defaultfont)
        menuBar.addmenuitem('Analysis', 'separator')
        menuBar.addmenuitem('Analysis', 'command', None, command=self.QC_prompt, label='QC plots', font=self.defaultfont) 
        if PLOT_error:
            menuBar.component('Analysis-menu').entryconfigure('end', state='disabled')

        menuBar.addmenu('Database', None)
        menuBar.addmenuitem('Database', 'command', None, command=self.database_prompt(0), label='New database', font=self.defaultfont)
        menuBar.addmenuitem('Database', 'command', None, command=self.database_prompt(1), label='Add to database', font=self.defaultfont)
        menuBar.addmenuitem('Database', 'command', None, command=self.database_prompt(2), label='Extract subset', font=self.defaultfont)
        menuBar.addmenuitem('Database', 'command', None, command=self.database_prompt(3), label='Search database', font=self.defaultfont)

        menuBar.addmenu('Settings', None)
        menuBar.addmenuitem('Settings', 'command', None, command=self.genelengths_prompt, label='Gene lengths', font=self.defaultfont)
        menuBar.addmenuitem('Settings', 'separator')
        menuBar.addmenuitem('Settings', 'command', None, label='Preferences', command=self.settingsPrompt, font=self.defaultfont)
        self.menuBar = menuBar

    def busy(self):
        self.busyManager.busy()
    
    def notbusy(self):
        self.busyManager.notbusy()
        
    def checkLoadedSamples(self, select, minimum=None, maximum=None, VF=True, filtered=True):
        '''Typically called imediately after a button press starting analysis. 
        Checks that there are sufficient loaded samples for analysis, and returns either indices or VF objects.
        If any problems occurs, a warning is displayed and False is returned.
        Select = either "selection" or "all"
        '''
        def plural_s(k):
            return '' if k==1 else 's'
        
        useMin, useMax = minimum is not None, maximum is not None
        try:
            if len(self.files)==0:
                raise IndexError("No samples are loaded")
            files = self.filteredFiles if filtered else self.files
            if len(files) == 0:
                files = self.filteredFiles = self.filteredFiles_initialcopy()
            
            if select =="all":
                if useMin and len(files) < minimum:
                    raise IndexError("This option requires at least %d loaded sample%s." % (minimum, plural_s(minimum)))
                return files if VF else range(len(files))
            
            if select == "selection":
                seleci = [int(i) for i in self.fileListbox.curselection()]
                nsel = len(seleci)
                if useMin and useMax and not minimum <= nsel <= maximum:
                    if minimum==maximum:
                        raise IndexError("Please select exactly %d sample%s in the 'Loaded samples' window" %(minimum, plural_s(minimum)))
                    else:
                        raise IndexError("Please select between %d and %d samples in the 'Loaded samples' window" %(minimum, maximum))
                elif useMin and nsel < minimum:
                    raise IndexError("Please select at least %d sample%s in the 'Loaded samples' window" %(minimum, plural_s(minimum)))
                elif useMax and nsel > maximum:
                    raise IndexError("Please select at most %d sample%s in the 'Loaded samples' window" %(maximum, plural_s(maximum)))
            
            return [files[i] for i in seleci] if VF else seleci    
        
        except Exception as e:
            FiltusUtils.warningMessage(e)
            return False
        
    def autozyg_prompt(self):
        if not self.checkLoadedSamples(select="selection", VF=False, minimum=1, maximum=1):
            return
        if not hasattr(self, 'autexgui'):
            self.autexgui = FiltusWidgets.AutEx_GUI(self)
        try:
            FiltusUtils.activateInCenter(self.parent, self.autexgui)
        except Exception as e:
            FiltusUtils.warningMessage("%s: %s" %(type(e).__name__, e))
    
    def plink_prompt(self):
        if not self.checkLoadedSamples(select="selection", VF=False, minimum=1, maximum=1):
            return
        if not hasattr(self, 'plinkgui'):
            self.plinkgui = FiltusWidgets.PLINK_GUI(self)
        FiltusUtils.activateInCenter(self.parent, self.plinkgui)

    def pedwriter_prompt(self):
        if not self.checkLoadedSamples(select="all"):
            return
        if not hasattr(self, 'pedwriter'):
            self.pedwriter = FiltusWidgets.PedWriter(self)
        try:
            FiltusUtils.activateInCenter(self.parent, self.pedwriter)
        except Exception as e:
            FiltusUtils.warningMessage("%s: %s" %(type(e).__name__, e))
    
    def advLoad_prompt(self):
        if not hasattr(self, 'advLoad'):
            self.advLoad = FiltusWidgets.AdvancedLoad(self)
        FiltusUtils.activateInCenter(self.parent, self.advLoad.dialog)

    def database_prompt(self, pgnr):
        def _f():
            if not hasattr(self, 'databaseTool'):
                self.databaseTool = FiltusDatabase.DatabaseWidget(self)
            self.databaseTool.notebook.selectpage(pgnr)
            FiltusUtils.activateInCenter(self.parent, self.databaseTool)
        return _f
        
    def QC_prompt(self):
        if not self.checkLoadedSamples(select="all", minimum=1):
            return
        if not hasattr(self, 'QC'):
            self.QC = FiltusQC.QC(self)
        try:
            FiltusUtils.activateInCenter(self.parent, self.QC.dialog)
        except Exception as e:
            print "Killing myself because of: %s"%e
            self.parent.destroy()
    
    def geneLookup_prompt(self):
        if not self.checkLoadedSamples(select="all"):
            return
        if all(VF.geneGetter is None for VF in self.filteredFiles):
            FiltusUtils.warningMessage("None of the loaded samples have known gene column.")
            return
        if not hasattr(self, 'geneLookup'):
            self.geneLookup = FiltusWidgets.GeneLookup(self.parent, self)
        FiltusUtils.activateInCenter(self.parent, self.geneLookup.prompt)

    def multiFilter_prompt(self):
        if not self.checkLoadedSamples(select="all", minimum=2, filtered=False):
            return
        if not hasattr(self, 'multiFilter'):
            self.multiFilter = FiltusWidgets.MultiFilter(self)
        self.multiFilter.showDialog()

    def genelengths_prompt(self):
        if not hasattr(self, 'genelengths'):
            self.genelengths = FiltusWidgets.GeneLengthFile(self)
        self.genelengths.browse()

    def find_genmapfile(self): # Find and parse file with recombination map
        decodeFile = "DecodeMap_thin.txt"
        decode_trypaths = [os.path.join(direc, decodeFile) for direc in ("..", "example_files", "../example_files")]
        return next((p for p in decode_trypaths if os.path.isfile(p)), None)
        
    def mergeCommand(self, collapse):
        def _m():
            VFlist = self.checkLoadedSamples(select="selection", minimum=2)
            if not VFlist: return
            try:
                mergedVF = FiltusAnalysis.merge(VFlist, collapse=collapse)
                self.text.prettyPrint(mergedVF)
            except Exception as e:
                FiltusUtils.warningMessage(e)
        return _m
        
    def pairwiseSharing(self):
        seleci = self.checkLoadedSamples(select="selection", minimum=2, VF=False)
        if seleci: FiltusAnalysis.pairwiseSharing(seleci, filtus=self)
        
    def browseFiles(self):
        longFileNameList = tkFileDialog.askopenfilenames(initialdir=self.currentFileDir, title="Load files")
        longFileNameList = list(self.parent.tk.splitlist(longFileNameList))
        if longFileNameList:
            longFileNameList = [os.path.normpath(p) for p in longFileNameList]
            self.currentFileDir = os.path.dirname(longFileNameList[0])
            self.loadFiles(longFileNameList)

    def setSettings(self, result):
        if result != 'OK':
            self.settingsDialog.deactivate()
            return
        sepDic = {'(automatic)':'auto', 'comma':',', 'tab':'\t', 'semicolon':';', 'space':' '}
        self.sepOutput = sepDic[self.sepOutputOM.getvalue().strip()]
        self.includePreamble = self.preambleCheck.getvalue()
        defaultsize = int(self.generalfontSizeEntry.getvalue())
        textsize = int(self.textfontSizeEntry.getvalue())
        self.setFontSizes(defaultsize, textsize)
        self.truncate = int(self.truncateEntry.getvalue())
        self.settingsDialog.deactivate()

    def setFontSizes(self, defaultsize, textsize, init=False):
        """Set font sizes used by various widget classes"""
        olddef, oldtext = self.defaultfont['size'], self.textfont['size']
        defaultsize, textsize = abs(defaultsize), abs(textsize)
        if init: # True during initializing of filtus
            defaultsize -= 1; textsize -= 1
        cpsgn = copysign
        if init or (defaultsize != abs(olddef)):
            self.defaultfont['size'] = int(cpsgn(defaultsize, olddef))
            self.titlefont['size'] = int(cpsgn(defaultsize+5, olddef))
            self.smallfont['size'] = int(cpsgn(defaultsize-1, olddef))
            self.smallbold['size'] = int(cpsgn(defaultsize-1, olddef))
            self.tinyfont['size'] = int(cpsgn(defaultsize-3, olddef))
        if init or (textsize != abs(oldtext)):
            self.textfont['size'] = int(cpsgn(textsize, oldtext))
            self.monofont['size'] = int(cpsgn(textsize - 1, oldtext))
            self.monobold['size'] = int(cpsgn(textsize - 1, oldtext))
        if not init:
            self.alignAll()
            self.scrollframe.reposition()

    def alignAll(self):
        for w_name in ['gs', 'fs', 'vs']:
            getattr(self, w_name).align()
        self.sharingNotebook.setnaturalsize()
        self.FM.align()
        for w_name in ['advLoad', 'plink', 'fileReader']: #'resamplingDialog',
            if hasattr(self, w_name):
                getattr(self, w_name).align()
        Pmw.alignlabels([self.sepOutputOM, self.generalfontSizeEntry, self.textfontSizeEntry]) # not in separate class...TODO.

    def clearAll(self):
        for w in [self.fileListbox, self.fileSummary1, self.fileSummary2, self.text, self.gs, self.vs, self.fs, self.dn]:
            w.clearAll()
        self.files, self.filteredFiles, self.longFileNameList, self.shortFileNameList, self.currentFileNameList = [], [], [], [], []
        self.fileListbox.settoptext('Loaded samples: 0')
        gc.collect()

    def loadFiles(self, longFileNameList, **kwargs):
        ''' Load & summarize variant files.'''
        if len(longFileNameList) == 0: return
        longFileNameList = sorted(longFileNameList)

        if not hasattr(self, 'fileReader'):
            self.fileReader = InputDialog.InputDialog(filtus=self)

        fileReader = self.fileReader
        fileReader.prompt, fileReader.guess = True, True
        for filename in longFileNameList:
            VF = fileReader.read(filename, **kwargs)
            
            if fileReader.stopLoading: break
            if fileReader.skipFile: continue
            
            try:
                self.addFile(VF)
            except AttributeError: # for vcf files there can be more than one
                for VFi in VF: self.addFile(VFi)

        if len(self.files) == 0:
            self.fileListbox.settoptext('Loaded samples: 0')
        self.filteredFiles = [] #self.filteredFiles_initialcopy()
        self.fileSummary1.setTotal()

        for w in [self.fileSummary2, self.text, self.gs, self.vs, self.fs]:
            w.clearAll()
        self.setAllColnames()
        
    def addFile(self, VF):
        i = len(self.files) + 1
        currentName = VF.shortName if self.shortFilenames else VF.longName
        self.longFileNameList.append(VF.longName)
        self.shortFileNameList.append(VF.shortName)
        self.currentFileNameList.append(currentName)
        self.files.append(VF)
        self.fileListbox.settoptext('Loaded samples: %d' %i)
        self.fileListbox.newline('%2d: %s' %(i, currentName))
        self.fileSummary1.addVF(VF)
    
    def toggleShortNames(self, choice=None): 
        '''argument 'choice' should be either "short", "long" or None (toggle)''' 
        if choice is None:
            self.shortFilenames = not self.shortFilenames
        else:
            self.shortFilenames = choice=="short"
            
        self.shortFileNameList = [VF.shortName for VF in self.files] if self.files else  []
        self.currentFileNameList = self.shortFileNameList if self.shortFilenames else self.longFileNameList
        self.fileListbox.setAll(self.currentFileNameList)
        
    def unloadFiles(self, seleci=[]):
        if len(seleci) == 0:
            seleci = [int(i) for i in self.fileListbox.curselection()]
        if len(seleci) == 0:
            return
        N = len(self.longFileNameList)
        self.files[:] = [self.files[i] for i in range(N) if not i in seleci]
        self.longFileNameList[:] = [self.longFileNameList[i] for i in range(N) if not i in seleci]
        self.toggleShortNames(choice=self.shortFilenames) # also fixes shortnamelist

        self.fileSummary1.setAll(self.files)
        if self.fileSummary2.size() > 0:
            self.filteredFiles[:] = [self.filteredFiles[i] for i in range(N) if not i in seleci]
            self.fileSummary2.setAll(self.filteredFiles)
        self.setAllColnames()
        self.text.clearAll()

    def filteredFiles_initialcopy(self):  #ad hoc deep copy of self.files
        return [VF.copyAttributes() for VF in self.files]

    def setAllColnames(self):
        uniqueCols = FiltusUtils.listUnique([head for VF in self.files for head in VF.columnNames])
        self.FM.setColnames(uniqueCols)

        menubar = self.menuBar
        viewmenu = menubar.component('View-menu')
        if not uniqueCols:
            viewmenu.entryconfigure(0, state='disabled')
            return

        viewmenu.entryconfigure(0, state='normal')
        summarymenu = menubar.component('columnsum-menu')
        L = summarymenu.index('end')
        if L is not None:
            menubar.deletemenuitems('columnsum', 0, L)

        summarizer = FiltusAnalysis.ColumnSummary()
        for col in uniqueCols:
            menubar.addmenuitem('columnsum', 'command', None, label=col, font=self.defaultfont,
                command=self._showSummary(summarizer, col))
        
        nUnique = len(uniqueCols)
        nCols = -(-nUnique/28) # What the hell... hvorfor? heltallsdivisjon?
        N = -(-nUnique/nCols)
        for i in range(1, nCols):
            summarymenu.entryconfigure(i * N, columnbreak=1)

    def _showSummary(self, summarizer, col): 
        def _f():
            files = self.checkLoadedSamples(select="all", VF=True, filtered=True)
            s = summarizer.summarize(files, col)
            self.text.prettyPrint(s)
        return _f
        
class BusyManager(object):

    def __init__(self, widget):
        self.toplevel = widget.winfo_toplevel()
        self.widgets = {}
        self.isBusy = False
        
    def busy(self, widget=None):
        # attach busy cursor to toplevel, plus all windows
        # that define their own cursor.
        if widget is None:
            #print 'setting busy (prior state = %s)'%self.isBusy
            if self.isBusy: return
            w = self.toplevel # myself
            self.isBusy = True
        else:
            w = widget

        if not self.widgets.has_key(str(w)):
            try:
                # attach cursor to this widget
                cursor = w.cget("cursor")
                if cursor != "watch":
                    self.widgets[str(w)] = (w, cursor)
                    w.config(cursor="watch")
            except TclError:
                pass

        for w in w.children.values():
            self.busy(w)

        
    def notbusy(self):
        # restore cursors
        if not self.isBusy: return
        for w, cursor in self.widgets.values():
            try:
                w.config(cursor=cursor)
            except TclError:
                pass
        self.widgets = {}
        self.isBusy = False

v = sys.version_info
if v[:2] != (2,7):
    print "Python 2.7 is needed to run Filtus. Your Python version is %d.%d.%d"%v[:3]
    sys.exit(0)

root = Tk()
Pmw.initialise(root)
filtus = Filtus(root, version="0.99-94")
root.title("FILTUS " + filtus.version)

if __name__ == "__main__":
    
    printVF = FiltusAnalysis._printVF
    test_csv = ["C:\\Projects\\FILTUS\\Testfiles\\test%d.csv" %i for i in (1,2)]
    test_vcf = "C:\\Projects\\FILTUS\\Testfiles\\vcf_example.vcf"
    triotest = "C:\\Projects\\FILTUS\\Testfiles\\triotest.txt" 
    
    def test_version():
        '''checks that the correct version number is used when saving output files'''
        import inspect
        a = inspect.getargspec(FiltusUtils.preambleNY) # version number is the last argument of this function
        save_version = a.defaults[-1]
        assert save_version == filtus.version
        
    def test_loading():
        filtus.clearAll()
        filtus.loadFiles([test_csv[0]], guess=1, prompt=0, geneCol="")
        assert len(filtus.files) == 1
        filtus.toggleShortNames()
        assert len(filtus.files) == 1
        filtus.clearAll()
        filtus.loadFiles(test_csv, guess=1, prompt=0, geneCol="")
        assert len(filtus.files) == 2
        filtus.toggleShortNames()
        assert len(filtus.files) == 2
        filtus.clearAll()
        filtus.loadFiles([test_vcf], guess=1, prompt=0, geneCol="")
        assert len(filtus.files) == 3
        filtus.toggleShortNames()
        assert len(filtus.files) == 3
        
    def test_sharing():
        filtus.clearAll()
        filtus.loadFiles(test_csv, guess=1, prompt=0, geneCol="Gene", splitAsInfo="INFO")
        gs, vs = filtus.gs, filtus.vs
        gs.cases.setvalue('1,2')
        gs.button.invoke()
        vs.fields[4].setvalue('1,2')
        vs.button.invoke()
    
    def test_denovo():
        filtus.clearAll()
        filtus.loadFiles([triotest], guess=1, prompt=0, splitAsInfo="", formatCol="vcf_format", keep00=1, geneCol='Gene.refGene')
        dn = filtus.dn
        dn.focus()
        dn.child.setvalue('ID CH')
        dn.father.setvalue('ID FA')
        dn.mother.setvalue('ID MO')
        dn.freqMenu.setvalue('1000g2012apr_all')
        dn.freqEntry.setvalue('0.1')
        dn.button.invoke()
        dn2 = FiltusAnalysis.DeNovoComputer()
        res1 = filtus.text.currentColDat
        ch, fa, mo = filtus.files
        res2 = dn2.analyze(VFch=ch, VFfa=fa, VFmo=mo, trioID=[0,1,2], mut=1e-8, defaultFreq=0.1, altFreqCol='1000g2012apr_all')
        assert res1.length == res2.length
        assert all(a[0]==b[0] for a,b in zip(res1.variants, res2.variants))
    
    def test_autex():
        filtus.clearAll()
        filtus.loadFiles([triotest], guess=1, prompt=0, splitAsInfo="", formatCol="vcf_format", keep00=1, geneCol='Gene.refGene')
        filtus.fileListbox.selection_set(0)
        filtus.autexgui = FiltusWidgets.AutEx_GUI(filtus)
        filtus.autexgui._prepare()
        altfreq, deffreq, thresh, length, unit, n = '1000g2012apr_all', 0.05, 0.1, 0.1, 'Mb', 20
        for w, val in zip(('_altFreqMenu', '_def_freq_entry', '_thresh_entry', '_minlength_entry', '_unitmenu', '_mincount_entry'), (altfreq, deffreq, thresh, length, unit, n)):
            getattr(filtus.autexgui, w).setvalue(str(val))
        
        filtus.autexgui.execute("Compute")
        res1 = filtus.text.currentColDat
        
        autex = FiltusAnalysis.AutExComputer(genmapfile="C:\\Projects\\FILTUS\\DecodeMap_thin.txt")
        res2 = autex.autex_segments(filtus.files[0], f=0.01, a=.5, error=0.005, defaultFreq=deffreq, altFreqCol=altfreq, threshold=thresh, minlength=length, unit=unit, mincount=n)
        
        assert res1.length == res2.length
        assert all(a[1]==b[1] for a,b in zip(res1.variants, res2.variants))
        filtus.autozyg_prompt()
        
    def test_qc():
        filtus.clearAll()
        filtus.loadFiles([triotest], guess=1, prompt=0, splitAsInfo="", formatCol="vcf_format", keep00=1, geneCol='Gene.refGene')
        QC = FiltusQC.QC(filtus)
        QC._prepare()
        QC.save_browser.setvalue("C:\\Projects\\FILTUS\\kast.txt")
        QC.save_browser.select()
        QC.scatter_x.setvalue("DP")
        QC.scatter_y.setvalue("vcf_qual")
        QC.histo_var.setvalue("caddgt10")
        QC._executeDialogButton("Ok")
        QC._scatterButtonExecute()
        QC._histogramButtonExecute()
        filtus.QC = QC
        filtus.QC_prompt()
        
    def test_summary():
        filtus.clearAll()
        filtus.loadFiles(test_csv, guess=1, prompt=0, splitAsInfo="INFO")
        summarymenu = filtus.menuBar.component('columnsum-menu')
        summarymenu.invoke(0)
    
    def test_geneloopkup(): #TODO
        filtus.clearAll()
        filtus.loadFiles(test_csv[:1], guess=1, prompt=0, geneCol="Gene", splitAsInfo="")
        variants = FiltusAnalysis.geneLookup(["KIAA1751", "_FOO_"], VFlist=filtus.files)
        variants.save("_kast.txt")
        assert variants.length == 4 == variants.collapse().length
        filtus.loadFiles(test_csv, guess=1, prompt=0, geneCol="Gene", splitAsInfo="", prefilter=("do not contain", "_foo_"))
        variants = FiltusAnalysis.geneLookup(["KIAA1751", "_FOO_"], VFlist=filtus.files)
        assert variants.length == 8 == 2*variants.collapse().length
        
    def test_advancedload():
        filtus.clearAll()
        filtus.advLoad = FiltusWidgets.AdvancedLoad(filtus)
        filtus.advLoad.dirEntry.setvalue("C:\\Projects\\FILTUS\\Testfiles\\")
        filtus.advLoad.endswithEntry.setvalue(".vcf")
        filtus.advLoad.endswithVar.set(1)
        filtus.advLoad_prompt()
    
    def test_plink():
        filtus.clearAll()
        filtus.loadFiles([triotest], guess=1, prompt=0, keep00=1, splitAsInfo="")
        filtus.fileListbox.selection_set(0)
        filtus.plinkgui = FiltusWidgets.PLINK_GUI(filtus)
        filtus.plinkgui._prepare()
        filtus.plinkgui.execute("Compute")
        filtus.plink_prompt()
    
    if len(sys.argv) > 1 and sys.argv[1]=="test":
        test_version()
        test_loading()
        test_advancedload()
        test_summary()
        test_geneloopkup()
        test_sharing()
        test_qc()
        test_denovo()
        test_plink()
        test_autex()
        print 'all tests passed'
        
root.mainloop()

