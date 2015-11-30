import Tkinter
import Pmw
import os.path
import math
import time
import tkFileDialog
import random

#import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.collections as mplcol
import matplotlib.transforms as mpltransforms

import FiltusWidgets
import FiltusUtils
import FiltusDatabase


class QC(object):
    def __init__(self, filtus):
        self.parent = filtus.parent
        self.filtus = filtus
        self.createDialog()
        
    def createDialog(self):
        filtus = self.filtus
        self.dialog = Pmw.Dialog(self.parent, title='Quality plots', buttons=('Close',), activatecommand=self._prepare, command=self._executeDialogButton, dialogchildsite_pady=5, buttonbox_pady=10)
        self.dialog.withdraw()
        interior0 = self.dialog.interior()
        fr = Tkinter.Frame(interior0) #self.dialog.interior()
        fr.columnconfigure(0, weight=1)
        fr.rowconfigure(1, weight=1)
        fr.grid(row=0, column=0, pady=10, sticky='news')
        FiltusWidgets.HelpButton(interior0, filtus=filtus, page="qcplots").grid(row=0, column=0, sticky="ne")
        
        #fr = self.dialog.interior()
        Tkinter.Label(fr, text="QUALITY CONTROL PLOTS", font=filtus.titlefont).grid(sticky='news', pady=8)
        
        button_OPTIONS = dict(menubutton_anchor = 'w', menubutton_padx=5, menubutton_pady=1, menubutton_width=9, labelmargin=5, menu_font = filtus.defaultfont)
        entry_OPTIONS = dict(labelpos="w", entry_width=5, labelmargin=5, entry_justify='center')
        grid_OPTIONS = dict(sticky='news', pady=5, padx=10)
        
        names_group = Pmw.Group(fr, tag_text = 'Select samples')
        names_interior = names_group.interior()
        names_interior.rowconfigure(0, weight=1)
        names = FiltusWidgets.LabeledListBox(names_interior, filtus=filtus, toptext="", width=30, 
                                            scrolledlistbox_selectioncommand=self._updateSelStatus,
                                            bottomtext = "Selected: 0", height=min(6, len(filtus.files)))
        names.component('scrolledlistbox_label').destroy()
        names.component('scrolledlistbox_listbox').bind('<Control-a>', self._selectall_and_update)
        names.component('scrolledlistbox_listbox').bind('<KeyRelease-Up>', self._updateSelStatus)
        names.component('scrolledlistbox_listbox').bind('<KeyRelease-Down>', self._updateSelStatus)
        names.grid(sticky='news')
        names_interior.grid(padx=5, pady=5)
        self.names = names
        
        comparative_group = Pmw.Group(fr, tag_text = 'Comparative plots')
        comparative_interior = comparative_group.interior()
        self.comparative_checks = Pmw.RadioSelect(comparative_interior, buttontype = 'checkbutton', orient = 'horizontal')
        for name, txt in zip(['gender', 'private', 'heterozygosity'], [' Gender', ' Private variants', ' Heterozygosity']):
            px = 15 if name=="private" else 0
            self.comparative_checks.add(name, text=txt, padx=px)
            self.comparative_checks.invoke(name)
        
        self.save_browser = FiltusWidgets.FileBrowser(comparative_interior, filtus=filtus, label="Write to text:", 
                checkbutton = True, labelpos='w', browsesticky='se', entryfield_entry_width=15, browsetitle="")
        self.save_browser.browsebutton.configure(command = self._browseSave)
        self.save_browser.entryfield.configure(command = None)
        
        self.comparative_checks.grid(row=0, column=0, **grid_OPTIONS)
        self.save_browser.grid(row=1, column=0, **grid_OPTIONS)
        Tkinter.Button(comparative_interior, text="Plot!", command=self._comparativeButtonExecute, bd=3, padx=5, pady=3).grid(row=0, column=1, rowspan=2, padx=30, pady=10)
        
        ### scatter plots
        scatter_group = Pmw.Group(fr, tag_text = 'Scatter plots')
        scatter_interior = scatter_group.interior()
        self.scatter_x = FiltusWidgets.OptionMenuExt(scatter_interior, labelpos="w", label_text="X variabel:", **button_OPTIONS)
        self.scatter_y = FiltusWidgets.OptionMenuExt(scatter_interior, labelpos="w", label_text="Y variabel:", **button_OPTIONS)
        self.scatter_alpha = Pmw.EntryField(scatter_interior, label_text="Transp.:", value='0.05', validate = {'validator':'real','min':0, 'max':1,'minstrict':0, 'maxstrict':0}, **entry_OPTIONS)
        self.scatter_thin = Pmw.EntryField(scatter_interior, label_text="Thin by:", value='1', validate = {'validator':'integer','min':1,'minstrict':0}, **entry_OPTIONS)
        
        self.scatter_x.grid(**grid_OPTIONS)
        self.scatter_y.grid(**grid_OPTIONS)
        self.scatter_thin.grid(row=0, column=1, **grid_OPTIONS)
        self.scatter_alpha.grid(row=1, column=1, **grid_OPTIONS)
        Tkinter.Button(scatter_interior, text="Plot!", command=self._scatterButtonExecute, bd=3, padx=5, pady=3).grid(row=0,column=2, rowspan=2, padx=30, pady=10)
        
        ### histograms
        histo_group = Pmw.Group(fr, tag_text = 'Histogram plots')
        histo_interior = histo_group.interior()
        #histo_interior.columnconfigure(2, weight=1)
        
        self.histo_var = FiltusWidgets.OptionMenuExt(histo_interior, labelpos="w", label_text="Variabel:", **button_OPTIONS)
        self.histo_bins = Pmw.EntryField(histo_interior, label_text="Bins:", value='20', validate = {'validator':'integer','min':1,'minstrict':0}, **entry_OPTIONS)
        self.histo_var.grid(**grid_OPTIONS)
        self.histo_bins.grid(row=0, column=1, **grid_OPTIONS)
        Tkinter.Button(histo_interior, text="Plot!", command=self._histogramButtonExecute, bd=3, padx=5, pady=3).grid(row=0,column=2, padx=30, pady=10)
        
        for g in (names_group, comparative_group, scatter_group, histo_group):
            g.interior().columnconfigure(0, weight=1)        
            g.configure(tag_font = filtus.smallbold)
            g.grid(ipady=5, **grid_OPTIONS)

        Pmw.alignlabels([self.scatter_x, self.scatter_y, self.histo_var])
        Pmw.alignlabels([self.scatter_alpha, self.scatter_thin, self.histo_bins])
        
    def _prepare(self):
        files = self.filtus.files
        self.names.setlist(['%2d: %s' %(i + 1, os.path.basename(VF.shortName)) for i, VF in enumerate(files)])
        self._selectall_and_update()
        
        cols = FiltusUtils.listUnique([head for VF in files for head in VF.columnNames])
        for colmenu in [self.scatter_x, self.scatter_y, self.histo_var]:
            colmenu.setItems(['']+cols)
        
    
    def _browseSave(self):
        fil = tkFileDialog.asksaveasfilename(initialdir=self.filtus.currentDir, title = "Save plot data as")
        if fil:
            self.filtus.currentDir = os.path.dirname(fil)
            self.save_browser.setvalue(os.path.normpath(fil))
             
    def _updateSelStatus(self, event=None):
        self.names.setbottomtext('Selected: %d' %len(self.names.curselection()))
        
    def _selectall_and_update(self, event=None):
        self.names.selectall()
        self._updateSelStatus()
        
    def clearAll(self):
        self.save_browser.deselect()
        
    def _scatterButtonExecute(self):
        try:   
            xcol, ycol = self.scatter_x.getvalue(), self.scatter_y.getvalue()
            if xcol=='': raise RuntimeError("X axis column not selected")
            if ycol=='': raise RuntimeError("Y axis column not selected")
            VFlist = self._validateInput(checkPresence=[xcol, ycol])
            alpha = float(self.scatter_alpha.getvalue())
            thin = int(self.scatter_thin.getvalue())
            scatterPlot(VFlist, xcol, ycol, alpha, thin)
        except Exception as e:
            FiltusUtils.warningMessage(e)
            
    def _histogramButtonExecute(self):
        try:   
            col = self.histo_var.getvalue()
            if col=='': raise RuntimeError("Column variable not selected")
            VFlist = self._validateInput(checkPresence=[col])
            bins = int(self.histo_bins.getvalue())
            histogramPlot(VFlist, col, bins)
        except Exception as e:
            FiltusUtils.warningMessage(e)
            
    def _comparativeButtonExecute(self):
        try:
            VFlist = self._validateInput()
            plotselect = self.comparative_checks.getvalue()
            p, h, g = (str in plotselect for str in ['private','heterozygosity', 'gender'])
            writetofile = self.save_browser.getvalue() if self.save_browser.on() else None
            QC_3plots(VFlist, private=p, heterozygosity=h, gender=g, writetofile=writetofile, save=None)
        except Exception as e:
            FiltusUtils.warningMessage(e)
            
    def _executeDialogButton(self, button):
        if button is None or button == 'Close':
            self.dialog.deactivate()
            return
    
    def _validateInput(self, checkPresence=[]):
        VFlist = self.filtus.checkLoadedSamples(select="all", VF=1, filtered=1)
        VFlist = [VFlist[int(i)] for i in self.names.curselection()]
        if not VFlist: 
            raise RuntimeError("No samples selected")
        for col in checkPresence:
            if any(col not in VF.columnNames for VF in VFlist):
                raise ValueError("The column '%s' does not appear in all the selected samples" %col)
        return VFlist   
    
        
def setPlotParams(ax, title, xlab, ylab, xlim=None, ylim=None):
    ax.set_title(title, size="medium")
    ax.set_xlabel(xlab, size="smaller")
    ax.set_ylabel(ylab, size="smaller")
    ax.tick_params(axis='both', labelsize="smaller")
    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)
    orig_xt = ax.get_xticks()
    ax.margins(0.05, 0.05, tight=True)
    if min(orig_xt) == 0 and max(orig_xt) == 1:
        new_xt = orig_xt
    else:
        xt = ax.get_xticks()
        new_xt = sorted(set(math.floor(x) for x in xt))
        if len(new_xt) < 4: new_xt = xt # undo
        if len(new_xt) > 5:  new_xt = new_xt[1::2]
    ax.xaxis.set_ticks(new_xt)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_aspect(float((xlim[1]-xlim[0])/(ylim[1]-ylim[0]))) # float to avoid matplotlib (or numpy?) unicodewarning
    
def scatterPlot(VFlist, xcol, ycol, alpha, thin, NA_vals = ('', 'NA', '.', '-'), GTlegend="upper left", save=None, show=True):
    N = len(VFlist)
    nrow = int(math.sqrt(N))
    ncol = math.ceil(float(N)/nrow)
    fig = plt.figure(figsize=(3.5*ncol, 3.5*nrow))
     
    for i, VF in enumerate(VFlist):
        getData = VF.columnGetter(xcol, ycol)
        GTnum, keep00 = VF.GTnum(), VF.keep00
        data = [(getData(v), GTnum(v)) for v in VF.variants[::thin]]
        try:
            floatdata = [(float(x), float(y), gt) for (x,y), gt in data if x not in NA_vals and y not in NA_vals] 
        except ValueError:
            raise ValueError("Cannot plot columns with non-numerical values.\nTip: Use filters to remove non-numerical values before plotting.")
        
        ax = fig.add_subplot(nrow, ncol, i+1, aspect=1)
        for num, col in zip([0,1,2], ['b','g','r']):
            if num==0 and not keep00: continue
            gt_subset = [(x,y) for x,y,gt in floatdata if gt==num]
            if gt_subset:
                X, Y = zip(*gt_subset)
                ax.plot(X, Y, alpha=alpha, color=col, ls='none', marker='o')
       
        ## GT legend
        txt = ['ref/ref', 'ref/alt', 'alt/alt'][not keep00:]
        fmt = ['bo', 'go', 'ro'][not keep00:]
        axes = [ax.plot([], [], y)[0] for y in fmt]
        ax.legend(axes, txt, numpoints=1, fontsize="small", fancybox=True, handletextpad=0.2, loc=GTlegend)
        
        setPlotParams(ax, VF.shortName, xcol, ycol)
    showAndSave(fig, show=show, save=save)
    return fig   
 
def histogramPlot(VFlist, column, bins, NA_vals = ('', 'NA', '.', '-'), save=None, show=True):
    N = len(VFlist)
    nrow = int(math.sqrt(N))
    ncol = math.ceil(float(N)/nrow)
    fig = plt.figure(figsize=(3.5*ncol, 3.5*nrow))
    for i, VF in enumerate(VFlist):
        getData = VF.columnGetter(column)
        stringdat = [getData(v) for v in VF.variants]
        floatdat = [float(x) for x in stringdat if x not in NA_vals] 
        ax = fig.add_subplot(nrow, ncol, i+1, aspect=1)
        ax.hist(floatdat, bins=bins, color='b')
        setPlotParams(ax, VF.shortName, column, '')
    showAndSave(fig, show=show, save=save)
    return fig   
    
def homPlot(pos, obs, scores, freqs, title='', segs=None, save=None, show=True):
    fig = plt.figure(figsize=(14, 6))
    ax = fig.add_subplot(1, 1, 1)
    obs_jit = [0.4*gt + 0.18*fr + random.random()*0.02 for gt,fr in zip(obs, freqs)]
    ax.plot(pos, obs_jit, 'bo')
    ax.plot(pos, scores, 'r-')
    ax.add_collection(mplcol.BrokenBarHCollection([(s[0], s[2]*1.0e6) for s in segs], (0,1), facecolor='green', alpha=0.2))
    ax.set_title(title, size="larger")
    ax.set_xlabel("Chromosomal position (bp)", size="medium")
    ax.set_ylabel("Posterior probabilities", size="medium")
    ax.set_ylim([-0.05, 1.02])
    ax.get_yaxis().tick_left()
    trans = mpltransforms.blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(1.02, .9, '1/1', color="blue", size=12, verticalalignment='center', transform=trans)
    ax.text(1.02, .5, '0/1', color="blue", size=12, verticalalignment='center', transform=trans)
    ax.text(1.02, .1, '0/0', color="blue", size=12, verticalalignment='center', transform=trans)
    ax.margins(0.01, 0.15, tight=True)
    showAndSave(fig, tight=False, show=show, save=save)
    return fig   
    
def homPlotSimple(pos, obs, title='', segs=None, save=None, show=True): # plink plot, without scores and frequencies
    fig = plt.figure(figsize=(14, 6))
    ax = fig.add_subplot(1, 1, 1)
    obs_jit = [0.4*gt + random.random()/5 for gt in obs]
    ax.plot(pos, obs_jit, 'bo')
    ax.add_collection(mplcol.BrokenBarHCollection([(s[0], s[2]*1.0e6) for s in segs], (0,1), facecolor='green', alpha=0.2))
    ax.set_title(title, size="larger")
    ax.set_xlabel("Chromosomal position (bp)", size="medium")
    ax.set_ylim([-0.02, 1.02])
    ax.get_yaxis().set_ticks([])
    trans = mpltransforms.blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(1.02, .9, '1/1', color="blue", size=12, verticalalignment='center', transform=trans)
    ax.text(1.02, .5, '0/1', color="blue", size=12, verticalalignment='center', transform=trans)
    ax.text(1.02, .1, '0/0', color="blue", size=12, verticalalignment='center', transform=trans)
    ax.margins(0.01, 0.15, tight=True)
    showAndSave(fig, tight=False, show=show, save=save)
    return fig   
    
    
def QC_3plots(VFlist, gender=True, private=True, heterozygosity=True, writetofile=None, save=None, show=True):
    if private + heterozygosity + gender == 0: return None
    N = len(VFlist)
    add_legend = N < 13
    Nplots = private + heterozygosity + gender + add_legend
    nrow = int(math.sqrt(Nplots))
    ncol = math.ceil(float(Nplots)/nrow)
    fig = plt.figure(figsize=(3.5*ncol, 3.5*nrow))
    
    if add_legend:
        markers = ['D','^','*','d','<','s','p','v','D','^','*','d']
        sizes = [6,8,8,8,8,8,8,8,6,8,8,8]
        cols = ['red', 'lime', 'cyan', 'brown', 'magenta', 'gold', 'pink', 'black', 'purple', 'gray', 'silver', 'green']
    else:
        markers, sizes, cols = ['o']*N, [6]*N, ['red']*N
    
    DB = FiltusDatabase.VariantDatabase.buildFromSamples(VFlist, "Extended")
    db_str = DB.variants
    
    if writetofile:
        sep = '\t'
        text_out = FiltusUtils.preambleNY(VFlist, analysis="QC PLOTS")
        
    plotnr = 0
    if gender:
        plotnr += 1
        ax_sex = fig.add_subplot(nrow, ncol, plotnr, aspect=1)
        XminusPAR = FiltusUtils.XminusPAR
        db_X_raw = [x[6:] for x in db_str if XminusPAR(x[:2])]
        if db_X_raw:
            db_X = zip(*db_X_raw)
            totals_X = [sum(map(bool, x)) for x in db_X]
            hets = [sum(g == 1 for g in sample)*100.0/tot if tot>0 else 0 for sample, tot in zip(db_X, totals_X)]
            for i in range(N): 
                ax_sex.plot(totals_X[i], hets[i], marker=markers[i], color=cols[i], markersize=sizes[i])
        else:
            totals_X, hets = [0]*N, [0]*N
            #print "Empty gender estimation plot.\n\nNo variants found on X \ PAR."
        setPlotParams(ax_sex, "Gender estimation", 'Variants on X (-PAR)', 'Heterozygosity (%)', ylim=(0,100))
        ax_sex.axhspan(0, 15, facecolor='blue', alpha=0.2)
        ax_sex.axhspan(15, 35, facecolor='red', alpha=0.2)
        ax_sex.axhspan(35, 100, facecolor='green', alpha=0.2)
        
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax_sex.text(0.05, 0.95, "FEMALE", transform=ax_sex.transAxes, fontsize="x-small", va="top", ha='left', bbox=props)
        ax_sex.text(0.05, 0.27, "? ", transform=ax_sex.transAxes, fontsize="x-small", va="center", ha='left', bbox=props)
        ax_sex.text(0.95, 0.05, "MALE", transform=ax_sex.transAxes, fontsize="x-small", va="bottom", ha='right', bbox=props)
    
        if writetofile:
            headers = sep.join(['Sample', 'Variants on X (-PAR)', 'Heterozygosity (%)', 'Gender'])
            genders = ['?' if tot==0 or 15<h<35 else 'Male' if h<=15 else 'Female' for tot, h in zip(totals_X, hets)]
            points = [sep.join([s, str(x), '%.2f'%y, g]) for s,x,y,g in zip(DB.sampleNames, totals_X, hets, genders)]
            text_out += "***Plot: Gender estimation***\n" + headers + '\n' + '\n'.join(points) + '\n\n'
       
    if private:
        plotnr += 1
        ax_priv = fig.add_subplot(nrow, ncol, plotnr, aspect=1)
        db_nonz = [map(bool, x) for x in zip(*db_str)[6:]]
        totals_all = map(sum, db_nonz)
        
        if max(totals_all)>2000:
            totals_all = [tot/1000.0 for tot in totals_all]
            xlab = '# variants/1000'
        else: xlab = '# variants'
        rowSums_nonz = map(sum, zip(*db_nonz))
        
        priv_ind = [i for i in range(len(rowSums_nonz)) if rowSums_nonz[i]==1]
        privates = [sum(sampl[i] for i in priv_ind) for sampl in db_nonz]
        for i in range(N):
            ax_priv.plot(totals_all[i], privates[i], marker=markers[i], color=cols[i], markersize=sizes[i])
        setPlotParams(ax_priv, "Private variants", xlab, 'Private')
             
        if writetofile:
            headers = sep.join(['Sample', xlab, 'Private'])
            points = [sep.join([s, str(x), str(y)]) for s,x,y in zip(DB.sampleNames, totals_all, privates)]
            text_out += "***Plot: Private variants***\n" + headers + '\n' + '\n'.join(points) + '\n\n'
            
    if heterozygosity:
        plotnr += 1
        ax_het = fig.add_subplot(nrow, ncol, plotnr, aspect=1)
        chromInt = FiltusUtils.chromInt
        db_AUT = zip(*[x[6:] for x in db_str if chromInt(x[0]) < 23])
        if not db_AUT:
            raise RuntimeError("Empty heterozygosity plot.\n\nNo autosomal variants found.")  
        totals_AUT = [sum(map(bool, x)) for x in db_AUT]
        hets = [sum(g == 1 for g in sample)*100.0/tot if tot>0 else 0 for sample, tot in zip(db_AUT, totals_AUT)]
        if max(totals_AUT) > 2000:
            totals_AUT = [tot/1000.0 for tot in totals_AUT]
            xlab = '# autosomal variants/1000'
        else: xlab = '# autosomal variants'
        
        for i in range(N): 
            ax_het.plot(totals_AUT[i], hets[i], marker=markers[i], color=cols[i], markersize=sizes[i])
        setPlotParams(ax_het, "Heterozygosity", xlab, 'Heterozygosity (%)', ylim=(-5,105))
        
        if writetofile:
            headers = sep.join(['Sample', 'A'+xlab[3:], 'Heterozygosity (%)'])
            points = [sep.join([s, str(x), '%.2f'%y]) for s,x,y in zip(DB.sampleNames, totals_AUT, hets)]
            text_out += "***Plot: Heterozygosity***\n" + headers + '\n' + '\n'.join(points) + '\n'
     
    if writetofile:
        with open(writetofile, 'w') as out:
            out.write(text_out)
            
    if add_legend:
        plotnr +=1 
        ax_legend = fig.add_subplot(nrow, ncol, plotnr, aspect=1)
        simplenames = [VF.shortName for VF in VFlist]
        ax_legend.set_frame_on(False)
        ax_legend.axis('off')
        for i in range(N):
            ax_legend.plot([], marker=markers[i], color=cols[i], markersize=sizes[i], label=simplenames[i], ls='None')
        ax_legend.legend(loc=2, numpoints=1, fontsize='small', frameon=False, title="Legend")
        
    showAndSave(fig, tight=True, show=show, save=save)
    return fig
    
    
def showAndSave(fig, tight=True, show=True, save=None):
    try: 
        if tight: fig.set_tight_layout(True)
        if show: fig.show()
        if save: fig.savefig(save)
    except ValueError as e:
        raise ValueError("Could not show plot.\n\nInternal error message:%s"%e)