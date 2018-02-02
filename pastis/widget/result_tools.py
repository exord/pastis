from Tkinter import *
import numpy as n

"""
def load_result_window(self):
	import os, pickle
	import tkFileDialog
	self.resultdir = tkFileDialog.askdirectory(initialdir = '/data/PASTIS/resultfiles/', title = "Select a PASTIS resultfiles directory")
	if self.resultdir <> '' :
		self.visu_result()
"""
def load_result_window(self):
	import tkFileDialog
	import tkMessageBox
	mcmcfilenames = tkFileDialog.askopenfilenames(initialdir = '/data/PASTIS/resultfiles/', defaultextension = '.mcmc', filetypes=[('PASTIS MCMC', '*.mcmc')], title = "Load PASTIS mcmc result files")
	if mcmcfilenames <> '' : 
		self.resultname = []
		self.resultrunid = []
		self.resultbeta = []
		for f in mcmcfilenames : 
			self.resultname.append(f.split('resultfiles/')[1].split('/')[0])
			self.resultrunid.append(f.split(self.resultname[-1])[2].split('Beta')[0][1:-1])
			self.resultbeta.append(float(f.split('Beta')[1].split('_rep')[0]))
		if len(set(self.resultname)) > 1: tkMessageBox.showwarning(title='Multiple targets', message='You should not import mcmc files from different targets')
		else : self.resultname2 = self.resultname[0]
		if len(set(self.resultrunid)) > 1: tkMessageBox.showwarning(title='Multiple runID', message='Be careful, you are going to load mcmc files with different runID/comments')
		else : self.resultrunid2 = self.resultrunid[0]
		if len(set(self.resultbeta)) > 1: tkMessageBox.showwarning(title='Multiple Betas', message='Be careful, you are going to load mcmc files with different beta values')
		else : self.resultbeta2 = self.resultbeta[0]
		self.chaindict = self.read_mcmc_chains(mcmcfilenames, self.resultname, self.resultrunid, self.resultbeta)
		self.update_list_jumpparams()
		self.visu_result()


def visu_result(self):
	import numpy as n
	self.w2 = Toplevel()

	self.rlabl = Label(self.w2, image = self.photo) # display the logo
	self.rlabl.pack()
		
	self.rlab = Label(self.w2, text="PASTIS result analysis tool", font=('Times', 14, 'bold')) # set a subtitle
	self.rlab.pack()

	Label(self.w2, text="PASTIS results for object:", font=('Times', 12)).pack() # set a subtitle
	self.rvname = StringVar()
	Label(self.w2, textvariable=self.rvname, font=('Times', 12)).pack() # set the target name
	self.rvname.set(self.resultname2)
	
	#self.resulttop = Frame(self.w2)

	self.resultpane = PanedWindow(self.w2, orient=HORIZONTAL)  # separate the main window in two sub-windows
	self.resultpane.pack(fill=BOTH, expand=1)

	self.resulttopleft = Frame(self.w2)  # define the left frame
	self.resulttopright = Frame(self.w2)  # define the left frame
	
	Label(self.resulttopleft, text="List of chains files", font=('Times', 12, 'bold')).pack()
	self.resultbox=Listbox(self.resulttopleft, bg='white', selectmode=EXTENDED, width = 70, height = 10)  # define the listbox
	self.resultbox.pack()
	self.update_mcmc_listbox()

	self.resultf1 = Frame(self.resulttopleft)
	#Button(self.resultf1, text='Load new chains', command=self.result_load_chains).pack(side='left')
	Button(self.resultf1, text='Remove chains', command=self.result_remove_chains).pack(side='left')
	Button(self.resultf1, text='Clear all chains', command=self.result_clean_chains).pack()
	self.resultf1.pack()
	
#	self.resulttopleft.pack()
	self.resultpane.add(self.resulttopleft)
	
	self.BIconverged = IntVar()
	self.BIconverged.set(90)
	self.result_sigma = DoubleVar()
	self.result_sigma.set(3.)
	self.BIlate = IntVar()
	self.BIlate.set(50)
	
	self.resulttoprightFrame1 = Frame(self.resulttopright)
	Scale(self.resulttoprightFrame1, variable=self.BIconverged, from_=0, to=100, length = 200, tickinterval = 30, showvalue = 'yes', orient='horizontal', label='Maximum <logL> BI').pack()
	Label(self.resulttoprightFrame1, text='Sigma:').pack(side='left')
	Entry(self.resulttoprightFrame1, textvariable=self.result_sigma, width=3, bg='white').pack(side='left')
	Button(self.resulttoprightFrame1, text='Remove outlier chains', command=self.remove_not_converged_chains).pack()
	self.resulttoprightFrame1.pack()
	
	self.resulttoprightFrame2 = Frame(self.resulttopright)
	Scale(self.resulttoprightFrame2, variable=self.BIlate, from_=0, to=100, length = 200, tickinterval = 30, showvalue = 'yes', orient='horizontal', label='Maximum convergence BI').pack()
	Button(self.resulttoprightFrame2, text='Remove late-convergence chains', command=self.remove_late_converged_chains).pack()
	self.resulttoprightFrame2.pack()
	self.resultpane.add(self.resulttopright)
	
	#self.resulttop.pack()
	############################
	self.resultbottom = Frame(self.w2)  # define the bottom frame

	Label(self.resultbottom, text='1. Chains visualisations tools', font=('Times', 12, 'bold')).grid(row=0, columnspan=7)
	
	self.resultvariable = StringVar()
	self.resultvariable.set('logL')
	self.BI = IntVar()
	self.BI.set(20)
	self.sampling = IntVar()
	self.sampling.set(100)
	
	Label(self.resultbottom, text='Parameter:').grid(row=1, column=0)
	OptionMenu(self.resultbottom, self.resultvariable, *self.list_jumpparams).grid(row=1, column=1)
	Label(self.resultbottom, text='Sampling:').grid(row=1, column=2)
	Entry(self.resultbottom, textvariable = self.sampling, width=10, bg='white').grid(row=1, column=3)
	Button(self.resultbottom, text = 'Plot chains !', command=self.plot_multichains).grid(row=1, column=6)
	Label(self.resultbottom, text='Burn-in phase (%):').grid(row=2, column=0)
	Scale(self.resultbottom, variable=self.BI, from_=0, to=100, length=300, tickinterval=25, showvalue='yes', orient='horizontal').grid(row=2, column=1, columnspan=5)
	
	Label(self.resultbottom, text='2. Compute correlation lenght', font=('Times', 12, 'bold')).grid(row=3, columnspan=7)
	
	self.corrstep = IntVar()
	self.corrstep.set(100)
	Label(self.resultbottom, text='Correlation step:').grid(row=4, column=0)
	Entry(self.resultbottom, textvariable=self.corrstep, width=10, bg='white').grid(row=4, column=1)
	Button(self.resultbottom, text = 'Parameter CorrLenght', command=self.compute_ParamCorrLenght).grid(row=4, column=2)
	
	self.status_corrlengt = StringVar()
	self.status_corrlengt.set('')
	Label(self.resultbottom, textvariable=self.status_corrlengt).grid(row=4, column=3)

	self.buttonChainCL = Button(self.resultbottom, text = 'Chain CorrLenght', command=self.compute_ChainCorrLenght, state = 'disabled')
	self.buttonChainCL.grid(row=4, column=4)
	
	self.status_corrlengt = StringVar()
	self.status_corrlengt.set('')
	Label(self.resultbottom, textvariable=self.status_corrlengt).grid(row=4, column=5)


	Label(self.resultbottom, text='3. Merge chains', font=('Times', 12, 'bold')).grid(row=5, columnspan=7)
	
	self.buttonMerged = Button(self.resultbottom, text = 'Merge chains', command=self.merge_chains, state = 'disabled', width = 30)
	self.buttonMerged.grid(row=6, column=0, columnspan=3)
	self.buttonMergedSave = Button(self.resultbottom, text = 'Save Merged chain', command=self.save_merged_chain, state = 'disabled', width = 30)
	self.buttonMergedSave.grid(row=6, column=4, columnspan=3)


	self.result_q = DoubleVar()
	self.result_q.set(0.6827)
	self.result_nb_bin = IntVar()
	self.result_nb_bin.set(50)
	self.result_difftol = DoubleVar()
	self.result_difftol.set(0.3)
	Label(self.resultbottom, text='4. Get confidence values', font=('Times', 12, 'bold')).grid(row=7, columnspan=7)
	
	Label(self.resultbottom, text='confidence limits').grid(row=8, column=0)
	Entry(self.resultbottom, textvariable=self.result_q, width=8, bg='white').grid(row=8, column=1)
	Label(self.resultbottom, text='nb bins').grid(row=8, column=2)
	Entry(self.resultbottom, textvariable=self.result_nb_bin, width=8, bg='white').grid(row=8, column=3)
	Label(self.resultbottom, text='diff tol').grid(row=8, column=4)
	Entry(self.resultbottom, textvariable=self.result_difftol, width=8, bg='white').grid(row=8, column=5)
	self.confidencebutton = Button(self.resultbottom, text='Confidence interval', command=self.get_confidence, state='disabled')
	self.confidencebutton.grid(row=8, column=6)
	
	Label(self.resultbottom, text='5. Export results', font=('Times', 12, 'bold')).grid(row=9, columnspan=7)
	self.pyramidbutton = Button(self.resultbottom, text='Make pyramid plot', command=self.pyramid_plot, state='disabled')
	self.pyramidbutton.grid(row=10, column=0)
	self.solutionbutton = Button(self.resultbottom, text='Make solution file', command=self.make_solution_file, state='disabled')
	self.solutionbutton.grid(row=10, column=2)
	self.plotbestbutton = Button(self.resultbottom, text='Plot best model', command=self.plot_best, state='active')
	self.plotbestbutton.grid(row=10, column=3)
	self.MCbutton = Button(self.resultbottom, text='Compare models', command=self.compare_model, state='active')
	self.MCbutton.grid(row=10, column=4)
	Button(self.resultbottom, text='Quit', command=self.w2.destroy).grid(row=10, column=6)

	self.resultbottom.pack()


def read_mcmc_chains(self, filenames, names, runids, betas):
#	import pickle
#	for f in filenames :
#		fname = f.split(self.resultname+'/')[1]
#		p = open(f, 'r')
#		vd[fname] = pickle.load(p)
	from ..MCMC import analysis
	#fnames, vdict = analysis.read_chains(self.resultname, self.resultname)
	fnames, vdict = analysis.read_chains2(self, filenames, betas, names, runids)
	vd = {}
	for i in xrange(len(fnames)) : 
		keydict = str(i)+'. '+names[i]+' runID='+runids[i]+' Beta='+str(betas[i])+' rep='+str(int(fnames[i].split('rep')[1].split('_')[0]))+' queued:'+fnames[i].split('queued')[1].split('_')[0]
		#vd[keydict] = {}
		vd[keydict] = vdict[i].get_value_dict()
	return vd

def update_list_jumpparams(self):
	import numpy as n
	from ..MCMC import analysis
	self.list_jumpparams = ['ALL']
	for f in self.chaindict.keys() :
		for key in self.chaindict[f].keys() :
			if key not in self.list_jumpparams :
				self.list_jumpparams.append(key)
	self.list_jumpparams = n.sort(self.list_jumpparams)

def update_mcmc_listbox(self):
	import numpy as n
	self.resultbox.delete(0, END)
	for f in n.sort(self.chaindict.keys()) :
		self.resultbox.insert(END, f)
	

def result_remove_chains(self):
	import tkMessageBox
	removeresult = []
	message = ''
	if len(self.resultbox.curselection()) == 0 :
		tkMessageBox.showerror(parent = self.w2, title = 'No chain selected', message = 'Please select a chain first !')
	else : 
		for i in self.resultbox.curselection():
			removeresult.append(self.resultbox.get(i))
			message+=removeresult[-1]+'\n'
		question = tkMessageBox.askyesno(master = self.w2, title = 'Are you sure ?', message = 'Are you sure to remove the following chains ?\n'+message)
		if question :
			for r in removeresult : 
				del self.chaindict[r]
		self.update_mcmc_listbox()

def result_remove_chains2(self):
	import tkMessageBox
	question = tkMessageBox.askyesno(parent = self.wplot, title = 'Are you sure ?', message = 'Are you sure to remove this chain ?')
	if question :
		r = n.sort(self.chaindict.keys())[self.highlightchain.get()]
		del self.chaindict[r]
	self.update_mcmc_listbox()
	self.wplot.destroy()
	self.plot_multichains()

def result_remove_chains3(self):
	import tkMessageBox
	nchains = len(self.chains_to_remove)
	question = tkMessageBox.askyesno(parent = self.w2, title = 'Are you sure ?', message = 'Are you sure to remove these ('+str(nchains)+')chains ?')
	if question :
		for i,r in enumerate(n.sort(self.chaindict.keys())):
			if i in self.chains_to_remove :
				del self.chaindict[r]
	self.update_mcmc_listbox()


def result_clean_chains(self):
	question = tkMessageBox.askyesno(parent = self.w2, title = 'Are you sure ?', message = 'Are you sure to remove all chains ?')
	if question : 
		for r in self.chaindict.keys() :
			del self.chaindict[r]
	self.update_mcmc_listbox()


def plot_multichains(self):
	from .. import plots
	from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
	from matplotlib.figure import Figure
	self.wplot = Toplevel()
	self.f = Figure(figsize=(7,5), dpi=100)
	self.canvas = FigureCanvasTkAgg(self.f, master=self.wplot)
	self.highlightchain = IntVar()
	self.highlightchain.set(0)
	self.get_plot(self.BI.get())
	self.canvas.get_tk_widget().pack(expand=1)
	toolbar = NavigationToolbar2TkAgg(self.canvas, self.wplot)
	toolbar.update()
	self.canvas._tkcanvas.pack(expand=1)
	#Label(self.resultbottom, text='Burn-in phase (%):').pack(side='left')
	Scale(self.wplot, variable=self.BI, from_=0, to=100, length=300, tickinterval=25, showvalue='yes', orient='horizontal', command=self.get_plot, label='Burn-in phase [%]').pack()
	#Label(self.resultbottom, text='Highlighted chain:').pack(side='left')
	Scale(self.wplot, variable=self.highlightchain, from_=0, to=len(self.chaindict.keys())-1, length=300, tickinterval=1, showvalue='yes', orient='horizontal', command=self.get_highlight, label='Highlighted chain').pack()
	Button(self.wplot, text='Remove chain', command=self.result_remove_chains2).pack()

def get_plot(self, BI):
	import numpy as n
	from .. import plots
	self.f.clf()
	
	if self.resultvariable.get() <> 'ALL' : 
		self.a = self.f.add_subplot(111)
		self.startindex = round(float(BI)/100.*len(self.chaindict[self.chaindict.keys()[0]][self.resultvariable.get()]))
		for k in n.sort(self.chaindict.keys()):
			self.a.plot(self.chaindict[k][self.resultvariable.get()][self.startindex::self.sampling.get()], alpha=0.1)
		self.get_highlight(self.highlightchain.get())
		self.a.set_xlabel('iterations %% %i' %self.sampling.get())
		self.a.set_ylabel(self.resultvariable.get())
	else : 
		for i, param in enumerate(self.list_jumpparams) :
			if param == 'ALL' : continue
			self.startindex = round(float(BI)/100.*len(self.chaindict[self.chaindict.keys()[0]][param]))
			self.a = self.f.add_subplot(len(self.list_jumpparams), 1, i)
			for k in n.sort(self.chaindict.keys()):
				self.a.plot(self.chaindict[k][param][self.startindex::self.sampling.get()], alpha=0.1)
			self.get_highlight(self.highlightchain.get())

def get_highlight(self, highlightchain):
	for l in xrange(len(self.a.lines)):
		self.a.lines[l].set_alpha(0.1)
	self.a.lines[int(highlightchain)].set_alpha(1.)
	self.canvas.show()

def compute_ParamCorrLenght(self):
	from ..MCMC import analysis
	vd = []
	for f in self.chaindict.keys(): vd.append(self.chaindict[f])
	try : self.ParamCorrLengt = analysis.corrlength_multichain2(vd, step=self.corrstep.get(), BI=self.maxBI, plot=False, widget = True)
	except : self.ParamCorrLengt = analysis.corrlength_multichain2(vd, step=self.corrstep.get(), BI=float(self.BI.get())/100., plot=False, widget = True)
	self.buttonChainCL.config(state = 'active')

def compute_ChainCorrLenght(self):
	from ..MCMC import analysis
	self.ChainCorrLengt = analysis.corrlenchain(self.ParamCorrLengt)
	self.buttonMerged.config(state = 'active')

def merge_chains(self):
	from ..MCMC import analysis
	vd = []
	for f in self.chaindict.keys(): vd.append(self.chaindict[f])
	try : self.mergedchain = analysis.merge_chains(vd, self.maxBI, self.ChainCorrLengt, N = 1e4, pickrandom = False, beta = None)
	except : self.mergedchain = analysis.merge_chains(vd, float(self.BI.get())/100., self.ChainCorrLengt, N = 1e4, pickrandom = False, beta = None)
	self.buttonMergedSave.config(state = 'active')
	self.solutionbutton.config(state = 'active')
	self.confidencebutton.config(state = 'active')
	self.pyramidbutton.config(state = 'active')

def save_merged_chain(self):
	from ..MCMC import analysis
	import tkMessageBox
	analysis.register_chain(self.mergedchain, beta = self.resultbeta2, target = self.resultname2, runid = self.resultrunid2)
	savedfile = os.path.join(resultpath,self.resultname2, '%s_%s_Beta%.6f_mergedchain.dat'%(self.resultname2, self.resultrunid2, self.resultbeta2) )
	tkMessageBox.showinfo(parent = self.w2, title='Merged chain saved', message='Merged chain saved in '+savedfile)

def get_confidence(self):
	from ..MCMC import analysis
	analysis.confidence_intervals(self.mergedchain, q = self.result_q.get(), nbins = self.result_nb_bin.get(), burnin=float(self.BI.get()) / 100., difftol = self.result_difftol.get())

def make_solution_file(self):
	from ..MCMC import analysis
	import tkFileDialog
	solutionfile = tkFileDialog.askopenfilename(initialdir = '/data/PASTIS/configfiles/', defaultextension = '.pastis', filetypes=[('PASTIS Configfile', '*.pastis')], title = "Select the PASTIS config file")
	if solutionfile <> '' : 
		analysis.make_solution_file(self.mergedchain, pastisfile = solutionfile, BI = 0., best = True)
	self.plotbestbutton.config(state = 'active')

def pyramid_plot(self):
	from .. import plots
	import tkFileDialog
	import tkMessageBox
	pyramidfile = tkFileDialog.asksaveasfilename(parent = self.w2, initialdir = '/data/PASTIS/resultfiles/', defaultextension = '.pdf', title = "Select the pyramid plot filename")
	if pyramidfile <> '' : 
		plots.pyramid_cont(self.mergedchain, BI = 0.0, Nlast = None, sample = 1, nbins = self.result_nb_bin.get(),label = self.resultname2, filename = pyramidfile)
	tkMessageBox.showinfo(parent = self.w2, title='Pyramid plot saved', message='Pyramid plot saved in '+pyramidfile)


def plot_best(self):
	import tkFileDialog
	return


def compare_model(self):
	import tkFileDialog
	return



def remove_not_converged_chains(self):
	"""
	need self.BIconverged, self.result_sigma
	"""
	import numpy as n
	n_chains = len(self.chaindict.keys())
	medLogL, stdLogL = n.zeros(n_chains, float), n.zeros(n_chains, float)
	self.startindex_converged = round(float(self.BIconverged.get())/100.*len(self.chaindict[self.chaindict.keys()[0]]['logL']))
	for i,k in enumerate(n.sort(self.chaindict.keys())):
		medLogL[i] = n.median(self.chaindict[k]['logL'][self.startindex_converged:])
		stdLogL[i] = n.std(self.chaindict[k]['logL'][self.startindex_converged:])
	maxLogL = max(medLogL)
	self.chains_to_remove = n.where(medLogL < maxLogL - self.result_sigma.get() * n.sqrt(stdLogL[n.argmax(medLogL)]**2. + stdLogL**2))[0]
	if len(self.chains_to_remove) > 0 : self.result_remove_chains3()
	else : print 'no chain to remove'


def remove_late_converged_chains(self):
	"""
	need self.BIlate
	"""
	import numpy as n
	from ..MCMC import analysis
	BI = []
	vd = []
	for f in n.sort(self.chaindict.keys()): vd.append(self.chaindict[f])
	
	for p in vd[0].keys():
		if p in ['logL', 'posterior'] : continue
		z, zz, BIfrac = analysis.find_BI(vd, param=p)
		BI.append(BIfrac)
	BI = n.array(BI)
	n_chains = len(vd)
	maxBI = n.zeros(n_chains, float)
	for i in xrange(n_chains):
		maxBI[i] = max(BI[:,i])
	self.chains_to_remove = n.where(maxBI > self.BIlate.get()/100.)[0]
	if len(self.chains_to_remove) > 0 : self.result_remove_chains3()
	else : print 'no chain to remove'
	self.maxBI = maxBI[n.where(maxBI <= self.BIlate.get()/100.)[0]]
