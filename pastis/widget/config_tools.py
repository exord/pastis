#from Tkinter import *
import smtplib
import tkFileDialog
import tkMessageBox

def save_prior(self):
	import numpy as n
	from scipy import stats 
	function = self.v6.get()
	#function = function.replace('x2', self.v2.get())
	#function = function.replace('x3', self.v3.get())
	#function = function.replace('x4', self.v4.get())
	#function = function.replace('x5', self.v5.get())
	#function = lambda x1, x2, x3, x4, x5 : eval(self.v6.get())
	#self.list_prior['prior'+str(self.nb_prior)] = { function : [self.v7.get(), self.v8.get(), self.v9.get(), self.v10.get(), self.v11.get(), self.v1.get(), self.v2.get(), self.v3.get(), self.v4.get(), self.v5.get()] }
	self.list_prior['prior'+str(self.nb_prior)] = { function : [self.v7.get(), self.v8.get(), self.v9.get(), self.v10.get(), self.v11.get(), self.v12.get(), self.v13.get(), self.v14.get()], 'variables' : [self.v1.get(), self.v2.get(), self.v3.get(), self.v4.get(), self.v5.get()] }
	self.nb_prior+=1
	self.refresh_list()
#	self.wprior.destroy()

def save_prior2(self):
	import numpy as n
	from scipy import stats 
	function = self.v6.get()
	#function = function.replace('x2', self.v2.get())
	#function = function.replace('x3', self.v3.get())
	#function = function.replace('x4', self.v4.get())
	#function = function.replace('x5', self.v5.get())
	#function = lambda x1, x2, x3, x4, x5 : eval(self.v6.get())
	#self.list_prior['prior'+str(self.nb_prior)] = { function : [self.v7.get(), self.v8.get(), self.v9.get(), self.v10.get(), self.v11.get(), self.v1.get(), self.v2.get(), self.v3.get(), self.v4.get(), self.v5.get()] }
	self.list_prior['prior'+str(self.nb_prior)] = { function : [self.v7.get(), self.v8.get(), self.v9.get(), self.v10.get(), self.v11.get(), self.v12.get(), self.v13.get(), self.v14.get()], 'variables' : [self.v1.get(), self.v2.get(), self.v3.get(), self.v4.get(), self.v5.get()] }
	self.refresh_list()
#	self.wprior.destroy()


def define_prior(self):
	self.wprior = Toplevel()
	Label(self.wprior, text='Define a custom prior',font=('Times', 16, 'bold')).grid(row=0, columnspan = 10)
	self.list_variables = []
	for o in self.list_objects.keys() :
		for p in self.list_objects[o].keys() : 
			self.list_variables.append(o+'_'+p)
	Label(self.wprior, text='x1').grid(row=1, column = 1)
	Label(self.wprior, text='x2').grid(row=2, column = 1)
	Label(self.wprior, text='x3').grid(row=3, column = 1)
	Label(self.wprior, text='x4').grid(row=4, column = 1)
	Label(self.wprior, text='x5').grid(row=5, column = 1)
	self.v1 = StringVar()
	self.v2 = StringVar()
	self.v3 = StringVar()
	self.v4 = StringVar()
	self.v5 = StringVar()
	OptionMenu(self.wprior, self.v1, *self.list_variables).grid(row=1, column=2, columnspan=8, sticky=W+E)
	OptionMenu(self.wprior, self.v2, *self.list_variables).grid(row=2, column=2, columnspan=8, sticky=W+E)
	OptionMenu(self.wprior, self.v3, *self.list_variables).grid(row=3, column=2, columnspan=8, sticky=W+E)
	OptionMenu(self.wprior, self.v4, *self.list_variables).grid(row=4, column=2, columnspan=8, sticky=W+E)
	OptionMenu(self.wprior, self.v5, *self.list_variables).grid(row=5, column=2, columnspan=8, sticky=W+E)

	Label(self.wprior, text='Function\n(x1, x2, x3, x4, x5)').grid(row=6, column = 1)
	self.v6 = StringVar()
	Entry(self.wprior, textvariable = self.v6, width=40, bg='white').grid(row=6, column=2, columnspan=8, sticky=W+E)
	
	self.v7 = DoubleVar()
	self.v8 = IntVar()
	self.v8.set(0)
	self.v9 = StringVar()
	self.v10 = DoubleVar()
	self.v11 = DoubleVar()
	self.v12 = DoubleVar()
	self.v13 = DoubleVar()
	self.v14 = StringVar()
	self.v9.set('Normal')
	
	Label(self.wprior, text='Prior parameters').grid(row=8, column = 1)
	Label(self.wprior, text='Value').grid(row=7, column = 2)
	Label(self.wprior, text='jump ?').grid(row=7, column = 3)
	Label(self.wprior, text='prior').grid(row=7, column = 4)
	Label(self.wprior, text='Value1').grid(row=7, column = 5)
	Label(self.wprior, text='Value2').grid(row=7, column = 6)
	Label(self.wprior, text='Value3').grid(row=7, column = 7)
	Label(self.wprior, text='Value4').grid(row=7, column = 8)
	Label(self.wprior, text='File').grid(row=7, column = 9)
	
	Entry(self.wprior, textvariable=self.v7, width=10, bg='white', state='disabled').grid(row=8, column=2)
	Checkbutton(self.wprior, variable=self.v8, state='disabled').grid(row=8, column=3)
	OptionMenu(self.wprior, self.v9, *self.list_prior_type).grid(row=8, column=4)
	Entry(self.wprior, textvariable=self.v10, width=10, bg='white').grid(row=8, column=5)
	Entry(self.wprior, textvariable=self.v11, width=10, bg='white').grid(row=8, column=6)
	Entry(self.wprior, textvariable=self.v12, width=10, bg='white').grid(row=8, column=7)
	Entry(self.wprior, textvariable=self.v13, width=10, bg='white').grid(row=8, column=8)
	Entry(self.wprior, textvariable=self.v14, width=10, bg='white').grid(row=8, column=9)

	Button(self.wprior, text='SAVE', state = 'active', command=self.save_prior).grid(row=9, column = 7)
	Button(self.wprior, text='QUIT', command=self.wprior.destroy).grid(row=9, column = 8)

def edit_prior(self):
	self.wprior = Toplevel()
	Label(self.wprior, text='Edit custom prior '+str(self.select),font=('Times', 16, 'bold')).grid(row=0, columnspan = 10)
	self.list_variables = []
	for o in self.list_objects.keys() :
		for p in self.list_objects[o].keys() : 
			self.list_variables.append(o+'_'+p)
	Label(self.wprior, text='x1').grid(row=1, column = 1)
	Label(self.wprior, text='x2').grid(row=2, column = 1)
	Label(self.wprior, text='x3').grid(row=3, column = 1)
	Label(self.wprior, text='x4').grid(row=4, column = 1)
	Label(self.wprior, text='x5').grid(row=5, column = 1)
	self.v1 = StringVar()
	self.v2 = StringVar()
	self.v3 = StringVar()
	self.v4 = StringVar()
	self.v5 = StringVar()

	self.v1.set(self.list_prior[self.select]['variables'][0])
	self.v2.set(self.list_prior[self.select]['variables'][1])
	self.v3.set(self.list_prior[self.select]['variables'][2])
	self.v4.set(self.list_prior[self.select]['variables'][3])
	self.v5.set(self.list_prior[self.select]['variables'][4])
	
	OptionMenu(self.wprior, self.v1, *self.list_variables).grid(row=1, column=2, columnspan=8, sticky=W+E)
	OptionMenu(self.wprior, self.v2, *self.list_variables).grid(row=2, column=2, columnspan=8, sticky=W+E)
	OptionMenu(self.wprior, self.v3, *self.list_variables).grid(row=3, column=2, columnspan=8, sticky=W+E)
	OptionMenu(self.wprior, self.v4, *self.list_variables).grid(row=4, column=2, columnspan=8, sticky=W+E)
	OptionMenu(self.wprior, self.v5, *self.list_variables).grid(row=5, column=2, columnspan=8, sticky=W+E)

	Label(self.wprior, text='Function\n(x1, x2, x3, x4, x5)').grid(row=6, column = 1)
	self.v6 = StringVar()
	for k in self.list_prior[self.select].keys() :
		if k <> 'variables' : 
			self.v6.set(k)
			label = k
	Entry(self.wprior, textvariable = self.v6, width=40, bg='white').grid(row=6, column=2, columnspan=8, sticky=W+E)
	
	self.v7 = DoubleVar()
	self.v8 = IntVar()
	self.v8.set(0)
	self.v9 = StringVar()
	self.v10 = DoubleVar()
	self.v11 = DoubleVar()
	self.v12 = DoubleVar()
	self.v13 = DoubleVar()
	self.v14 = StringVar()
	self.v9 .set(self.list_prior[self.select][label][2])
	self.v10.set(self.list_prior[self.select][label][3])
	self.v11.set(self.list_prior[self.select][label][4])
	self.v12.set(self.list_prior[self.select][label][5])
	self.v13.set(self.list_prior[self.select][label][6])
	self.v14.set(self.list_prior[self.select][label][7])
	
	Label(self.wprior, text='Prior parameters').grid(row=8, column = 1)
	Label(self.wprior, text='Value').grid(row=7, column = 2)
	Label(self.wprior, text='jump ?').grid(row=7, column = 3)
	Label(self.wprior, text='prior').grid(row=7, column = 4)
	Label(self.wprior, text='Value1').grid(row=7, column = 5)
	Label(self.wprior, text='Value2').grid(row=7, column = 6)
	Label(self.wprior, text='Value3').grid(row=7, column = 7)
	Label(self.wprior, text='Value4').grid(row=7, column = 8)
	Label(self.wprior, text='File').grid(row=7, column = 9)
	
	Entry(self.wprior, textvariable=self.v7, width=10, bg='white', state='disabled').grid(row=8, column=2)
	Checkbutton(self.wprior, variable=self.v8, state='disabled').grid(row=8, column=3)
	OptionMenu(self.wprior, self.v9, *self.list_prior_type).grid(row=8, column=4)
	Entry(self.wprior, textvariable=self.v10, width=10, bg='white').grid(row=8, column=5)
	Entry(self.wprior, textvariable=self.v11, width=10, bg='white').grid(row=8, column=6)
	Entry(self.wprior, textvariable=self.v12, width=10, bg='white').grid(row=8, column=7)
	Entry(self.wprior, textvariable=self.v13, width=10, bg='white').grid(row=8, column=8)
	Entry(self.wprior, textvariable=self.v14, width=10, bg='white').grid(row=8, column=9)

	Button(self.wprior, text='SAVE', state = 'active', command=self.save_prior2).grid(row=9, column = 7)
	Button(self.wprior, text='QUIT', command=self.wprior.destroy).grid(row=9, column = 8)



################################################################################################################
################################################################################################################
################################################################################################################

def save_info(self):
	self.infodict['name'] = self.v1.get()
	self.vname.set(self.infodict['name'])
	self.infodict['comment'] = (self.v1b.get()).replace (" ", "_")
	self.infodict['alpha'] = self.v2.get()
	if self.infodict['alpha'] < 0. or self.infodict['alpha'] > 360. :
		tkMessageBox.showerror("Error", "Right Ascension must be between 0 and 360 degrees !\nReturn to default value : 0")
		self.infodict['alpha'] = 0.
	self.infodict['delta'] = self.v3.get()
	if self.infodict['delta'] < -90. or self.infodict['delta'] > 90. :
		tkMessageBox.showerror("Error", "Declination must be between -90 and +90 degrees !\nReturn to default value : 0")
		self.infodict['delta'] = 0.
	self.infodict['MaxDist'] = self.v4.get()
	self.infodict['EXTstep'] = self.v5.get()
	if self.infodict['EXTstep'] > self.infodict['MaxDist'] : 
		tkMessageBox.showerror("Error", "Extinction step must be smaller than the maximum distance !\nReturn to default value : 100")
		self.infodict['EXTstep'] = 100.
	self.infodict['Nmax'] = self.v6.get()
	if self.infodict['Nmax'] == 0 :
		tkMessageBox.showerror("Error", "Number of iterations must be greater than zero !\nReturn to minimum value : 1")
		self.infodict['Nmax'] = 1
	self.infodict['Nchain'] = self.v7.get()
	if self.infodict['Nchain'] == 0 :
		tkMessageBox.showerror("Error", "Number of chains must be greater than zero !\nReturn to minimum value : 1")
		self.infodict['Nchain'] = 1
	self.infodict['beta'] = self.v8.get()
	if self.infodict['beta'] == 0 :
		tkMessageBox.showerror("Error", "Beta parameter must be greater than zero !\nReturn to minimum value : 0.01")
		self.infodict['beta'] = 0.01
	elif self.infodict['beta'] > 1 :
		tkMessageBox.showerror("Error", "Beta parameter must be smaller than one !\nReturn to maximum value : 1")
		self.infodict['beta'] = 1
	self.infodict['Nbeta'] = self.v9.get()
	if self.infodict['Nbeta'] == 0 :
		tkMessageBox.showerror("Error", "Number of different Beta must be greater than zero !\nReturn to minimum value : 1")
		self.infodict['Nbeta'] = 1
	self.infodict['PCA'] = self.v10.get()
	self.infodict['N_update_PCA'] = round(self.v11.get())
	if self.infodict['PCA'] and self.infodict['N_update_PCA'] > self.infodict['Nmax'] :
		tkMessageBox.showerror("Error", "PCA should be updated more than once !\nReturn to default value : 5%")
		self.infodict['N_update_PCA'] = 0.05*self.infodict['Nmax']
	if self.infodict['PCA'] and self.infodict['Min_PCA'] < 1 :
		tkMessageBox.showerror("Error", "PCA can not be updated after zero iteration !\nPCA disabled")
		self.infodict['PCA'] = 0
	self.infodict['Min_PCA'] = round(self.v12.get())
	if self.infodict['PCA'] and self.infodict['Min_PCA'] > self.infodict['Nmax'] :
		tkMessageBox.showerror("Error", "PCA should be started before the end of MCMC !\nReturn to default value : 3%")
		self.infodict['Min_PCA'] = 0.03*self.infodict['Nmax']
	if self.infodict['PCA'] and self.infodict['Min_PCA'] < 1 :
		tkMessageBox.showerror("Error", "PCA can not be started after zero iteration !\nPCA disabled")
		self.infodict['PCA'] = 0
	self.infodict['Max_PCA'] = round(self.v13.get())
	if self.infodict['PCA'] and self.infodict['Min_PCA'] > self.infodict['Max_PCA'] :
		tkMessageBox.showerror("Error", "PCA should be started before stopped !\nReturn to default value : infinity")
		self.infodict['Max_PCA'] = n.inf
	if self.infodict['PCA'] and self.infodict['Max_PCA'] < 1 :
		tkMessageBox.showerror("Error", "PCA can not be stopped after zero iteration !\nPCA disabled")
		self.infodict['PCA'] = 0
	self.infodict['BI_PCA'] = round(self.v14.get())
	self.infodict['randomstart'] = self.v15.get()
	self.infodict['email'] = self.v16.get()
	if len(self.infodict['email']) > 0 and self.infodict['email'].find('@') == -1 :
		tkMessageBox.showerror("Error", "Please fill a valid email address !")
		self.infodict['email'] = ''
	self.infodict['LDC'] = self.v17.get()
	self.infodict['EvolModel'] = self.v18.get()
	self.infodict['SAM'] = self.v19.get()
	self.infodict['walltime'] = self.v20.get()
	self.infodict['save_chain'] = self.v21.get()
	self.infodict['qsub'] = self.v22.get()
	for subdir in [configpath, datapath, resultpath, runpath] :
		if os.path.isdir(os.path.join(subdir, self.infodict['name'])) == 0 :
			os.mkdir(os.path.join(subdir, self.infodict['name']))
		else : continue
	self.winfo.destroy()

def info_box(self):
	"""
	name str
	alpha [deg]
	delta [deg]
	LDCFile str
	ISOPath str
	ZEROMAGFile str
	maxdist [pc]
	extinction step [pc]
	N_max (MCMC)
	Beta parameter (MCMC)
	"""
	self.winfo=Toplevel()
	Label(self.winfo, text='General configuration informations',font=('Times', 16, 'bold')).grid(row=0, columnspan = 6)
	Label(self.winfo, text='Object name').grid(row=1, column = 1)
	Label(self.winfo, text='User comments').grid(row=2, column = 1)
	Label(self.winfo, text='Right Ascension [deg]').grid(row=3, column = 1)
	Label(self.winfo, text='Declination [deg]').grid(row=3, column = 3)
	Label(self.winfo, text='Maximun distance [pc]').grid(row=4, column = 1)
	Label(self.winfo, text='Extinction step [pc]').grid(row=4, column = 3)
	Label(self.winfo, text='# MCMC iterations').grid(row=5, column = 1)
	Label(self.winfo, text='Maximum number of chains').grid(row=5, column = 3)
	Label(self.winfo, text='Beta parameter').grid(row=6, column = 1)
	Label(self.winfo, text='Number of different Beta').grid(row=6, column = 3)
	Label(self.winfo, text='Use Principal Component Analysis ?').grid(row=7, column = 1)
	Label(self.winfo, text='Update PCA after # iterations ?').grid(row=7, column = 3)
	Label(self.winfo, text='Start PCA after # iterations').grid(row=8, column = 1)
	Label(self.winfo, text='End PCA after # iterations').grid(row=8, column = 3)
	Label(self.winfo, text='PCA Burning iterations').grid(row=9, column = 1)
	Label(self.winfo, text='Random starting points in priors ?').grid(row=9, column = 3)
	Label(self.winfo, text='Save .chain file ?').grid(row=10, column = 1)
	Label(self.winfo, text='Use QSUB ?').grid(row=10, column = 3)
	Label(self.winfo, text='MCMC Wall time [h]').grid(row=11, column = 1)
	Label(self.winfo, text='User email address').grid(row=12, column = 1)
	Label(self.winfo, text='Limb Darkening Coefficients').grid(row=13, column = 1)
	Label(self.winfo, text='Evolution Models').grid(row=14, column = 1)
	Label(self.winfo, text='Stellar Atmospheric Models').grid(row=15, column = 1)

	self.v1 = StringVar() #Object name
	self.v1b = StringVar() #Comments
	self.v2 = DoubleVar() #RA
	self.v3 = DoubleVar() #DEC
	self.v4 = DoubleVar() #maxt dist
	self.v5 = DoubleVar() #dist step
	self.v6 = DoubleVar()    #MCMC iterations
	self.v7 = DoubleVar()    #chain number
	self.v8 = DoubleVar() #Beta
	self.v9 = DoubleVar()    #number of different Beta
	self.v10 = IntVar()    #Use PCA ?
	self.v11 = DoubleVar() #PCA update
	self.v12 = DoubleVar() #PCA start
	self.v13 = DoubleVar() #PCA stop
	self.v14 = DoubleVar() #PCA BI
	self.v15 = IntVar()    #random starting point ?
	self.v20 = IntVar()    #Wall time
	self.v21 = IntVar()    #save chain file ?
	self.v22 = IntVar()    #use qsub ?
	self.v16 = StringVar() #email
	self.v17 = StringVar() #LDC
	self.v18 = StringVar() #EvolModel
	self.v19 = StringVar() #SAM

	self.v1.set(self.infodict['name'])
	self.v1b.set(self.infodict['comment'])
	self.v2.set(self.infodict['alpha'])
	self.v3.set(self.infodict['delta'])
	self.v4.set(self.infodict['MaxDist'])
	self.v5.set(self.infodict['EXTstep'])
	self.v6.set(self.infodict['Nmax'])
	self.v7.set(self.infodict['Nchain'])
	self.v8.set(self.infodict['beta'])
	self.v9.set(self.infodict['Nbeta'])
	self.v10.set(self.infodict['PCA'])
	self.v11.set(self.infodict['N_update_PCA'])
	self.v12.set(self.infodict['Min_PCA'])
	self.v13.set(self.infodict['Max_PCA'])
	self.v14.set(self.infodict['BI_PCA'])
	self.v15.set(self.infodict['randomstart'])
	self.v16.set(self.infodict['email'])
	self.v17.set(self.infodict['LDC'])
	self.v18.set(self.infodict['EvolModel'])
	self.v19.set(self.infodict['SAM'])
	try : self.v20.set(self.infodict['walltime'])
	except KeyError : self.v20.set(720.)
	try : self.v21.set(self.infodict['save_chain'])
	except KeyError : self.v21.set(0)
	try : self.v22.set(self.infodict['qsub'])
	except KeyError : self.v22.set(1)

	Entry(self.winfo, textvariable = self.v1, width = 40).grid(row=1, column=2, columnspan = 5, sticky=W+E)
	Entry(self.winfo, textvariable = self.v1b, width = 40).grid(row=2, column=2, columnspan = 5, sticky=W+E)
	Entry(self.winfo, textvariable = self.v2, width = 10).grid(row=3, column=2)
	Entry(self.winfo, textvariable = self.v3, width = 10).grid(row=3, column=4)
	Entry(self.winfo, textvariable = self.v4, width = 10).grid(row=4, column=2)
	Entry(self.winfo, textvariable = self.v5, width = 10).grid(row=4, column=4)
	Entry(self.winfo, textvariable = self.v6, width = 10).grid(row=5, column=2)
	Entry(self.winfo, textvariable = self.v7, width = 10).grid(row=5, column=4)
	Entry(self.winfo, textvariable = self.v8, width = 10).grid(row=6, column=2)
	Entry(self.winfo, textvariable = self.v9, width = 10).grid(row=6, column=4)
	Checkbutton(self.winfo, variable=self.v10).grid(row=7, column=2)
	Entry(self.winfo, textvariable = self.v11, width = 10).grid(row=7, column=4)
	Entry(self.winfo, textvariable = self.v12, width = 10).grid(row=8, column=2)
	Entry(self.winfo, textvariable = self.v13, width = 10).grid(row=8, column=4)
	Entry(self.winfo, textvariable = self.v14, width = 10).grid(row=9, column=2)
	Checkbutton(self.winfo, variable=self.v15).grid(row=9, column=4)
	Checkbutton(self.winfo, variable=self.v21).grid(row=10, column=2)
	Checkbutton(self.winfo, variable=self.v22).grid(row=10, column=4)
	Scale(self.winfo, variable=self.v20, from_=0, to=720, length=300, tickinterval=72, showvalue='yes', orient='horizontal').grid(row=11, column=2, columnspan = 5, sticky=W+E)
	Entry(self.winfo, textvariable = self.v16, width = 10).grid(row=12, column=2, columnspan = 5, sticky=W+E)
	OptionMenu(self.winfo, self.v17, *self.list_LDC).grid(row=13, column=2, columnspan = 5, sticky=W+E)
	OptionMenu(self.winfo, self.v18, *self.list_EvolModel).grid(row=14, column=2, columnspan = 5, sticky=W+E)
	OptionMenu(self.winfo, self.v19, *self.list_SAM).grid(row=15, column=2, columnspan = 5, sticky=W+E)

	Button(self.winfo, text='Save & Quit', command=self.save_info).grid(row=16, column = 3)


################################################################################################################
################################################################################################################
################################################################################################################

def refresh_list(self):
	self.listb.delete(0, END)
	for obj in self.list_objects :
		self.listb.insert(END,obj)
	for obj in self.list_prior :
		self.listb.insert(END,obj)
	self.merge_object_list()

def delete_list(self):
	#self.select=str(self.listb.get(self.listb.curselection()[0]))
	try : del self.list_prior[self.select]
	except : del self.list_objects[self.select]
	try : 
		del self.datafile_list[self.select]
		del self.list_objects_data[self.select]
		del self.list_objects[self.select]
	except KeyError:1
	self.refresh_list()
	self.check_del.destroy()

def clear_list(self):
	self.listb.delete(0, END)
	self.list_objects.clear()
	self.raz_nb()
	self.check_clr.destroy()
	self.refresh_state()

def check_delete(self):
	self.check_del = Toplevel()
	self.select=str(self.listb.get(self.listb.curselection()[0]))
	Label(self.check_del, text='You are going to delete object: '+self.select+'\nAre you sure ?', font=('Times', 14)).pack()
	self.f_del = Frame(self.check_del)
	Button(self.f_del, text='Yes, I am', command=self.delete_list).pack(side='left')
	Button(self.f_del, text='No, I am not', command=self.check_del.destroy).pack()
	self.f_del.pack()
	
def check_clear(self):
	self.check_clr = Toplevel()
	Label(self.check_clr, text='You are going to delete all objects.\nAre you sure ?', font=('Times', 14)).pack()
	self.f_clr = Frame(self.check_clr)
	Button(self.f_clr, text='Yes, I am', command=self.clear_list).pack(side='left')
	Button(self.f_clr, text='No, I am not', command=self.check_clr.destroy).pack()
	self.f_clr.pack()
	

def raz_nb(self):
	self.nb_target=1
	self.nb_star=1
	self.nb_planet=1
	self.nb_binary=1
	self.nb_triple=1
	self.nb_plansys=1
	self.nb_fitobs=1

def merge_object_list(self):
	self.list_objects.update(self.list_objects_data)
	

def backward_compatibility(self):
	if not self.infodict.has_key('walltime') : self.infodict['walltime'] = 720.
	if not self.infodict.has_key('save_chain') : self.infodict['save_chain'] = 0
	if not self.infodict.has_key('qsub') : self.infodict['qsub'] = 1
	for key in self.list_objects_data.keys():
		if key == 'SED' : continue
		if not self.list_objects_data[key].has_key('jitter') : self.list_objects_data[key]['jitter'] = [0., 'Jeffreys', 0.,0.,0.,0.,0.,'']



################################################################################################################
################################################################################################################
################################################################################################################

def refresh_obj_lists(self) :
	for o in self.list_objects.keys() :
		if o[:-1] in ['target', 'star', 'binary', 'plansys'] : self.object_list.append(o)
		if o[:-1] in ['binary', 'plansys'] : self.complex_object_list.append(o)
		if o[:-1] in ['target', 'star']  : self.star_list.append(o)
		if o[:-1] == 'planet' :  self.planet_list.append(o)


def save_list_window(self):
	import os, pickle
	import tkFileDialog

	self.list_filename = tkFileDialog.asksaveasfilename(defaultextension = ".pastis.obj", initialdir = os.path.join(configpath, self.infodict['name']), initialfile = self.infodict['name']+'_'+self.infodict['comment']+'.pastis.obj', filetypes=[('PASTIS Objects', '*.pastis.obj')], title = "Save a PASTIS Objects list file")
	if self.list_filename <> '' :
		ask = True
		#if os.path.isfile(self.list_filename) : ask = tkMessageBox.askokcancel("Overwrite existing file ?", "You are going to overwrite file '"+str(self.list_filename).split('/')[-1]+"'.\n Are you sure ??")
		if os.path.isfile(self.list_filename) == 0 or ask :
			self.fileopen = open(self.list_filename, 'w')
			pickle.dump(self.list_objects, self.fileopen)
			self.fileopen.close()
			tkMessageBox.showinfo("File saved", "File '"+str(self.list_filename).split('/')[-1]+"' saved successfully.\n Great job !")
		else : tkMessageBox.showwarning("File not saved", "File '"+str(self.list_filename).split('/')[-1]+"' NOT saved.")


def load_list_window(self):
	import os, pickle
	import tkFileDialog

	self.list_filename = tkFileDialog.askopenfilename(defaultextension = ".pastis.obj", initialdir = os.path.join(configpath, self.infodict['name']), filetypes=[('PASTIS Objects', '*.pastis.obj')], title = "Load a PASTIS Objects list file")
	if self.pastisfile <> '' :
		self.fileopen = open(self.list_filename, 'r')
		tmp_list = pickle.load(self.fileopen)
		self.list_objects.update(tmp_list)
		self.fileopen.close()
		self.refresh_list()
		self.refresh_obj_lists()
		self.refresh_nb_objects()
		self.refresh_state()
		tkMessageBox.showinfo("File saved", "File '"+str(self.pastis_filename).split('/')[-1]+"' loaded successfully.\n Great job !")



def refresh_state(self) :
	if len(self.star_list) > 1 : 
		self.but4.config(state='active') #active binary button
	if len(self.star_list) > 0 : 
		self.but5c.config(state='disable') #disable fitobs button
	if self.nb_planet > 1 : 
		self.but5c.config(state='disable') #disable fitobs button
	if len(self.list_objects) >= 1 : 
		self.but5b.config(state='active') #active custom prior button
	if len(self.star_list) >= 1 and len(self.planet_list) > 1 : 
		self.but3.config(state='active') #active planetary system button
	if len(self.complex_object_list) >= 1 and len(self.object_list) >= 1 : 
		self.but5.config(state='active') #active multiple system button
	if self.nb_target == 1 and self.nb_star == 1 and self.nb_planet == 1 :
		self.but5c.config(state='active') #active fitobs button
	if self.nb_target == 1 and self.nb_star == 1 and self.nb_planet == 1 and self.nb_fitobs == 1 :
		self.but5b.config(state='disable') #active fitobs button
	self.close_state()
	
def close_state(self):
	if self.nb_fitobs > 1 :
		self.but0.config(state='disable')
		self.but1.config(state='disable')
		self.but2.config(state='disable')
		self.islight = 1
	elif self.nb_fitobs == 1 :
		self.but0.config(state='active')
		self.but1.config(state='active')
		self.but2.config(state='active')
		self.islight = 0

def refresh_nb_objects(self) : 
	self.nb_target=1
	self.nb_star=1
	self.nb_planet=1
	self.nb_binary=1
	self.nb_triple=1
	self.nb_plansys=1
	self.nb_prior=1
	self.nb_fitobs=1
	for o in self.list_objects.keys() : 
		if o[:-1] == 'target' : self.nb_target+=1
		elif o[:-1] == 'star' : self.nb_star+=1   
		elif o[:-1] == 'planet' : self.nb_planet+=1 
		elif o[:-1] == 'binary' : self.nb_binary+=1 
		elif o[:-1] == 'triple' : self.nb_triple+=1 
		elif o[:-1] == 'plansys' : self.nb_plansys+=1
		elif o[:-1] == 'prior' : self.nb_prior+=1  
		elif o[:-1] == 'fitobs' : self.nb_fitobs+=1  


################################################################################################################
################################################################################################################
################################################################################################################


def save_PASTIS_window(self):
	import os, pickle
	import tkFileDialog

	self.pastis_filename = tkFileDialog.asksaveasfilename(defaultextension = ".pastis", initialdir = os.path.join(configpath, self.infodict['name']), initialfile = self.infodict['name']+'_'+self.infodict['comment']+'.pastis', filetypes=[('PASTIS', '*.pastis')], title = "Save a PASTIS configuration file")
	if self.pastis_filename <> '' :
		ask = True
		#if os.path.isfile(self.pastis_filename) : ask = tkMessageBox.askokcancel("Overwrite existing file ?", "You are going to overwrite file '"+str(self.pastis_filename).split('/')[-1]+"'.\n Are you sure ??")
		if os.path.isfile(self.pastis_filename) == 0 or ask :
			self.fileopenp = open(self.pastis_filename, 'w')
			pickle.dump([self.infodict, self.list_objects, self.datafile_list, self.list_prior], self.fileopenp)
			self.fileopenp.close()
			tkMessageBox.showinfo("File saved", "File '"+str(self.pastis_filename).split('/')[-1]+"' saved successfully.\n Great job !")
		else : tkMessageBox.showwarning("File not saved", "File '"+str(self.pastis_filename).split('/')[-1]+"' NOT saved.")


def load_PASTIS_window(self):
	import pickle
	import tkFileDialog
	
	self.pastisfile = tkFileDialog.askopenfilename(defaultextension = ".pastis", initialdir = configpath, filetypes=[('PASTIS', '*.pastis')], title = "Load a PASTIS configuration file")
	if self.pastisfile <> '' :
		self.fileopenpastis = open(self.pastisfile, 'r')
		self.infodict, self.list_objects, self.datafile_list, self.list_prior = pickle.load(self.fileopenpastis)
		self.fileopenpastis.close()
		self.refresh_obj_lists()
		self.refresh_nb_objects()
		self.refresh_state()
		self.refresh_list()
		self.refresh_data_list()
		self.backward_compatibility()
		tkMessageBox.showinfo("File loaded", "File '"+str(self.pastisfile).split('/')[-1]+"' loaded successfully.\n Great job !")
		self.vname.set(self.infodict['name'])
		self.pastis_filename = self.pastisfile


def run_PASTIS_check(self) :
	import tkFileDialog
	import tkMessageBox

	check = tkMessageBox.askyesno(title="save file ?", message='Do you want to save the config file before starting PASTIS ?')
	if check :
		self.save_PASTIS_window()
		self.run_PASTIS()
	else : self.run_PASTIS()

def run_PASTIS(self) :

	from .. import version
	from .. import hostname
	PASTIS = __import__('PASTIS_'+version)
	import time
	import datetime
	import pickle
	import numpy as n
	import os
	
	# Read dictionaries
	#objectdict = self.list_objects
	
	#custompriordict = self.list_prior
	#dd = self.infodict

	# Initialize PASTIS
	try : 
		toto = len(self.pastis_filename)
	except AttributeError :
		try : 
			self.pastis_filename = self.pastisfile
		except AttributeError : 
			self.load_PASTIS_window()
			self.pastis_filename = self.pastisfile
	
	if self.infodict['Nchain'] == 1 and self.infodict['Nbeta'] == 1 and hostname <> 'cluster' and self.infodict['qsub'] == 0 :
		datadict, spectrographs, lightcurves = PASTIS.DataReader(self.datafile_list)
		#if not self.islight : 
		PASTIS.initialize(self.infodict, datadict, self.list_objects)
		#else :
		#	PASTIS.inputdicts=[self.infodict, self.list_objects, self.datafile_list]
		self.w1.destroy()

		from ..MCMC import PASTIS_MCMC

		sufix=datetime.datetime.isoformat(datetime.datetime.now())

		C = PASTIS_MCMC.mcmc(self.list_objects, datadict, self.list_prior, self.infodict['Nmax'],
			     Npca = self.infodict['Min_PCA'], NupdatePCA = self.infodict['N_update_PCA'], BIpca = self.infodict['BI_PCA'],
			     beta = self.infodict['beta'], randomstart = self.infodict['randomstart'], Nlastupdate = self.infodict['Max_PCA'], usePCA = self.infodict['PCA'])
		vd = C.get_value_dict()
		vd['logL'] = C.get_logL()
		vd['posterior'] = C.get_posterior()
		
		sufix=datetime.datetime.isoformat(datetime.datetime.now())
		filename=self.infodict['name']+'_'+(self.infodict['comment']).replace (" ", "_")+'_Beta%.6f_rep'%self.infodict['beta']+str(0).zfill(4)+'_queued'+sufix
		# Define files to save the chain
		chainfile = os.path.join(resultpath, self.infodict['name'], filename+'_started'+sufix+'.chain')
		mcmcfile = chainfile.replace('.chain', '.mcmc')

		# Save value dict to file
		p = open(mcmcfile, 'w')
		pickle.dump(vd, p)
		p.close()

		# Pickle chain to file
		if self.infodict['save_chain'] : C.savetofile(chainfile)

	else : 
		for obj in self.list_objects.keys():
			if self.list_objects[obj].has_key('ebmv') and self.list_objects[obj]['ebmv'][8] == 1:
				from .. import extinction
				ra = self.infodict['alpha']
				dec = self.infodict['delta']
				extinction.initialize_extinction(ra, dec, self.infodict['MaxDist'],
									self.infodict['EXTstep'], Rv = 3.1)
				break


		sufix=datetime.datetime.isoformat(datetime.datetime.now())
		#self.w1.destroy()
		
		betas = n.zeros(self.infodict['Nbeta'], float)
		chains = n.zeros(self.infodict['Nbeta'], float)
		if self.infodict['Nbeta'] == 1 :
			betas += self.infodict['beta']
			chains += self.infodict['Nchain']
		else : 
			betas = n.logspace(-2., 0., self.infodict['Nbeta'])
			chains = n.linspace(1., self.infodict['Nchain'], self.infodict['Nbeta'])
		for beta, nc in zip(betas, chains):
			for c in xrange(round(nc)):

				time.sleep(0.2)
				filename=self.infodict['name']+'_'+(self.infodict['comment']).replace (" ", "_")+'_Beta%.6f_rep'%beta+str(c).zfill(4)+'_queued'+sufix
				f = open(os.path.join(runpath, self.infodict['name'], filename+'.py'), 'w')
				f.write('import pickle\n')
				f.write('import PASTIS_%s as PASTIS\n'%version)
				f.write('import time\n')
				f.write('import os\n')
				f.write('import datetime\n')
				f.write('from numpy import *\n')
				f.write('import smtplib\n\n')
				f.write('time.sleep(4*'+str(c)+')\n') # to avoid file access problems in the cluster 
				f.write("sufix=datetime.datetime.isoformat(datetime.datetime.now())\n\n")
				f.write("f = open('"+self.pastis_filename+"', 'r')\n")
				f.write("dd = pickle.load(f)\n")
				f.write("f.close()\n\n")
				f.write("datadict, sp, lc = PASTIS.DataReader(dd[2])\n\n")
				f.write("input_dict = dd[1].copy()\n\n")
				f.write("PASTIS.initialize(dd[0], datadict, input_dict)\n\n")
				f.write("from PASTIS_%s.MCMC import PASTIS_MCMC\n"%version)
				f.write("C = PASTIS_MCMC.MCMC(input_dict, datadict, dd[3], "+str(self.infodict['Nmax'])+", beta = "+str(beta)+", Npca = "+str(self.infodict['Min_PCA'])+", NupdatePCA = "+str(self.infodict['N_update_PCA'])+", Nlastupdate = "+str(self.infodict['Max_PCA'])+", BIpca = "+str(self.infodict['BI_PCA'])+", randomstart = "+str(self.infodict['randomstart'])+", usePCA = "+str(self.infodict['PCA'])+")\n\n" )
				f.write("vd = C.get_value_dict()\n")
				f.write("vd['logL'] = C.get_logL()\n")
				f.write("vd['posterior'] = C.get_posterior()\n")
				
				f.write("# Save value dict to file\n")
				f.write("p = open(os.path.join(%s, %s, %s+sufix+'.mcmc'), 'w')\n"%( "'"+resultpath+"'", "'"+self.infodict['name']+"'", "'"+filename+"_started'"))
				f.write("pickle.dump(vd, p)\n")
				f.write("p.close()\n")
				f.write("# Pickle chain to file\n")
				if self.infodict['save_chain'] : f.write("C.savetofile(os.path.join(%s, %s, %s+sufix+'.chain'))\n\n"%( "'"+resultpath+"'", "'"+self.infodict['name']+"'", "'"+filename+"_started'"))
				
				if beta == 1. and c == round(nc)-1 and self.infodict['email'] <> '' :
					f.write("sender = 'no-reply@oamp.fr'\n")
					f.write("receivers = ['"+self.infodict['email']+"']\n")
					f.write('message = """From: PASTIS\n')
					f.write('To: '+self.infodict['email']+'\n')
					f.write('Subject: PASTIS simulations status\n')
					f.write('The PASTIS run is successfully terminated.\n')
					f.write('\n')
					f.write('Object name: '+self.infodict['name']+'\n')
					f.write('Comments: '+self.infodict['comment']+'\n')
					f.write('Configuration file: '+self.pastis_filename+'\n')
					f.write('\n')
					f.write('Thank you for using PASTIS - use with(out) moderation\n')
					f.write('\n\n')
					f.write('This email is generated automatically by PASTIS, please do not reply."""\n')
					f.write('try :\n')
					f.write("   smtpObj = smtplib.SMTP('smtps.oamp.fr')\n")
					f.write("   smtpObj.sendmail(sender, receivers, message)\n")
					f.write("except SMTPException:\n")
					f.write('   print "Error: unable to send email"\n')
					f.close()
				else : 
					f.close()
				
				qsubname=(self.infodict['name']+' '+self.infodict['comment']).replace (" ", "_")
				qsubsubmit='qsub '+os.path.join(runpath,'run_PASTIS_qsub.sh')+' -v file="'+self.infodict['name']+"/"+filename+'.py" -N '+qsubname+' -o '+os.path.join(runpath,self.infodict['name'],filename+'.output')+' -e '+os.path.join(runpath,self.infodict['name'],filename+'.error'+' -l walltime='+str(int(self.infodict['walltime']))+':00:00')
				#print qsubsubmit 
				os.system(qsubsubmit)
	return 



def compute_prior(self, priortype, params, x):
	'''
	prior = priors.Prior(priortype, *params)
	return prior.prob(x)
	'''
	return 1

"""
	from scipy import stats
	if priortype == 'File':
		# Read file containing distribution
		x, p = n.loadtxt(params[-1], usecols = (0, 1), unpack = True)
		
		# Normalise prior
		p = p/sum(p*(x[1] - x[0]))

		# Interpolate values
		prob = interpolate.interp1d(x, p, 'linear', bounds_error = False, fill_value = 0.0)
		
		return self.prob(x)
	
	if priortype == 'Uniform': return stats.uniform.pdf(x, params[0], params[1]-params[0])
		
	if priortype == 'Normal': return stats.norm.pdf(x, params[0], params[1])
		
	if priortype == 'LogNormal': return stats.norm.pdf(n.log10(x), params[0], params[1])
		
	if priortype == 'Jeffreys': return 1.0/(x*log(params[1]/params[0]))
		
	if priortype == 'Binormal': return 0.5*(stats.norm.pdf(x, params[0], params[1]) + stats.norm.pdf(x, params[2], params[3]))
"""
