'''
2012-04-25 Code debugging
'''

from Tkinter import *
import tkMessageBox
import tkFileDialog


def save_rv_data(self, quit = None):
	"""
	Save RV data into datafile dictionnary
	"""
	import tkMessageBox

	if self.spectrograph_flag.get() == 0 : self.spectrograph = 'HARPS'
	elif self.spectrograph_flag.get() == 1 : self.spectrograph = 'SOPHIE HE'
	elif self.spectrograph_flag.get() == 2 : self.spectrograph = 'SOPHIE HR'
	elif self.spectrograph_flag.get() == 3 : self.spectrograph = 'RotProfile'
	elif self.spectrograph_flag.get() == 4 : 
		self.spectrograph = 'Other'+str(self.nb_spectro)
		self.nb_spectro+=1
	
	if self.mask_flag.get() == 0 : self.mask='G2'
	elif self.mask_flag.get() == 1 : self.mask='K5'

	if self.v11.get() and self.compute_prior(self.v12.get(), [self.v13.get(), self.v14.get(), self.v15.get(), self.v16.get(), self.v17.get()], self.v10.get()) == 0 :
		tkMessageBox.showerror("Error", "First instrument offset value out of prior.\n Please put a correct value within the prior range")
	else : 
		self.datafile_list[self.spectrograph+'-'+self.mask] = {'datafile' : self.rvdatafile.get(), 'instrument' : self.spectrograph, 'is_RV' : self.rv.get(), 'is_BIS' : self.bis.get(), 'is_FWHM' : self.fwhm.get(), 'is_CONT' : self.cont.get(), 'type' : 'RV'}
		self.list_objects_data[self.spectrograph+'-'+self.mask] = {'offset' : [self.v10.get(), self.v11.get(), self.v12.get(), self.v13.get(), self.v14.get(), self.v15.get(), self.v16.get(), self.v17.get()], 'jitter' : [self.v20.get(), self.v21.get(), self.v22.get(), self.v23.get(), self.v24.get(), self.v25.get(), self.v26.get(), self.v27.get()]}
		self.refresh_data_list()
		self.refresh_list()
		self.merge_object_list()
		
	if quit: self.wrv.destroy()

def save_phot_data(self, quit = None, save_edit = None):
	"""
	Save Phot data into datafile dictionnary
	"""
	import tkMessageBox

	if self.v4.get()%1 <> 0 : tkMessageBox.showwarning("Warning", "Oversampling factor must be an integer.\n New value : "+str(int(self.v4.get())))
	if self.v11.get() and self.compute_prior(self.v12.get(), [self.v13.get(), self.v14.get(), self.v15.get(), self.v16.get(), self.v17.get()], self.v10.get()) == 0 :
		tkMessageBox.showerror("Error", "First contamination value out of prior.\n Please put a correct value within the prior range")
	elif self.v21.get() and self.compute_prior(self.v22.get(), [self.v23.get(), self.v24.get(), self.v25.get(), self.v26.get(), self.v27.get()], self.v20.get()) == 0 :
		tkMessageBox.showerror("Error", "First Out-of-Transit flux value out of prior.\n Please put a correct value within the prior range")
	else : 
		if save_edit:
			self.datafile_list[self.select] = {'datafile' : self.photdatafile.get(), 'filter' : self.instrument.get(), 'is_phase' : self.time.get(), 'MeanFlux' : self.v3.get(), 'type' : 'PHOT', 'sampling' : int(self.v4.get())}
			self.list_objects_data[self.select] = {'contamination' : [self.v10.get(), self.v11.get(), self.v12.get(), self.v13.get(), self.v14.get(), self.v15.get(), self.v16.get(), self.v17.get()], 'foot' : [self.v20.get(), self.v21.get(), self.v22.get(), self.v23.get(), self.v24.get(), self.v25.get(), self.v26.get(), self.v27.get()], 'jitter' : [self.v30.get(), self.v31.get(), self.v32.get(), self.v33.get(), self.v34.get(), self.v35.get(), self.v36.get(), self.v37.get()]}
					
			self.refresh_data_list()
			self.refresh_list()
			self.merge_object_list()
		else:
			self.datafile_list[self.instrument.get()+str(self.nb_photband[self.instrument.get()])] = {'datafile' : self.photdatafile.get(), 'filter' : self.instrument.get(), 'is_phase' : self.time.get(), 'MeanFlux' : self.v3.get(), 'type' : 'PHOT', 'sampling' : int(self.v4.get())}
			self.list_objects_data[self.instrument.get()+str(self.nb_photband[self.instrument.get()])] = {'contamination' : [self.v10.get(), self.v11.get(), self.v12.get(), self.v13.get(), self.v14.get(), self.v15.get(), self.v16.get(), self.v17.get()], 'foot' : [self.v20.get(), self.v21.get(), self.v22.get(), self.v23.get(), self.v24.get(), self.v25.get(), self.v26.get(), self.v27.get()], 'jitter' : [self.v30.get(), self.v31.get(), self.v32.get(), self.v33.get(), self.v34.get(), self.v35.get(), self.v36.get(), self.v37.get()]}
			
			self.nb_photband[self.instrument.get()]+=1
			self.refresh_data_list()
			self.refresh_list()
			self.merge_object_list()
			
		if quit: self.wphot.destroy()
	
def save_sed_data(self):
	"""
	Save SED data into datafile dictionnary
	"""
	self.datafile_list['SED'] = {'datafile' : self.seddatafile.get(), 'type' : 'SED'}
	self.refresh_data_list()
	self.wsed.destroy()

def save_ccf_data(self, save_edit = None):
	"""
	Save SED data into datafile dictionnary
	"""
	if self.spectrograph_flag.get() == 0 : self.spectrograph = 'HARPS'
	elif self.spectrograph_flag.get() == 1 : self.spectrograph = 'SOPHIE HE'
	elif self.spectrograph_flag.get() == 2 : self.spectrograph = 'SOPHIE HR'
	elif self.spectrograph_flag.get() == 3 : self.spectrograph = 'RotProfile'
	elif self.spectrograph_flag.get() == 4 : 
		self.spectrograph = 'Other'+str(self.nb_spectro)
		self.nb_spectro+=1
	
	if self.mask_flag.get() == 0 : self.mask='G2'
	elif self.mask_flag.get() == 1 : self.mask='K5'

	self.datafile_list['CCF'+self.nb_ccf] = {'datafile' : self.ccfdatafile.get(), 'instrument' : self.spectrograph, 'mask' : self.mask, 'type' : 'CCF'}
	if save_edit: 
		self.datafile_list['CCF'+self.nb_ccf] = {'datafile' : self.ccfdatafile.get(), 'instrument' : self.spectrograph, 'mask' : self.mask, 'type' : 'CCF'}
	else:
		self.nb_ccf+=1
	
	self.refresh_data_list()
	self.wccf.destroy()

def open_rv_file(self):
	import tkFileDialog
	self.rvdatafile.set(tkFileDialog.askopenfilename(initialdir = os.path.join(datapath,self.infodict['name']), title = "Load a RV data file"))

def open_phot_file(self):
	import tkFileDialog
	self.photdatafile.set(tkFileDialog.askopenfilename(initialdir = os.path.join(datapath,self.infodict['name']), title = "Load a photometric data file"))

def open_sed_file(self):
	import tkFileDialog
	self.seddatafile.set(tkFileDialog.askopenfilename(initialdir = os.path.join(datapath,self.infodict['name']), title = "Load a SED data file"))

def open_ccf_file(self):
	import tkFileDialog
	self.ccfdatafile.set(tkFileDialog.askopenfilename(initialdir = os.path.join(datapath,self.infodict['name']), title = "Load a CCF data file"))

################################################################################################################
################################################################################################################
################################################################################################################


def rv_data_box(self, edit = None):
	"""
	Display RV data file widget
	"""
	import tkFileDialog, tkMessageBox
	self.wrv=Toplevel()
	Label(self.wrv, text='Add a Radial Velocity data file',font=('Times', 16, 'bold')).grid(row=0, columnspan = 9)
	Label(self.wrv, text='Filename').grid(row=1, column=1)
	self.rvdatafile = StringVar()
	if edit: 
		self.rvdatafile.set(self.datafile_list[self.select]['datafile'])
	Entry(self.wrv, textvariable=self.rvdatafile, width=40, bg='white').grid(row=1, column=2, columnspan=6, sticky=W+E)
	Button(self.wrv, image=self.image_file, command=self.open_rv_file).grid(row=1, column=8)
	Label(self.wrv, text='Spectrograph').grid(row=2, column=1)
	self.spectrograph_flag = IntVar()
	if edit:
		if self.datafile_list[self.select]['instrument'] == 'HARPS' : self.spectrograph_flag.set(0)
		elif self.datafile_list[self.select]['instrument'] == 'SOPHIE HE' : self.spectrograph_flag.set(1)
		elif self.datafile_list[self.select]['instrument'] == 'SOPHIE HR' : self.spectrograph_flag.set(2)
		elif self.datafile_list[self.select]['instrument'] == 'RotProfile' : self.spectrograph_flag.set(3)
		elif self.datafile_list[self.select]['instrument'] == 'Other' : self.spectrograph_flag.set(4)
	Radiobutton(self.wrv, text = 'HARPS', variable=self.spectrograph_flag, value=0).grid(row=2, column=2)
	Radiobutton(self.wrv, text = 'SOPHIE HE', variable=self.spectrograph_flag, value=1).grid(row=2, column=3)
	Radiobutton(self.wrv, text = 'SOPHIE HR', variable=self.spectrograph_flag, value=2).grid(row=2, column=4)
	Radiobutton(self.wrv, text = 'RotProfile', variable=self.spectrograph_flag, state = 'disabled', value=3).grid(row=2, column=5)
	Radiobutton(self.wrv, text = 'Other', variable=self.spectrograph_flag, state = 'disabled', value=4).grid(row=2, column=6)
	Label(self.wrv, text='Correlation mask').grid(row=3, column=1)
	self.mask_flag = IntVar()
	if edit:
		if self.select.split('-')[1] == 'G2' : self.mask_flag.set(0)
		elif self.select.split('-')[1] == 'K5' : self.mask_flag.set(1)
	Radiobutton(self.wrv, text = 'G2', variable=self.mask_flag, value=0).grid(row=3, column=2)
	Radiobutton(self.wrv, text = 'K5', variable=self.mask_flag, value=1).grid(row=3, column=3)
	Label(self.wrv, text='Data type').grid(row=4, column=1)
	self.rv = IntVar()
	self.bis = IntVar()
	self.fwhm = IntVar()
	self.cont = IntVar()
	if edit:
		self.rv.set(self.datafile_list[self.select]['is_RV'])
		self.bis.set(self.datafile_list[self.select]['is_BIS'])
		self.fwhm.set(self.datafile_list[self.select]['is_FWHM'])
		self.cont.set(self.datafile_list[self.select]['is_CONT'])
	Checkbutton(self.wrv, text='RV', variable=self.rv).grid(row=4, column=2)
	Checkbutton(self.wrv, text='BIS', variable=self.bis).grid(row=4, column=3)
	Checkbutton(self.wrv, text='FWHM', variable=self.fwhm).grid(row=4, column=4)
	Checkbutton(self.wrv, text='Cont', variable=self.cont).grid(row=4, column=5)

	self.v10 = DoubleVar()
	self.v11 = IntVar()
	self.v12 = StringVar()
	self.v13 = DoubleVar()
	self.v14 = DoubleVar()
	self.v15 = DoubleVar()
	self.v16 = DoubleVar()
	self.v17 = StringVar()

	self.v20 = DoubleVar()
	self.v21 = IntVar()
	self.v22 = StringVar()
	self.v23 = DoubleVar()
	self.v24 = DoubleVar()
	self.v25 = DoubleVar()
	self.v26 = DoubleVar()
	self.v27 = StringVar()

	Label(self.wrv, text='Parameter').grid(row=5, column=1)
	Label(self.wrv, text='Value').grid(row=5, column=2)
	Label(self.wrv, text='jump ?').grid(row=5, column=3)
	Label(self.wrv, text='prior').grid(row=5, column=4)
	Label(self.wrv, text='Value1').grid(row=5, column=5)
	Label(self.wrv, text='Value2').grid(row=5, column=6)
	Label(self.wrv, text='Value3').grid(row=5, column=7)
	Label(self.wrv, text='Value4').grid(row=5, column=8)
	Label(self.wrv, text='File').grid(row=5, column=9)

	if edit:
		self.v10.set(self.list_objects[self.select]['offset'][0])
		self.v11.set(self.list_objects[self.select]['offset'][1])
		self.v12.set(self.list_objects[self.select]['offset'][2])
		self.v13.set(self.list_objects[self.select]['offset'][3])
		self.v14.set(self.list_objects[self.select]['offset'][4])
		self.v15.set(self.list_objects[self.select]['offset'][5])
		self.v16.set(self.list_objects[self.select]['offset'][6])
		self.v17.set(self.list_objects[self.select]['offset'][7])
		try : 
			self.v20.set(self.list_objects[self.select]['jitter'][0])
			self.v21.set(self.list_objects[self.select]['jitter'][1])
			self.v22.set(self.list_objects[self.select]['jitter'][2])
			self.v23.set(self.list_objects[self.select]['jitter'][3])
			self.v24.set(self.list_objects[self.select]['jitter'][4])
			self.v25.set(self.list_objects[self.select]['jitter'][5])
			self.v26.set(self.list_objects[self.select]['jitter'][6])
			self.v27.set(self.list_objects[self.select]['jitter'][7])
		except KeyError :
			tkMessageBox.showwarning(title = 'Problem with jitter Value', message = 'An error occurs with the jitter parameter, please check its value and prior')
			self.v22.set('Jeffreys')

	else : 
		self.v12.set('Normal')
		self.v22.set('Jeffreys')

	Label(self.wrv, text='Instr. offset [km/s]').grid(row=6, column=1)
	Entry(self.wrv, textvariable=self.v10, width=10, bg='white').grid(row=6, column=2)
	Checkbutton(self.wrv, variable=self.v11).grid(row=6, column=3)
	OptionMenu(self.wrv, self.v12, *self.list_prior_type).grid(row=6, column=4)
	Entry(self.wrv, textvariable=self.v13, width=10, bg='white').grid(row=6, column=5)
	Entry(self.wrv, textvariable=self.v14, width=10, bg='white').grid(row=6, column=6)
	Entry(self.wrv, textvariable=self.v15, width=10, bg='white').grid(row=6, column=7)
	Entry(self.wrv, textvariable=self.v16, width=10, bg='white').grid(row=6, column=8)
	Entry(self.wrv, textvariable=self.v17, width=10, bg='white').grid(row=6, column=9)

	Label(self.wrv, text='Jitter [km/s]').grid(row=7, column=1)
	Entry(self.wrv, textvariable=self.v20, width=10, bg='white').grid(row=7, column=2)
	Checkbutton(self.wrv, variable=self.v21).grid(row=7, column=3)
	OptionMenu(self.wrv, self.v22, *self.list_prior_type).grid(row=7, column=4)
	Entry(self.wrv, textvariable=self.v23, width=10, bg='white').grid(row=7, column=5)
	Entry(self.wrv, textvariable=self.v24, width=10, bg='white').grid(row=7, column=6)
	Entry(self.wrv, textvariable=self.v25, width=10, bg='white').grid(row=7, column=7)
	Entry(self.wrv, textvariable=self.v26, width=10, bg='white').grid(row=7, column=8)
	Entry(self.wrv, textvariable=self.v27, width=10, bg='white').grid(row=7, column=9)

	Button(self.wrv, text='SAVE & QUIT', command=lambda: self.save_rv_data(quit=1)).grid(row=8, column = 6)
	Button(self.wrv, text='SAVE', command=self.save_rv_data).grid(row=8, column = 7)
	Button(self.wrv, text='QUIT', command=self.wrv.destroy).grid(row=8, column = 8)


def phot_data_box(self, edit = None):
	"""
	Display Phot data file widget
	"""
	import tkFileDialog
	self.wphot = Toplevel()
	Label(self.wphot, text='Add a Light Curve data file',font=('Times', 16, 'bold')).grid(row=0, columnspan = 10)
	Label(self.wphot, text='Filename').grid(row=1, column=1)
	self.photdatafile = StringVar()
	if edit: 
		self.photdatafile.set(self.datafile_list[self.select]['datafile'])
	Entry(self.wphot, textvariable=self.photdatafile, width=40, bg='white').grid(row=1, column=2, columnspan=7, sticky=W+E)
	Button(self.wphot, image=self.image_file, command=self.open_phot_file).grid(row=1, column=9)
	Label(self.wphot, text='Instrument').grid(row=2, column=1)
	self.instrument = StringVar()
	if edit: 
		self.instrument.set(self.datafile_list[self.select]['filter'])
	else: 
		self.instrument.set('CoRoT-W')

	OptionMenu(self.wphot, self.instrument, *self.list_phot).grid(row=2, column=2, columnspan = 9, sticky=W+E)
	
	self.time = IntVar()
	if edit: 
		self.time.set(self.datafile_list[self.select]['is_phase'])
	Label(self.wphot, text='Data type').grid(row=3, column=1)
	Radiobutton(self.wphot, text = 'phase', variable=self.time, value=1).grid(row=3, column=2, columnspan=2)
	Radiobutton(self.wphot, text = 'time', variable=self.time, value=0).grid(row=3, column=4, columnspan=2)

	Label(self.wphot, text='Parameter').grid(row=5, column=1)
	Label(self.wphot, text='Value').grid(row=5, column=2)
	Label(self.wphot, text='jump ?').grid(row=5, column=3)
	Label(self.wphot, text='prior').grid(row=5, column=4)
	Label(self.wphot, text='Value1').grid(row=5, column=5)
	Label(self.wphot, text='Value2').grid(row=5, column=6)
	Label(self.wphot, text='Value3').grid(row=5, column=7)
	Label(self.wphot, text='Value4').grid(row=5, column=8)
	Label(self.wphot, text='File').grid(row=5, column=9)
	Label(self.wphot, text='Contamination').grid(row=6, column=1)
	Label(self.wphot, text='OutOfTransit flux').grid(row=7, column=1)
	Label(self.wphot, text='Jitter').grid(row=8, column=1)
	Label(self.wphot, text='Mean CoRoT-RGB flux').grid(row=9, column=1)
	Label(self.wphot, text='Oversampling factor').grid(row=9, column=4)
	
	self.v10 = DoubleVar()
	self.v11 = IntVar()
	self.v12 = StringVar()
	self.v13 = DoubleVar()
	self.v14 = DoubleVar()
	self.v15 = DoubleVar()
	self.v16 = DoubleVar()
	self.v17 = StringVar()
	self.v20 = DoubleVar()
	self.v21 = IntVar()
	self.v22 = StringVar()
	self.v23 = DoubleVar()
	self.v24 = DoubleVar()
	self.v25 = DoubleVar()
	self.v26 = DoubleVar()
	self.v27 = StringVar()
	self.v30 = DoubleVar()
	self.v31 = IntVar()
	self.v32 = StringVar()
	self.v33 = DoubleVar()
	self.v34 = DoubleVar()
	self.v35 = DoubleVar()
	self.v36 = DoubleVar()
	self.v37 = StringVar()
	self.v3 = DoubleVar()
	self.v4 = IntVar()

	if edit:
		self.v10.set(self.list_objects[self.select]['contamination'][0])
		self.v11.set(self.list_objects[self.select]['contamination'][1])
		self.v12.set(self.list_objects[self.select]['contamination'][2])
		self.v13.set(self.list_objects[self.select]['contamination'][3])
		self.v14.set(self.list_objects[self.select]['contamination'][4])
		self.v15.set(self.list_objects[self.select]['contamination'][5])
		self.v16.set(self.list_objects[self.select]['contamination'][6])
		self.v17.set(self.list_objects[self.select]['contamination'][7])
		self.v20.set(self.list_objects[self.select]['foot'][0])
		self.v21.set(self.list_objects[self.select]['foot'][1])
		self.v22.set(self.list_objects[self.select]['foot'][2])
		self.v23.set(self.list_objects[self.select]['foot'][3])
		self.v24.set(self.list_objects[self.select]['foot'][4])
		self.v25.set(self.list_objects[self.select]['foot'][5])
		self.v26.set(self.list_objects[self.select]['foot'][6])
		self.v27.set(self.list_objects[self.select]['foot'][7])
		self.v3.set(self.datafile_list[self.select]['MeanFlux'])
		self.v4.set(self.datafile_list[self.select]['sampling'])
		try : 
			self.v30.set(self.list_objects[self.select]['jitter'][0])
			self.v31.set(self.list_objects[self.select]['jitter'][1])
			self.v32.set(self.list_objects[self.select]['jitter'][2])
			self.v33.set(self.list_objects[self.select]['jitter'][3])
			self.v34.set(self.list_objects[self.select]['jitter'][4])
			self.v35.set(self.list_objects[self.select]['jitter'][5])
			self.v36.set(self.list_objects[self.select]['jitter'][6])
			self.v37.set(self.list_objects[self.select]['jitter'][7])
		except KeyError :
			tkMessageBox.showwarning(title = 'Problem with jitter Value', message = 'An error occurs with the jitter parameter, please check its value and prior')
			self.v32.set('Jeffreys')

	else:
		self.v12.set('Normal')
		self.v22.set('Normal')
		self.v4.set(1)
		self.v32.set('Jeffreys')


	Entry(self.wphot, textvariable=self.v10, width=10, bg='white').grid(row=6, column=2)
	Checkbutton(self.wphot, variable=self.v11).grid(row=6, column=3)
	OptionMenu(self.wphot, self.v12, *self.list_prior_type).grid(row=6, column=4)
	Entry(self.wphot, textvariable=self.v13, width=10, bg='white').grid(row=6, column=5)
	Entry(self.wphot, textvariable=self.v14, width=10, bg='white').grid(row=6, column=6)
	Entry(self.wphot, textvariable=self.v15, width=10, bg='white').grid(row=6, column=7)
	Entry(self.wphot, textvariable=self.v16, width=10, bg='white').grid(row=6, column=8)
	Entry(self.wphot, textvariable=self.v17, width=10, bg='white').grid(row=6, column=9)
	
	Entry(self.wphot, textvariable=self.v20, width=10, bg='white').grid(row=7, column=2)
	Checkbutton(self.wphot, variable=self.v21).grid(row=7, column=3)
	OptionMenu(self.wphot, self.v22, *self.list_prior_type).grid(row=7, column=4)
	Entry(self.wphot, textvariable=self.v23, width=10, bg='white').grid(row=7, column=5)
	Entry(self.wphot, textvariable=self.v24, width=10, bg='white').grid(row=7, column=6)
	Entry(self.wphot, textvariable=self.v25, width=10, bg='white').grid(row=7, column=7)
	Entry(self.wphot, textvariable=self.v26, width=10, bg='white').grid(row=7, column=8)
	Entry(self.wphot, textvariable=self.v27, width=10, bg='white').grid(row=7, column=9)
	
	Entry(self.wphot, textvariable=self.v30, width=10, bg='white').grid(row=8, column=2)
	Checkbutton(self.wphot, variable=self.v31).grid(row=8, column=3)
	OptionMenu(self.wphot, self.v32, *self.list_prior_type).grid(row=8, column=4)
	Entry(self.wphot, textvariable=self.v33, width=10, bg='white').grid(row=8, column=5)
	Entry(self.wphot, textvariable=self.v34, width=10, bg='white').grid(row=8, column=6)
	Entry(self.wphot, textvariable=self.v35, width=10, bg='white').grid(row=8, column=7)
	Entry(self.wphot, textvariable=self.v36, width=10, bg='white').grid(row=8, column=8)
	Entry(self.wphot, textvariable=self.v37, width=10, bg='white').grid(row=8, column=9)
	
	Entry(self.wphot, textvariable=self.v3, width=10, bg='white').grid(row=9, column=2)
	Entry(self.wphot, textvariable=self.v4, width=10, bg='white').grid(row=9, column=5)

	
	Button(self.wphot, text='SAVE & QUIT', command=lambda: self.save_phot_data(quit=1,save_edit=edit)).grid(row=13, column = 6)
	Button(self.wphot, text='SAVE', command=lambda: self.save_phot_data(save_edit=edit)).grid(row=13, column = 7)
	Button(self.wphot, text='QUIT', command=self.wphot.destroy).grid(row=13, column = 8)
	       

def sed_data_box(self, edit = None):
	"""
	Display SED data file widget
	"""
	import tkFileDialog

	self.wsed = Toplevel()
	Label(self.wsed, text='Add a SED data file',font=('Times', 16, 'bold')).grid(row=0, columnspan = 6)
	Label(self.wsed, text='Filename').grid(row=1, column=1)
	self.seddatafile = StringVar()
	if edit: self.seddatafile.set(self.datafile_list['SED']['datafile'])
	Entry(self.wsed, textvariable=self.seddatafile, width=40, bg='white').grid(row=1, column=2, columnspan=3, sticky=W+E)
	Button(self.wsed, image=self.image_file, command=self.open_sed_file).grid(row=1, column=5)
	Button(self.wsed, text='SAVE & QUIT', command=self.save_sed_data).grid(row=2, column = 3, columnspan=2)
	#Button(self.wsed, text='QUIT', command=self.wsed.destroy).grid(row=2, column = 4)


def ccf_data_box(self, edit = None):
	"""
	Display CCF data file widget
	"""
	import tkFileDialog

	self.wccf = Toplevel()
	if edit:
		Label(self.wccf, text='Edit a CCF data file',font=('Times', 16, 'bold')).grid(row=0, columnspan = 6)
	else:
		Label(self.wccf, text='Add a CCF data file',font=('Times', 16, 'bold')).grid(row=0, columnspan = 6)
	Label(self.wccf, text='Filename').grid(row=1, column=1)
	self.ccfdatafile = StringVar()
	if edit: self.ccfdatafile.set(self.datafile_list[self.select]['datafile'])
	Entry(self.wccf, textvariable=self.ccfdatafile, width=40, bg='white').grid(row=1, column=2, columnspan=3, sticky=W+E)
	Button(self.wccf, image=self.image_file, command=self.open_ccf_file).grid(row=1, column=5)
	self.spectrograph_flag = IntVar()
	Radiobutton(self.wccf, text = 'HARPS', variable=self.spectrograph_flag, value=0).grid(row=2, column=2)
	Radiobutton(self.wccf, text = 'SOPHIE HE', variable=self.spectrograph_flag, value=1).grid(row=2, column=3)
	Radiobutton(self.wccf, text = 'SOPHIE HR', variable=self.spectrograph_flag, value=2).grid(row=2, column=4)
	Radiobutton(self.wccf, text = 'RotProfile', variable=self.spectrograph_flag, state = 'disabled', value=3).grid(row=2, column=5)
	Radiobutton(self.wccf, text = 'Other', variable=self.spectrograph_flag, state = 'disabled', value=4).grid(row=2, column=6)
	Label(self.wccf, text='Correlation mask').grid(row=3, column=1)
	self.mask_flag = IntVar()
	Radiobutton(self.wccf, text = 'G2', variable=self.mask_flag, value=0).grid(row=3, column=2)
	Radiobutton(self.wccf, text = 'K5', variable=self.mask_flag, value=1).grid(row=3, column=3)

	Label(self.wccf, text='Parameter').grid(row=4, column=1)
	Label(self.wccf, text='Value').grid(row=4, column=2)
	Label(self.wccf, text='jump ?').grid(row=4, column=3)
	Label(self.wccf, text='prior').grid(row=4, column=4)
	Label(self.wccf, text='Value1').grid(row=4, column=5)
	Label(self.wccf, text='Value2').grid(row=4, column=6)
	Label(self.wccf, text='Value3').grid(row=4, column=7)
	Label(self.wccf, text='Value4').grid(row=4, column=8)
	Label(self.wccf, text='File').grid(row=4, column=9)

	Label(self.wccf, text='Instr. offset [km/s]').grid(row=5, column=1)
	Entry(self.wccf, textvariable=self.v10, width=10, bg='white').grid(row=5, column=2)
	Checkbutton(self.wccf, variable=self.v11).grid(row=5, column=3)
	OptionMenu(self.wccf, self.v12, *self.list_prior_type).grid(row=5, column=4)
	Entry(self.wccf, textvariable=self.v13, width=10, bg='white').grid(row=5, column=5)
	Entry(self.wccf, textvariable=self.v14, width=10, bg='white').grid(row=5, column=6)
	Entry(self.wccf, textvariable=self.v15, width=10, bg='white').grid(row=5, column=7)
	Entry(self.wccf, textvariable=self.v16, width=10, bg='white').grid(row=5, column=8)
	Entry(self.wccf, textvariable=self.v17, width=10, bg='white').grid(row=5, column=9)


	Button(self.wccf, text='SAVE & QUIT', command=lambda: self.save_ccf_data(save_edit=edit)).grid(row=6, column = 3, columnspan=2)


################################################################################################################
################################################################################################################
################################################################################################################

def edit_data_list(self, e):
	"""
	Edition datafile function
	"""
	
	self.select=str(self.listbdata.get(self.listbdata.curselection()[0]))
	if self.select == 'SED' : self.sed_data_box(edit=1)
	elif self.select.split('-')[0] in ['HARPS', 'SOPHIE HE', 'SOPHIE HR', 'RotProfil', 'Other'] : self.rv_data_box(edit=1)
	elif self.select[:-1] in self.list_phot : self.phot_data_box(edit=1,)
	elif self.select[:-1] == 'CCF' : self.ccf_data_box(edit=1)


################################################################################################################
################################################################################################################
################################################################################################################


def refresh_data_list(self):
	"""
	Update the datafile listbox
	"""
	self.listbdata.delete(0, END)
	for obj in self.datafile_list :
		self.listbdata.insert(END,obj)

def delete_data_list(self):
	"""
	Delete a datafile / function
	"""
	del self.datafile_list[self.select]
	try : 
		del self.list_objects[self.select]
		del self.list_objects_data[self.select]
	except KeyError:1
	self.refresh_data_list()
	self.refresh_list()
	self.check_data_del.destroy()


def check_data_delete(self):
	"""
	Delete a datafile / widget
	"""
	self.check_data_del = Toplevel()
	self.select=str(self.listbdata.get(self.listbdata.curselection()[0]))
	Label(self.check_data_del, text='You are going to delete data file: '+self.select+'\nAre you sure ?', font=('Times', 14)).pack()
	self.f_data_del = Frame(self.check_data_del)
	Button(self.f_data_del, text='Yes, I am', command=self.delete_data_list).pack(side='left')
	Button(self.f_data_del, text='No, I am not', command=self.check_data_del.destroy).pack()
	self.f_data_del.pack()
	

def clear_data_list(self):
	"""
	Clear datafiles / function
	"""
	self.listbdata.delete(0, END)
	self.datafile_list.clear()
	self.check_data_clr.destroy()

def check_data_clear(self):
	"""
	Clear datafiles / widget
	"""
	self.check_data_clr = Toplevel()
	Label(self.check_data_clr, text='You are going to delete all data files.\nAre you sure ?', font=('Times', 14)).pack()
	self.f_data_clr = Frame(self.check_data_clr)
	Button(self.f_data_clr, text='Yes, I am', command=self.clear_data_list).pack(side='left')
	Button(self.f_data_clr, text='No, I am not', command=self.check_data_clr.destroy).pack()
	self.f_data_clr.pack()

################################################################################################################

def save_data_list_window(self):
	"""
	Save datafiles / widget
	"""
	import os, pickle
	import tkFileDialog

	self.list_data_filename = tkFileDialog.asksaveasfilename(defaultextension = ".pastis.dat", initialdir = '/data/PASTIS/configfiles/'+self.infodict['name'], initialfile = self.infodict['name']+'_'+self.infodict['comment']+'.pastis.dat', filetypes=[('PASTIS data', '*.pastis.dat')], title = "Save a PASTIS datafile list file")
	if self.list_data_filename <> '' :
		ask = True
		#if os.path.isfile(self.list_data_filename) : ask = tkMessageBox.askokcancel("Overwrite existing file ?", "You are going to overwrite file '"+str(self.list_data_filename).split('/')[-1]+"'.\n Are you sure ??")
		if os.path.isfile(self.list_data_filename) == 0 or ask :
			self.fileopend = open(self.list_data_filename, 'w')
			pickle.dump([self.datafile_list, self.list_objects_data], self.fileopend)
			self.fileopend.close()
			tkMessageBox.showinfo("File saved", "File '"+str(self.list_data_filename).split('/')[-1]+"' saved successfully.\n Great job !")
		else : tkMessageBox.showwarning("File not saved", "File '"+str(self.list_data_filename).split('/')[-1]+"' NOT saved.")


def load_data_list_window(self):
	import os, pickle
	import tkFileDialog

	self.list_data_filename = tkFileDialog.askopenfilename(defaultextension = ".pastis.dat", initialdir = '/data/PASTIS/configfiles/'+self.infodict['name'], filetypes=[('PASTIS data', '*.pastis.dat')], title = "Load a PASTIS datafile list file")
	if self.list_data_filename <> '' :
		self.fileopen_data = open(self.list_data_filename, 'r')
		tmp_list1, tmp_list2 = pickle.load(self.fileopen_data)
		self.datafile_list.update(tmp_list1)
		self.list_objects_data.update(tmp_list2)
		self.fileopen_data.close()
		self.refresh_data_list()
		self.merge_object_list()
		tkMessageBox.showinfo("File saved", "File '"+str(self.list_data_filename).split('/')[-1]+"' loaded successfully.\n Great job !")

