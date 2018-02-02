from Tkinter import *
import tkMessageBox
import matplotlib
import numpy as n
matplotlib.use('TkAgg')
import os, os.path

from .. import widgetpath, libpath, configpath, datapath, resultpath, runpath, PriorTypeList



#PASTIS = __import__('PASTIS_'+version)
from ..MCMC import priors as priors
from .. import photometry
#from .. import plots

imgpath = os.path.join(libpath, 'Images')

#from object_tools import define_target, define_star, define_planet

class Application() :
    execfile(os.path.join(widgetpath, 'object_tools.py'))
    execfile(os.path.join(widgetpath, 'config_tools.py'))
    execfile(os.path.join(widgetpath, 'data_tools.py'))
    execfile(os.path.join(widgetpath, 'result_tools.py'))

    def __init__(self) :
        # Init objects counter
	self.nb_target = 1
	self.nb_star = 1
	self.nb_planet = 1
	self.nb_binary = 1
	self.nb_triple = 1
	self.nb_plansys = 1
	self.nb_fitobs = 1
	self.nb_prior = 1
	self.nb_plot = 0
	self.nb_spectro = 1
	self.nb_ccf = 1
		
	# Init objects lists
	self.list_objects = {}
	self.list_objects_data = {}
	self.object_list = []
	self.complex_object_list = []
	self.star_list = []
	self.planet_list = []
	self.fitobs_list = []
	self.planet_list.append('None')
	
	# Define photometric instruments
	self.list_phot = photometry.Allpbands 
	self.list_phot.append('CoRoT-R')
	self.list_phot.append('CoRoT-G')
	self.list_phot.append('CoRoT-B')
	#['CoRoT-W', 'CoRoT-R', 'CoRoT-G', 'CoRoT-B', 'Kepler', 'IRAC-I1', 'IRAC-I2', 'IRAC-I3', 'IRAC-I4', 'SDSS-U', 'SDSS-G', 'SDSS-R', 'SDSS-I', 'SDSS-Z', 'Johnson-U', 'Johnson-B', 'Johnson-V', 'Johnson-R', 'Johnson-I', '2MASS-J', '2MASS-H', '2MASS-Ks', 'STROMGREN-u', 'STROMGREN-v', 'STROMGREN-b', 'STROMGREN-y', 'Bessell-U', 'Bessell-B', 'Bessell-V', 'Bessell-R', 'Bessell-I', 'CAMELOT-U', 'CAMELOT-B', 'CAMELOT-V', 'CAMELOT-R', 'CAMELOT-I']
	self.nb_photband = {}
	for i in self.list_phot : 
		self.nb_photband[i] = 1
	self.list_prior_type = PriorTypeList

	# Init datafile list
	self.datafile_list = {}

	# Init datafile list
	self.list_prior = {}
		
	# Init info dictionnary
	self.infodict = {}
	self.infodict['name'] = 'CoRoT-42'
	self.infodict['comment'] = ''
	self.infodict['alpha'] = 0.0000
	self.infodict['delta'] = 0.0000
	self.infodict['MaxDist'] = 20000.
	self.infodict['EXTstep'] = 100.
	self.infodict['Nmax'] = 1e6
	self.infodict['Nchain'] = 1 #Integer
	self.infodict['beta'] = 1. # >0 <=1
	self.infodict['Nbeta'] = 1 #Integer
	self.infodict['PCA'] = 1 #boolean
	self.infodict['N_update_PCA'] = 50000
	self.infodict['Min_PCA'] = 30000
	self.infodict['Max_PCA'] = n.inf
	self.infodict['BI_PCA'] = 10000
	self.infodict['randomstart'] = 0 #boolean
	self.infodict['email'] = ''
	self.infodict['walltime'] = 720
	self.infodict['save_chain'] = 0
	self.infodict['qsub'] = 1
	self.infodict['LDC'] = 'Claret2011'
	self.infodict['EvolModel'] = 'Dartmouth'
	self.infodict['SAM'] = 'BT-settl'
	
	#Initiate info lists
	self.list_LDC = ['Claret2011']
	self.list_EvolModel = ["Padova", "Dartmouth", "Geneva", "StarEvol"]
	self.list_SAM = ['CastelliKurucz', 'BT-settl']
	
	
	#Check for PASTISlight
	self.islight = False
	
	self.w1 = Tk()  # Start main window
	self.w1.title('PASTIS user tool')  #set main window title
	self.image_file = PhotoImage(file=imgpath+"/toolbarOpenFile.gif")
	self.photo = PhotoImage(file=imgpath+"/PASTIS_logo.gif")  # define the logo file
	self.labl = Label(self.w1, image = self.photo) # display the logo
	self.labl.pack()
		
	self.lab = Label(self.w1, text="PASTIS configuration file tool", font=('Times', 14, 'bold')) # set a subtitle
	self.lab.pack()
		
	Label(self.w1, text="General configuration information:", font=('Times', 12, 'bold')).pack() # set a subtitle
	Label(self.w1, text="Object of Interest:", font=('Times', 12)).pack() # set a subtitle
	self.vname = StringVar()
	Label(self.w1, textvariable=self.vname, font=('Times', 12)).pack() # set the target name
	self.vname.set(self.infodict['name'])
		
	self.f0 = Frame(self.w1)  # define general information frame
	Button(self.f0, text='PASTIS Configuration', command=self.info_box).pack()  # add a button
	self.f0.pack()
		
################################################################################################################
	self.pane = PanedWindow(orient=HORIZONTAL)  # separate the main window in two sub-windows
	self.pane.pack(fill=BOTH, expand=1)
	self.left_frame = Frame(self.w1)  # define the left frame
	self.right_frame = Frame(self.w1)  # define the right frame
	
	Label(self.left_frame, text="Define scenario", font=('Times', 12, 'bold')).pack()
	self.f1 = Frame(self.left_frame)  # define the stars frame
	self.f2 = Frame(self.left_frame)  # define the planet frame
	self.f3 = Frame(self.left_frame)  # define the system frame
	self.f3b = Frame(self.left_frame)  # define the advanced prior frame
	
	self.but0 = Button(self.f1,text = "Add a target star", state = 'active',command = self.define_target)  # add the target button
	self.but1 = Button(self.f1,text="Add a star", state = 'active', command=self.define_star)  # add the star button
	self.but2 = Button(self.f2,text="Add a planet", state = 'active', command=self.define_planet)  # add the planet button
	self.but3 = Button(self.f2,text="Make a planetary system", state = 'disabled', command=self.define_plansys)  # add plansys button
	self.but4 = Button(self.f3,text="Make a binary", state = 'disable', command=self.define_binary)  # add the binary button
	self.but5 = Button(self.f3,text="Make a multiple system", state = 'disabled', command=self.define_triple)  # add the triple button
	self.but5b = Button(self.f3b,text="Add a custom prior", state = 'disabled', command=self.define_prior)  # add the prior button
	self.but5c = Button(self.f3b,text="Fit obs module", state = 'active', command=self.define_fitobs)  # add the prior button
		
	# put all previous Buttons in respective frames
	
	self.but0.pack(side="left")
	self.but1.pack()
	self.but2.pack(side="left")
	self.but3.pack()
	self.but4.pack(side="left")
	self.but5.pack()
	self.but5b.pack(side='left')
	self.but5c.pack()

	# display frames
		
	self.f1.pack()
	self.f2.pack()
	self.f3.pack()
	self.f3b.pack()
	
	# Add the Listbox containing all objects
		
	Label(self.left_frame,text="List of objects:\n(double-click to edit objects)").pack()
	self.listb=Listbox(self.left_frame, bg='white')  # define the listbox
	self.listb.pack()
		
	# Add listbox options buttons
	
	self.f4 = Frame(self.left_frame)
	Button(self.f4, text='Refresh', command=self.refresh_list).pack(side='left')
	Button(self.f4, text='Delete', command=self.check_delete).pack(side='left')
	Button(self.f4, text='Clear', command=self.check_clear).pack()
	self.f4.pack()
		
	# Add save/load buttons
		
	self.f5 = Frame(self.left_frame)
	Button(self.f5, text='Save object list', command=self.save_list_window).pack(side='left')
	Button(self.f5, text='Load object list', command=self.load_list_window).pack(side='left')
	self.f5.pack()
		
	# define double-click command
		
	self.listb.bind('<Double-1>', self.edit_list)
		
	# display objects sub-window
		
	self.pane.add(self.left_frame)
	self.pane.paneconfigure(self.left_frame, minsize=300)

################################################################################################################
		
	# Add Data file button
		
	Label(self.right_frame, text="Add data files", font=('Times', 12, 'bold')).pack()
	self.brv = Button(self.right_frame, text='RV Data File', command=self.rv_data_box)
	self.bphot = Button(self.right_frame, text='Phot Data File', command=self.phot_data_box)
	self.bsed = Button(self.right_frame, text='SED Data File', command=self.sed_data_box)
	self.bccf = Button(self.right_frame, text='CCF Data File', state = 'active', command=self.ccf_data_box)
		
	self.brv.pack()
	self.bphot.pack()
	self.bsed.pack()
	self.bccf.pack()
	
	# Add listbox containing data files
		
	Label(self.right_frame,text="List of data files:\n(double-click to edit data files)").pack()
	self.listbdata=Listbox(self.right_frame, bg='white')
	self.listbdata.pack()
		
	# Add listbox options buttons
		
	self.f7 = Frame(self.right_frame)
	Button(self.f7, text='Refresh', command=self.refresh_data_list).pack(side='left')
	Button(self.f7, text='Delete', command=self.check_data_delete).pack(side='left')
	Button(self.f7, text='Clear', command=self.check_data_clear).pack()
	self.f7.pack()
		
	# Add save/load buttons
		
	self.f8 = Frame(self.right_frame)
	Button(self.f8, text='Save data list', command=self.save_data_list_window).pack(side='left')
	Button(self.f8, text='Load data list', command=self.load_data_list_window).pack()
	self.f8.pack()
		
	# define double-click command
		
	self.listbdata.bind('<Double-1>', self.edit_data_list)
		
	# Display datafiles sub-windows
		
	self.pane.add(self.right_frame)
	self.pane.paneconfigure(self.right_frame, minsize=300)

################################################################################################################
		
	# Add PASTIS save/load buttons
		
	Label(self.w1, text="PASTIS tools", font=('Times', 12, 'bold')).pack()
	self.fpastis = Frame(self.w1)
	self.but6 = Button(self.fpastis, text = 'Save Config File', command = self.save_PASTIS_window)
	self.but6.pack(side='left')

	self.but7 = Button(self.fpastis, text = 'Load Config File', state = 'active', command = self.load_PASTIS_window)
	self.but7.pack(side='left')

	self.fpastis.pack()
		
	# Add PASTIS run/exit buttons
		
	self.fexit = Frame(self.w1)
	self.r=Button(self.fexit,text="Run PASTIS", font=('Times', 12, 'bold'), state = 'active', width = 20, command=self.run_PASTIS_check)
	self.r.pack(side='left')
	self.result = Button(self.fexit, text='Plot results', font=('Times', 12, 'bold'), width = 20, state = 'active', command = self.load_result_window)
	self.result.pack(side='left')
	self.q=Button(self.fexit,text="Exit PASTIS", font=('Times', 12, 'bold'), width = 20, command=self.w1.destroy)
	self.q.pack()
	self.fexit.pack()
		
		
	# star main window main loop
		
	self.w1.mainloop()

