import os
import pickle
from .paths import runpath, resultpath

def run_sim(pastisfile, pastisversion, submit=True, profiler=False,
            fullsave=False):

    # leer pastisfile
    f = open(pastisfile, 'rb')
    dd = pickle.load(f)
    f.close()

    infodict = dd[0]
    list_objects = dd[1]
    datafile_list = dd[2]
    list_prior = dd[3]

    PASTIS = __import__(pastisversion)
    import time
    import datetime
    import numpy as n

    sufix = datetime.datetime.isoformat(datetime.datetime.now())

    betas = n.zeros(infodict['Nbeta'], float)
    chains = n.zeros(infodict['Nbeta'], float)
    if infodict['Nbeta'] == 1:
        betas += infodict['beta']
        chains += infodict['Nchain']
    else:
        # betas[1]=1e-6
        # betas[2]=1e-4
        betas[0] = 1e-4
        betas[1:] = n.logspace(-2., 0., infodict['Nbeta'] - 1.)
        chains = n.linspace(1., infodict['Nchain'], infodict['Nbeta'])
    for beta, nc in zip(betas, chains):
        for c in range(int(round(nc))):

            time.sleep(0.5)
            filename = infodict['name'] + '_' + (infodict['comment']).replace(
                " ", "_")
            f = open(os.path.join(runpath, infodict['name'], filename + '.py'),
                     'w')
            f.write('\'\'\'\n')
            f.write('This file was automatically generated by PASTIS.\n')
            f.write('Do not modify it unless you know what you are doing.\n')
            f.write('\'\'\'\n')
            f.write('import time\n')
            f.write('import os\n')
            f.write('import datetime\n')
            f.write('import pickle\n')
            f.write('import %s as PASTIS\n' % pastisversion)
            f.write('from {}.DataTools import readdata'
                    '\n'.format(pastisversion))
            #f.write('from numpy import *\n')
            f.write('time.sleep(10*' + str(
                c) + ')\n')  # to avoid file access problems in the cluster
            
            f.write('PASTISroot = os.getenv(\'PASTISPATH\')\n')
            f.write(
                "sufix = datetime.datetime.isoformat(datetime.datetime.now())\n\n")
            f.write("outfile = os.path.join(PASTISroot, \'resultfiles\', "
                    "\'{}\',\'{}\'+sufix+'.mcmc')\n\n".format(infodict['name'], 
                                                              filename+"_started"))
            
            f.write("f = open('" + pastisfile + "', 'rb')\n")
            f.write("dd = pickle.load(f)\n")
            f.write("f.close()\n\n")
            f.write("datadict, lc = readdata(dd[2])\n\n")
            f.write("input_dict = dd[1].copy()\n\n")
            f.write("PASTIS.initialize(dd[0], datadict, input_dict)\n\n")
            f.write("from %s.MCMC import PASTIS_MCMC\n" % pastisversion)

            if infodict['Nsave'] is not None:
                # Prepare command with preliminary dumping of results
                f.write(
                    "outfilepart = outfile.replace('.mcmc', '.part.mcmc')\n\n")
                if profiler:
                    f.write("import cProfile\n\n")
                    f.write("outproffile = outfile.replace"
                            "('.mcmc','.prof')\n\n")
                    f.write(
                        "cProfile.runctx('C = PASTIS_MCMC.mcmc(input_dict, datadict, dd[3], " + str(
                            infodict['Nmcmc']) + ", beta=" + str(
                            beta) + ", Npca=" + str(
                            infodict['Min_PCA']) + ", NupdatePCA=" + str(
                            infodict[
                                'N_update_PCA']) + ", Nlastupdate=" + str(
                            infodict['Max_PCA']) + ", BIpca=" + str(
                            infodict['BI_PCA']) + ", randomstart=" + str(
                            infodict['randomstart']) + ", usePCA=" + str(
                            infodict['PCA']) + ", Nsave=" + str(infodict[
                                                                      'Nsave']) + ", outputfile=outfilepart)', globals(), locals(),  filename = outproffile)\n\n")
                else:
                    runparams = {'steps': int(infodict['Nmcmc']),
                                 'beta': beta,
                                 'PCAstart': int(infodict['Min_PCA']),
                                 'PCAupdate': int(infodict['N_update_PCA']),
                                 'PCAstop': int(infodict['Max_PCA']),
                                 'BI_PCA': int(infodict['BI_PCA']),
                                 'randomstart': infodict['randomstart'],
                                 'usePCA': infodict['PCA'],
                                 'Nsave': int(infodict['Nsave'])
                                 }
                    
                    f.write("C = PASTIS_MCMC.mcmc(input_dict, datadict, dd[3],"
                            " {steps:d}, beta={beta}, Npca={PCAstart:d}, "
                            "NupdatePCA={PCAupdate:d}, Nlastupdate={PCAstop:d}"
                            ", BIpca={BI_PCA:d}, randomstart={randomstart}, "
                            "usePCA={usePCA}, Nsave={Nsave:d}, "
                            "outputfile=outfilepart)\n\n".format(**runparams))

            else:
                if profiler:
                    f.write("import cProfile\n\n")
                    f.write(
                        "outproffile = outfile.replace('.mcmc','.prof')\n\n")
                    f.write(
                        "cProfile.runctx('C = PASTIS_MCMC.MCMC(input_dict, datadict, dd[3], " + str(
                            infodict['Nmcmc']) + ", beta = " + str(
                            beta) + ", Npca = " + str(
                            infodict['Min_PCA']) + ", NupdatePCA = " + str(
                            infodict[
                                'N_update_PCA']) + ", Nlastupdate = " + str(
                            infodict['Max_PCA']) + ", BIpca = " + str(
                            infodict['BI_PCA']) + ", randomstart = " + str(
                            infodict['randomstart']) + ", usePCA = " + str(
                            infodict[
                                'PCA']) + ")', globals(), locals(),  filename = outproffile)\n\n")
                else:
                    runparams = {'steps': int(infodict['Nmcmc']),
                                 'beta': beta,
                                 'PCAstart': int(infodict['Min_PCA']),
                                 'PCAupdate': int(infodict['N_update_PCA']),
                                 'PCAstop': int(infodict['Max_PCA']),
                                 'BI_PCA': int(infodict['BI_PCA']),
                                 'randomstart': infodict['randomstart'],
                                 'usePCA': infodict['PCA'],
                                 }
                    f.write("C = PASTIS_MCMC.mcmc(input_dict, datadict, dd[3],"
                            " {steps:d}, beta={beta}, Npca={PCAstart:d}, "
                            "NupdatePCA={PCAupdate:d}, Nlastupdate={PCAstop:d}"
                            ", BIpca={BI_PCA:d}, randomstart={randomstart}, "
                            "usePCA={usePCA}\n\n".format(**runparams))

            f.write("vd = C.get_value_dict()\n")
            f.write("vd['logL'] = C.get_logL()\n")
            f.write("vd['posterior'] = C.get_posterior()\n")

            f.write("# Save value dict to file\n")
            if 'gzip' in infodict and infodict['gzip']:
                f.write("import gzip\n")
                f.write("p = gzip.open(outfile+'.gz', 'w')\n")
            else:
                f.write("p = open(outfile, 'wb')\n")

            f.write("pickle.dump(vd, p)\n")
            f.write("p.close()\n")

            if infodict['Nsave'] != None and infodict['Nmcmc'] >= infodict[
                'Nsave']:
                ## Erase preliminary file
                f.write("os.remove(outfilepart)")

            if fullsave:
                f.write("# Pickle full chain to file\n")
                f.write(
                    "C.savetofile(os.path.join(%s, %s, %s+sufix+'.chain'))\n\n" % (
                    "'" + resultpath + "'", "'" + infodict['name'] + "'",
                    "'" + filename + "_started'"))

            f.close()

            if submit:
                qsubname = (infodict['name'] + ' ' + infodict[
                    'comment']).replace(" ", "_") + '_rep' % beta + str(
                    c).zfill(4)
                # qsubname=(infodict['name']+' '+infodict['comment']).replace (" ", "_")
                # qsubsubmit='qsub -V '+os.path.join(runpath,'run_PASTIS_qsub.sh')+' -v file="'+infodict['name']+"/"+filename+'.py" -N '+qsubname+' -o '+os.path.join(runpath,infodict['name'],filename+'.output')+' -e '+os.path.join(runpath,infodict['name'],filename+'.error'+' -l walltime=%d:00:00'%infodict['walltime'])
                qsubsubmit = 'sbatch -p {queue} -J {jobname} -t {walltime} ' \
                             '-o {pyfile}.out -e {pyfile}.error {cshfile} ' \
                             '{pyfile}.py'.format(
                    cshfile=os.path.join(runpath, 'run_PASTIS_qsub.sh'),
                    pyfile=os.path.join(runpath, infodict['name'], filename),
                    jobname=qsubname,
                    walltime=int(infodict.pop('walltime', 7 * 24) * 60.0),
                    queue=infodict['queue'])

                os.system(qsubsubmit)

    return
