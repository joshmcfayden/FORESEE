import numpy as np
from src.foresee import Foresee, Utility, Model
import os
from matplotlib import pyplot as plt
import random, math

def unweight_events(energies, thetas, weights, number):
    
    #initialize arrays and random number generator
    random.seed()
    unweighted_thetas=[]
    unweighted_energies=[]
    unweighted_weights=[]
    total_weight = np.sum(weights)
    event_weight = total_weight/float(number)
    
    #unweighting
    for irand in range (number):
        stopweight = random.random()*total_weight
        partid, weightsum = 0, 0
        while weightsum < stopweight: 
            weightsum+=weights[partid]
            partid+=1
        unweighted_thetas.append(thetas[partid])
        unweighted_energies.append(energies[partid])
        unweighted_weights.append(event_weight)
        
    return  np.array(unweighted_energies), np.array(unweighted_thetas), np.array(unweighted_weights)


run_plotpions=False
run_setupmodel=True
run_LLPspectra=False
run_rateexample=True
run_setupscans=False
run_plotseparations=False
run_runscans=False
run_plotreach=False

scan_search=None
scan_search=["F2","F-default-sep"]
scan_search=["S1","F-default-sep"]
scan_search=["S2","F-default-sep"]
#scan_search=["S3","F-default-sep"]
scan_search=["F-default-sep-station0"]
#scan_search=["TEST_1m0T","TEST_1m1T","TEST_1p5m0p5T"]
scan_search=["TEST_1m0T","TEST_1p5m0p5T"]

#############
# Initialization
print("INFO   : Initialise FORESEE")
foresee = Foresee()


#############
# Plot pions
if run_plotpions:
    print("INFO   : Plot pions")
    plot=foresee.get_spectrumplot(pid="111", generator="EPOSLHC", energy="14")
    #plot.show()
    plot.savefig("NewConfigs_v3-DarkPhoton-EPOSLHC_Pion-Angle_vs_Momentum.pdf")

    
#############
# Specifying the Model: Dark Photons

if run_setupmodel:    
    print("INFO   : Setting up Dark Photon model")
    energy = "14"
    modelname="DarkPhoton"
    model = Model(modelname)
    
    print("INFO   :   - Adding production modes")
    
    ## pion decay production
    model.add_production_2bodydecay(
        pid0 = "111",
        pid1 = "22",
        br = "2.*0.99 * coupling**2 * pow(1.-pow(mass/self.masses('111'),2),3)",
        generator = "EPOSLHC",
        energy = energy,
        nsample = 10
    )
    
    ## eta decay production
    model.add_production_2bodydecay(
        pid0 = "221",
        pid1 = "22",
        br = "2.*0.39 * coupling**2 * pow(1.-pow(mass/self.masses('221'),2),3)",
        generator = "EPOSLHC",
        energy = energy,
        nsample = 10, 
    )
    
    ## Resonant mixing with SM V bosons production
    model.add_production_mixing(
        pid = "113",
        mixing = "coupling * 0.3/5. * 0.77545**2/abs(mass**2-0.77545**2+0.77545*0.147*1j)",
        generator = "EPOSLHC",
        energy = energy,
    )
    
    ## Bremsstrahlung production
    model.add_production_direct(
        label = "Brem",
        energy = energy,
        condition = "p.pt<1",
        coupling_ref=1,
    )
    
    
    ## Drell-Yan production
    model.add_production_direct(
        label = "DY",
        energy = energy,
        coupling_ref=1,
        massrange=[1.5, 10.]
    )
    
    print("INFO   :   - Setting lifetimes")
    ## Lifetime
    model.set_ctau_1d(
        filename="files/models/"+modelname+"/ctau.txt", 
        coupling_ref=1
    )
    
    print("INFO   :   - Setting branching fractions")
    ## Branching ratio
    model.set_br_1d(
        modes=["e_e", "mu_mu"],
        filenames=["files/models/"+modelname+"/br/e_e.txt","files/models/"+modelname+"/br/mu_mu.txt"]
    )
    
    ## Set model just created
    foresee.set_model(model=model)
        

masses = [ 
    0.01  ,  0.0126,  0.0158,  0.02  ,  0.0251,  0.0316,  0.0398,
    0.0501,  0.0631,  0.0794,  0.1   ,  0.1122,  0.1259,  0.1413,
    0.1585,  0.1778,  0.1995,  0.2239,  0.2512,  0.2818,  0.3162,
    0.3548,  0.3981,  0.4467,  0.5012,  0.5623,  0.6026,  0.631 ,
    0.6457,  0.6607,  0.6761,  0.6918,  0.7079,  0.7244,  0.7413,
    0.7586,  0.7762,  0.7943,  0.8128,  0.8318,  0.8511,  0.871 ,
    0.8913,  0.912 ,  0.9333,  0.955 ,  0.9772,  1.    ,  1.122 ,
    1.2589,  1.4125,  1.5849,  1.7783,  1.9953,  2.2387,  2.5119,
    2.8184,  3.1623,  3.9811,  5.0119,  6.3096,  7.9433, 10.    
]

masses = [ 0.01  ,  0.1,  1.0, 10. ]

    
##########
# Generate LLP Spectra
if run_LLPspectra:    
    print("INFO   : Generating LLP Spectra")
    

    
    print("INFO   :   - Generating for m_A=100 MeV and epsilon=10^-5")
    ## Look at benchmark scenario with m_A=100 MeV and epsilon=10^-5
    plt = foresee.get_llp_spectrum(0.1, coupling=10**(-5), do_plot=True)
    #plt.show()
    plt.savefig("NewConfigs_v3-DarkPhoton-EPOSLHC_LLP_m100_e10m5-Angle_vs_Momentum.pdf")
    
    print("INFO   :   - Generating for all masses")
    for mass in masses:
        foresee.get_llp_spectrum(mass=mass,coupling=1)
    
    
    
#################
# Count Event Rate in Detector
if run_rateexample:    
    print("INFO   : Count Event Rate in Detector")
    
    ## To count the #decays in detector volume need detector geometry
    ## These are FASER2 defaults
    distance, selection, length, luminosity, channels, geo = 480, "np.sqrt(x.x**2 + x.y**2)< 1", 5, 3000, None, None #[(0.5,1.5)]
    foresee.set_detector(distance=distance, selection=selection, length=length, luminosity=luminosity, channels=channels,geo=geo)
    
    
    print("INFO   :   - For dark photon (m_A'=100 MeV) check #events in decay volume")
    ## For one dark photon (m_A'=100 MeV) look at how many particles decay inside the decay volume.
    mass = 0.1
    print(np.logspace(-8,-3,6))
    output = foresee.get_events(mass=0.1, energy=energy, couplings=np.logspace(-8,-3,3),nsample=500 )
    coups, ctaus, nsigs, energies, weights, seps, thetas = output
    print("TEST1",len(energies))
    print("TEST2",len(energies[0]))
    
    weighted_energies, weighted_thetas, weighted_weights = energies[0], thetas[0], weights[0]
    unweighted_energies, unweighted_thetas, unweighted_weights = unweight_events (energies[0], thetas[0], weights[0], int(len(energies[0])/10))
   
    print("TEST3",len(unweighted_energies))
    
    datax=[]
    datay=[]
    dataz=[]
    datae=[]
    dataw=[]
    #print "INFO:   - Unweighted energies":
    #print unweighted_energies
    
    # get events
    f= open("output.hepmc","w")
    f.write("\nHepMC::Version 2.06.09\n")
    f.write("HepMC::IO_GenEvent-START_EVENT_LISTING\n")
    for i,(theta,energy,weight) in enumerate(zip(unweighted_thetas, unweighted_energies, unweighted_weights)):
        
        #Event Info
        f.write("E %s -1 -1.0000000000000000e+00 -1.0000000000000000e+00 -1.0000000000000000e+00 20 0 1 1 2 0 0\n"%(i+1))
        #f.write("N 1 \"0\" \n")
        f.write("U GEV MM\n")
        #f.write("C "+str(weight)+" 0 \n")
        
        #vertex
        posz = random.uniform(0,length)*1000.
        phi  = random.uniform(-math.pi,math.pi)
        posy = theta*480.*1000*np.sin(phi)
        posx = theta*480.*1000*np.cos(phi) 
        post = 3.*10**8 * np.sqrt(posz**2 + posy**2 + posz**2)
        f.write("V -1 0 ")
        f.write(str(posx)+" ")
        f.write(str(posy)+" ")
        f.write(str(posz)+" ")
        f.write(str(post)+" ")
        f.write("0 2 0\n")
        datax.append(posx)
        datay.append(posy)
        dataz.append(posz)
        datae.append(energy)
        dataw.append(weight)
        
        #particles
        particles = foresee.decay_llp(mass=0.1, energy=energy)
        for n,particle in enumerate(particles):  
            if n == 0 :
                f.write("P "+str(n+1)+" 11 ")
            else:
                f.write("P "+str(n+1)+" -11 ")
            f.write(str(particle.px)+" ")
            f.write(str(particle.py)+" ")
            f.write(str(particle.pz)+" ")
            f.write(str(particle.e)+" ")
            f.write("0 1 0 0 0 0\n")
    f.write('HepMC::IO_GenEvent-END_EVENT_LISTING\n')
    f.close()

    
    for coup,ctau,nsig in zip(coups, ctaus, nsigs):
        print ("epsilon =", '{:5.3e}'.format(coup), ": nsignal =", '{:5.3e}'.format(nsig))
    
    print("INFO   :   - Plot energy distribution for different couplings")
    ## Look at energy distribution of the dark photons which decay inside the detector
    fig = plt.figure(figsize=(7,5))
    ax = plt.subplot(1,1,1)
    for coup,en,weight in zip(coups,energies,weights):
        print("en:",len(en),en[0])
        print("coup:",coup)
        print("weight:",len(weight),weight[0])
        
        if sum(weight)<10**-5 : continue
        ax.hist(en, weights=weight, bins=np.logspace(2,4, 20), histtype='step', label=r"$\epsilon=$"+str(coup)) 
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(10**-7,10**5) 
        ax.set_xlabel("E(A') [GeV]") 
        ax.set_ylabel("Number of Events per bin") 
        ax.legend(frameon=False, labelspacing=0)
    plt.tight_layout()
    #plt.show()
    plt.savefig("NewConfigs_v3-DarkPhoton-EPOSLHC_LLP_m100-Energy.pdf")

    #fig = plt.figure(figsize=(7,5))
    #ax = plt.subplot(1,1,1)
    #for coup,sep,weight in zip(coups,seps,weights):
    #    tmp_s=[]
    #    tmp_w=[]
    #    #sep_cut=300.e-6
    #    sep_cut=0 # no cut
    #    for n,s in enumerate(sep):
    #        tmp_s.extend([e for e in s if e>sep_cut])
    #        tmp_w.extend([weight[n] for e in s if e>sep_cut])
    #    #sep=tmp_s
    #    #weight=tmp_w
    #    print(len(tmp_s),len(tmp_w))
    #    print("sep:",len(sep),sep[0])
    #    print("coup:",coup)
    #    print("weight:",len(weight),weight[0])
    #    
    #    if sum(weight)<10**-5 : continue
    #    ax.hist(tmp_s, weights=tmp_w, bins=np.logspace(-5,-1, 40), histtype='step', label=r"$\epsilon=$"+str(coup)) 
    #    ax.set_xscale("log")
    #    ax.set_yscale("linear")
    #    ax.set_ylim(0,7000) 
    #    #ax.set_yscale("log")
    #    #ax.set_ylim(10**-7,10**5) 
    #    ax.set_xlabel("Separation [m]") 
    #    ax.set_ylabel("Number of Events per bin") 
    #    ax.legend(frameon=False, labelspacing=0)
    #plt.tight_layout()
    ##plt.show()
    #plt.savefig("NewConfigs_v3-DarkPhoton-EPOSLHC_LLP_m100-Separation.pdf")





##########
# Parameter Scan
setup_dict={}
if run_setupscans:
    print("INFO   : Setup parameter scan for different detector configs")
    
    ## Get the LLP sensitivity reach for different detector configuraions. Just need to loop over different masses and use the previously introduced funtion get_events

    setup_dict={
##        "F-default":{
##            "name":"F - No cut",
##            "color":"maroon",
##            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.1",
##            "length":1.5,
##            "distance":480,
##            "lumi":140,
##            "channels": None,
##            "geo": None,
##            "separation":0
##        },
##        "F-default-test":{
##            "name":"F - test",
##            "color":"red",
##            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.1",
##            "length":1.5,
##            "distance":480,
##            "lumi":140,
##            "channels": None,
##            "geo": [(0.55,1.5),(0,0.3),(0.55,1),(0,0.3),(0.55,1)],
##            "separation":0
##        },
#        "F-default-sep":{
#            "name":"F - 300 um",
#            "color":"firebrick",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.1",
#            "length":1.5,
#            "distance":480,
#            "lumi":140,
#            "channels": None,
#            "geo":[(0.55,1.5),(0,0.3),(0.55,1),(0,0.3),(0.55,1)],
#            "separation":300e-6,
#            #"refsetup":"F-default-test"
##        },
        "TEST_1m0T":{
            "name":"1m 0T",
            "color":"firebrick",
            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.1",
            "length":1,
            "distance":480,
            "lumi":140,
            "channels": None,
            "geo":[(0.,1)],
            "separation":300e-6,
            #"refsetup":"F-default-test"
        },
#        "TEST_1m1T":{
#            "name":"1m 1T",
#            "color":"red",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.1",
#            "length":1,
#            "distance":480,
#            "lumi":140,
#            "channels": None,
#            "geo":[(1.,1)],
#            "separation":300e-6,
#            #"refsetup":"F-default-test"
#        },
#        "TEST_1p5m0p5T":{
#            "name":"1.5m 0.5T",
#            "color":"red",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.1",
#            "length":1.5,
#            "distance":480,
#            "lumi":140,
#            "channels": None,
#            "geo":[(0.5,1.5)],
#            "separation":300e-6,
#            #"refsetup":"F-default-test"
#        },
#        "F-default-sep-station0":{
#            "name":"F - Station1",
#            "color":"firebrick",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.1",
#            "length":1.5,
#            "distance":480,
#            "lumi":140,
#            "channels": None,
#            "geo":[(0.5,1.5)],
#            "separation":300e-6,
#            #"refsetup":"F-default-test"
#        },


##        "F-500um":{
##            "name":"F - 500 um",
##            "color":"coral",
##            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.1",
##            "length":1.5,
##            "distance":480,
##            "lumi":140,
##            "channels": None,
##            "geo":[(0.55,1.5),(0,0.3),(0.55,1),(0,0.3),(0.55,1)],
##            "separation":500e-6,
##            "refsetup":"F-default-test"
##        },
#        "F2_ee_default":{
#            "name":"F2 no cut",
#            "color":"maroon",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#            "length":5,
#            "distance":480,
#            "channels": ["e_e"],
#            "lumi":3000,
#            "geo": None,
#            "separation":0
#        },
#        
#        "F2_0p55T5m-0T20m_300um":{
#            "name":"F2 0.55T 5m",
#            "color":"firebrick",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#            "length":5,
#            "distance":480,
#            "channels": ["e_e"],
#            "lumi":3000,
#            "geo": [(0.55,5),(0,20)],
#            "separation":300e-6
#        },
#
#        "F2_0T25m-300um":{
#            "name":"F2 0T 25m",
#            "color":"coral",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#            "length":5,
#            "distance":480,
#            "channels": ["e_e"],
#            "lumi":3000,
#            "geo": [(0.,25)],
#            "separation":300e-6
#        },
#
#        "S1_L1p5-D2_0p5T1p5m":{
#            "name":"S1 0.5T 1.5m",
#            "color":"rebeccapurple",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#            "length":1.5,
#            "distance":497,
#            "channels": ["e_e"],
#            "lumi":3000,
#            "geo": [(0.5,1.5),(0,4.5)],
#            "separation":300e-6
#        },
#
#        "S1_L1p5-D2_1T1m":{
#            "name":"S1 1T 1m",
#            "color":"darkorchid",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#            "length":1.5,
#            "distance":497,
#            "channels": ["e_e"],
#            "lumi":3000,
#            "geo": [(1,1),(0,5)],
#            "separation":300e-6
#        },
#        
#        "S2_L2-D2_0p5T2m":{
#            "name":"S2 0.5T 2m",
#            "color":"darkorange",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#            "length":2,
#            "distance":500,
#            "channels": ["e_e"],
#            "lumi":3000,
#            "geo": [(0.5,2),(0,5)],
#            "separation":300e-6
#        },
#
#        "S2_L1p5-D2_1T1m":{
#            "name":"S2 1T 2m",
#            "color":"orange",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#            "length":2,
#            "distance":500,
#            "channels": ["e_e"],
#            "lumi":3000,
#            "geo": [(1,2),(0,5)],
#            "separation":300e-6
#        },
#
#        "S3_L1p5-D2_0p5T10m":{
#            "name":"S3 0.5T 10m",
#            "color":"darkgreen",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#            "length":10,
#            "distance":615,
#            "channels": ["e_e"],
#            "lumi":3000,
#            "geo": [(0.5,10),(0,15)],
#            "separation":300e-6
#        },
#
#        "S3_L1p5-D2_1T10m":{
#            "name":"S3 1T 10m",
#            "color":"forestgreen",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#            "length":10,
#            "distance":615,
#            "channels": ["e_e"],
#            "lumi":3000,
#            "geo": [(1,10),(0,15)],
#            "separation":300e-6
#        },


    }



    #################
# Plot separations
if run_plotseparations:    
    print("INFO   : Plot particle separations")
    
    ## Look at separation distributions of the dark photons which decay inside the detector
    fig = plt.figure(figsize=(7,5))
    ax = plt.subplot(1,1,1)

    maxy=0
    for setup in setup_dict:
        if scan_search and not sum([1 if s in setup else 0 for s in scan_search]):
            continue
        
        if not setup_dict[setup]["geo"]:
            continue
        
        print(f"INFO   :     - Running config: {setup}")
        
        setup, selection, length, channels, distance, luminosity, geo = setup, setup_dict[setup]["selection"], setup_dict[setup]["length"], setup_dict[setup]["channels"], setup_dict[setup]["distance"], setup_dict[setup]["lumi"], setup_dict[setup]["geo"]
        foresee.set_detector(selection=selection, channels=channels, length=length, distance=distance, luminosity=luminosity, geo=geo)

        #### Get reach 
        list_nevents = []
        list_seps = []
        #for mass in [0.01,0.1,1.0]:
        for mass in [0.1]:
            print(f"INFO   :       - Mass {mass}")

            coups=[1e-6]
            linestyle="solid"
            if mass == 0.1:
                coups=[1e-5]
                linestyle="dashed"
            if mass == 0.01:
                coups=[1e-4]
                linestyle="dotted"
                
            couplings, _, nevents, _, weights, seps, _ = foresee.get_events(mass=mass, energy=energy, couplings=coups)

            for sep,weight in zip(seps,weights):
                tmp_s=[]
                tmp_w=[]
                #sep_cut=300.e-6
                sep_cut=0 # no cut
                for n,s in enumerate(sep):
                    tmp_s.extend([e for e in s if e>sep_cut])
                    tmp_w.extend([weight[n] for e in s if e>sep_cut])
                
                y, x, _ = ax.hist(tmp_s, weights=tmp_w/sum(tmp_w), bins=np.logspace(-5,-1, 40), histtype='step', label=f'{setup_dict[setup]["name"]}, m={mass} GeV', linestyle=linestyle,color=setup_dict[setup]["color"])#, density=True) 
                ax.set_xscale("log")
                ax.set_yscale("linear")
                if max(y)*1.2 > maxy:
                    maxy=max(y)*1.2
                ax.set_ylim(0,maxy)
                #ax.set_ylim(0,1)
                #ax.set_yscale("log")
                #ax.set_ylim(10**-7,10**5) 
                ax.set_xlabel("Separation [m]") 
                ax.set_ylabel("Number of Events per bin") 
                ax.legend(frameon=False, labelspacing=0)
    plt.tight_layout()
    #plt.show()
    plt.savefig("NewConfigs_v3-DarkPhoton-EPOSLHC_LLP_m100-Separation%s.pdf"%("_"+scan_search[0] if scan_search else ""))




if run_runscans:
    print("INFO   :   - Run each config")

    done_scans=[]
    for setup in setup_dict:
        print(f"INFO   :     - Running config: {setup}")
        #### Check if run already processed
        
        if "refsetup" in setup_dict[setup]:
            refsetup=setup_dict[setup]["refsetup"]
            print(f"Using reference setup: {refsetup}. Skipping production")
            continue

        outfile="files/models/"+modelname+"/results/"+energy+"TeV_"+setup+".npy"
        #### Specify setup 
        
        setup, selection, length, channels, distance, luminosity, geo = setup, setup_dict[setup]["selection"], setup_dict[setup]["length"], setup_dict[setup]["channels"], setup_dict[setup]["distance"], setup_dict[setup]["lumi"], setup_dict[setup]["geo"]
        foresee.set_detector(selection=selection, channels=channels, length=length, distance=distance, luminosity=luminosity, geo=geo)

        
        #### Get reach 
        list_nevents = []
        list_seps = []
        for mass in masses:
            print(f"INFO   :       - Mass {mass}")
            couplings, _, nevents, _, _, seps, _ = foresee.get_events(mass=mass, energy=energy)
            #if setup_dict[setup]["separation"]:
            #    for n,nevent in enumerate(nevents):
                    
            list_nevents.append(nevents)
            list_seps.append(seps)
        
            
        #### Save results
        np.save(outfile,[masses,couplings,list_nevents,list_seps])
    
    

##########
# Plot the Results - Configurations
if run_plotreach:
    
    print("INFO   :   - Plotting reach")
    
    #Now let's plot the results. We first specify all detector setups for which we want to show result (filename in model/results directory, label, color, linestyle, opacity alpha for filled contours, required number of events).
    
    setups=[]
    if not scan_search:
        setups = [ ["14TeV_%s.npy"%(setup_dict[setup]["refsetup"] if "refsetup" in setup_dict[setup] else setup), setup_dict[setup]["name"],setup_dict[setup]["color"], "solid", 0., 3, setup_dict[setup]["separation"]] for setup in setup_dict ]
    else:
        setups = [ ["14TeV_%s.npy"%(setup_dict[setup]["refsetup"] if "refsetup" in setup_dict[setup] else setup), setup_dict[setup]["name"],setup_dict[setup]["color"], "solid", 0., 3, setup_dict[setup]["separation"]] for setup in setup_dict if scan_search in setup]

    
    print(f"INFO   :     - Found {len(setups)} setups")
    #Then we specify all the existing bounds (filename in model/bounds directory, label, label position x, label position y, label rotation)
    
    bounds = [ 
        ["bounds_LHCb1.txt", "LHCb",  0.220, 2.8*10**-4, 0  ],
        ["bounds_LHCb2.txt",  None  , 0    , 0         , 0  ],
        ["bounds_LHCb3.txt",  None  , 0    , 0         , 0  ],
        ["bounds_E137.txt",  "E137",  0.015, 1.2*10**-7, 0  ],
        ["bounds_CHARM.txt", "CHARM", 0.120, 1.4*10**-7, -5 ],
        ["bounds_NuCal.txt", "NuCal", 0.042, 1.3*10**-5, -28],
        ["bounds_E141.txt",  "E141",  0.011, 5.8*10**-5, 33 ],
        ["bounds_NA64.txt",  "NA64",  0.014, 3.6*10**-4, -35],
        ["bounds_BaBar.txt", "BaBar", 0.360, 1.4*10**-3, 0  ],
        ["bounds_NA48.txt",  "NA48",  0.040, 1.4*10**-3, 0  ],
    ]
    
    
    #We then specify other projected sensitivitities (filename in model/bounds directory, color, label, label position x, label position y, label rotation)
    
    projections = [
        ["limits_SeaQuest.txt",   "lime",         "SeaQuest", 1.350, 2.1*10**-7, 0  ],
        ["limits_NA62.txt",       "limegreen",    "NA62"    , 0.999, 1.1*10**-7, 0  ],
        ["limits_SHiP.txt",       "forestgreen",  "SHiP"    , 1.750, 8.2*10**-7, 0  ],
        ["limits_HPS.txt",        "deepskyblue",  "HPS"     , 0.050, 1.5*10**-4, 0  ],
        ["limits_HPS-1.txt",      "deepskyblue",  None      , 0    , 0         , 0  ],
        ["limits_Belle2.txt",     "blue",         "Belle2"  , 0.570, 1.3*10**-4, 0  ],
        ["limits_LHCb.txt",       "dodgerblue",   "LHCb"    , 0.135, 2.8*10**-4, 0  ],
        ["limits_LHCb-mumu1.txt", "dodgerblue",   None      , 0    , 0         , 0  ],
        ["limits_LHCb-mumu2.txt", "dodgerblue",   None      , 0    , 0         , 0  ],
    ]
    
    # Finally, we can plot everything using foresee.plot_reach(). It returns a matplotlib instance, to which we can add further lines and which we can show or save. Below, we add the dark matter relict target line for a specific benchmark.
    plot = foresee.plot_reach(
        setups=setups,
        bounds=bounds,
        projections=projections,
        title="Dark Photons", 
        xlims = [0.01,3], 
        ylims=[10**-7,0.002],
        xlabel=r"Dark Photon Mass $m_{A'}$ [GeV]", 
        ylabel=r"Kinetic Mixing $\epsilon$",
        legendloc=(1,0.68),
        figsize=(8,6),
    )
    
    data = foresee.readfile("files/models/"+modelname+"/lines/scalar_DM_Oh2_intermediate_eps_vs_mAprime.txt")
    plot.plot(data.T[0], data.T[1], color="k", lw=2)
    plot.text(0.010, 3.40*10**-5, "relic target",  fontsize=13,color="k",rotation=25)
    plot.text(0.011, 2.15*10**-5, r"$m_\chi\!=\!0.6 m_{A'}$",fontsize=13,color="k",rotation=25)
    plot.text(0.013, 1.20*10**-5, r"$\alpha_D\!=\!0.6$",fontsize=13,color="k",rotation=25)
    
    plot.subplots_adjust(left=0.12, right=0.97, bottom=0.10, top=0.95)
    plot.savefig("NewConfigs_v3-DarkPhoton-EPOSLHC-Reach%s.pdf"%("_"+scan_search[0] if scan_search else ""))
    #plot.show()




