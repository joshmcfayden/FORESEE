import numpy as np
from src.foresee import Foresee, Utility, Model
import os
import random, math

def unweight_events(energies, thetas, weights, number):
    
    #initialize arrays and random number generator
    random.seed()
    unweighted_thetas=[]
    unweighted_energies=[]
    unweighted_weights=[]
    total_weight = np.sum(weights)
    event_weight = total_weight/float(number)

    print("Total weight",total_weight)
    #unweighting
    for irand in range (number):
        stopweight = random.random()*total_weight
        partid, weightsum = 0, 0
        
        #while weightsum < stopweight:
        #    weightsum+=weights[partid]
        #    partid+=1

        while True:
            weightsum+=weights[partid]
            if weightsum >= stopweight:
                break
            partid+=1
        
            
        #print("JOSH",partid,len(thetas))
        unweighted_thetas.append(thetas[partid])
        unweighted_energies.append(energies[partid])
        unweighted_weights.append(event_weight)
        
    return  np.array(unweighted_energies), np.array(unweighted_thetas), np.array(unweighted_weights)


run_plotpions=False
run_setupmodel=True
run_LLPspectra=False
run_rateexample=False
run_setupscans=True
run_runscans=False
run_plotreach=True

scan_search=None
scan_search={"F2":["F2"]}
scan_search={"S1":["S1"]}
scan_search={"S2":["S2"]}
scan_search={"S3":["S3"]}

scan_name="pick"
scan_search={"pick":["F2-default","S1-L1p5-D2","S2-L2-D2","S3-L10-D1","S3-L10-D2"]}
scan_search={"pick":["F2-default","S1-L1p5-D1","S1-L1p5-D2","S2-L10-D1","S2-L10-D2"]}
scan_search=None

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
    plot.savefig("NewConfigs_v2-GenHepMC-DarkPhoton-EPOSLHC_Pion-Angle_vs_Momentum.pdf")

    
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

masses = [ 0.01  , 0.05, 0.1,  1.,  10. ]

masses = [ 
    0.01  ,
    0.0501,
    0.1585,
    0.3548,
    0.6457,
    0.7586,
    0.8913,
    1.2589,
    2.8184
]

##########
# Generate LLP Spectra
if run_LLPspectra:    
    print("INFO   : Generating LLP Spectra")
    

    
    print("INFO   :   - Generating for m_A=100 MeV and epsilon=10^-5")
    ## Look at benchmark scenario with m_A=100 MeV and epsilon=10^-5
    plt = foresee.get_llp_spectrum(0.1, coupling=10**(-5), do_plot=True)
    #plt.show()
    plt.savefig("NewConfigs_v2-GenHepMC-DarkPhoton-EPOSLHC_LLP_m100_e10m5-Angle_vs_Momentum.pdf")
    
    print("INFO   :   - Generating for all masses")
    for mass in masses:
        foresee.get_llp_spectrum(mass=mass,coupling=1)
    
    
    
#################
# Count Event Rate in Detector
if run_rateexample:    
    print("INFO   : Count Event Rate in Detector")
    
    ## To count the #decays in detector volume need detector geometry
    ## These are FASER2 defaults
    distance, selection, length, luminosity, channels = 480, "np.sqrt(x.x**2 + x.y**2)< 1", 5, 3000, None
    foresee.set_detector(distance=distance, selection=selection, length=length, luminosity=luminosity, channels=channels)
    
    
    print("INFO   :   - For dark photon (m_A'=100 MeV) check #events in decay volume")
    ## For one dark photon (m_A'=100 MeV) look at how many particles decay inside the decay volume.
    mass = 0.1 
    output = foresee.get_events(mass=0.1, energy=energy, couplings=np.logspace(-8,-3,6), )
    coups, ctaus, nsigs, energies, weights, _ = output
    for coup,ctau,nsig in zip(coups, ctaus, nsigs):
        print ("epsilon =", '{:5.3e}'.format(coup), ": nsignal =", '{:5.3e}'.format(nsig))
    
    print("INFO   :   - Plot energy distribution for different couplings")
    ## Look at energy distribution of the dark photons which decay inside the detector
    fig = plt.figure(figsize=(7,5))
    ax = plt.subplot(1,1,1)
    for coup,en,weight in zip(coups,energies,weights):
        if sum(weight)<10**-5 : continue
        ax.hist(en, weights=weight, bins=np.logspace(2,4, 20), histtype='step', label=r"$\epsilon=$"+str(coup)) 
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(10**-8,10**5) 
        ax.set_xlabel("E(A') [GeV]") 
        ax.set_ylabel("Number of Events per bin") 
        ax.legend(frameon=False, labelspacing=0)
    plt.tight_layout()
    #plt.show()
    plt.savefig("NewConfigs_v2-GenHepMC-DarkPhoton-EPOSLHC_LLP_m100-Energy.pdf")




##########
# Parameter Scan
setup_dict={}
if run_setupscans:
    print("INFO   : Setup parameter scan for different detector configs")
    
    ## Get the LLP sensitivity reach for different detector configuraions. Just need to loop over different masses and use the previously introduced funtion get_events

    setup_dict={
#        "F2-default":{
#            "name":"F2 L=5m D=2m",# (Default)",
#            #"color":"maroon",
#            "color":"firebrick",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#            "length":5,
#            "distance":480,
#            "channels": None
#        },
#        "F2-L5-D1":{
#            "name":"F2 L=5m D=1m",
#            "color":"firebrick",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.5",
#            "length":5,
#            "distance":480,
#            "channels": None
#        },
#        "F2-L5-D0p5":{
#            "name":"F2 L=5m D=0.5m",
#            "color":"indianred",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.25",
#            "length":5,
#            "distance":480,
#            "channels": None
#        },
#        "S1-L1p5-D2":{
#            "name":"S1 L=1.5m D=2m",
#            "color":"rebeccapurple",
#            #"color":"darkorchid",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#            "length":1.5,
#            "distance":497,
#            "channels": None
#        },
#        "S1-L1p5-D1":{
#            "name":"S1 L=1.5m D=1m",
#            "color":"darkorchid",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.5",
#            "length":1.5,
#            "distance":497,
#            "channels": None
#        },
#        "S1-L1p5-D0p5":{
#            "name":"S1 L=1.5m D=0.5m",
#            "color":"mediumorchid",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.25",
#            "length":1.5,
#            "distance":497,
#            "channels": None
#        },
#        "S2-L2-D2":{
#            "name":"S2 L=2m D=2m",
#            "color":"darkorange",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
#            "length":2,
#            "distance":500,
#            "channels": None
#        },
#        "S2-L2-D1":{
#            "name":"S2 L=2m D=1m",
#            "color":"orange",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.5",
#            "length":2,
#            "distance":500,
#            "channels": None
#        },
#        "S2-L2-D0p5":{
#            "name":"S2 L=2m D=0.5m",
#            "color":"navajowhite",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.25",
#            "length":2,
#            "distance":500,
#            "channels": None
#        },
        "S2-L10-D2":{
            "name":"S2 L=10m D=2m",
            ##"color":"royalblue",
            #"color":"darkgreen",
            "color":"darkgray",
            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
            "length":10,
            "distance":615,
            "channels": None
        },
        "S2-L10-D2_sep1":{
            "name":"S2 L=10m D=2m (sep=1 mm)",
            ##"color":"royalblue",
            #"color":"darkgreen",
            "color":"darkred",
            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
            "length":10,
            "distance":615,
            "channels": None,
            "separation":1
        },
        "S2-L10-D2_sep10":{
            "name":"S2 L=10m D=2m (sep=10 mm)",
            ##"color":"royalblue",
            #"color":"darkgreen",
            "color":"red",
            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
            "length":10,
            "distance":615,
            "channels": None,
            "separation":10
        },
        "S2-L10-D2_sep50":{
            "name":"S2 L=10m D=2m (sep=50 mm)",
            ##"color":"royalblue",
            #"color":"darkgreen",
            "color":"orange",
            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
            "length":10,
            "distance":615,
            "channels": None,
            "separation":50
        },
        "S2-L10-D2_sep100":{
            "name":"S2 L=10m D=2m (sep=100 mm)",
            ##"color":"royalblue",
            #"color":"darkgreen",
            "color":"gold",
            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
            "length":10,
            "distance":615,
            "channels": None,
            "separation":100
        },
#        "S2-L10-D1":{
#            "name":"S2 L=10m D=1m",
#            ##"color":"cornflowerblue",
#            #"color":"forestgreen",
#            "color":"navajowhite",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.5",
#            "length":10,
#            "distance":615,
#            "channels": None
#        },
#        "S2-L10-D0p5":{
#            "name":"S2 L=10m D=0.5m",
#            #"color":"lightsteelblue",
#            "color":"limegreen",
#            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.25",
#            "length":10,
#            "distance":615,
#            "channels": None
#        },
    }


if run_runscans:
    print("INFO   :   - Run each config")
    
    for setup in setup_dict:
        print(f"INFO   :     - Running config: {setup}")
        #### Check if run already processed
        outfile="files/models/"+modelname+"/results/"+energy+"TeV_"+setup+".npy"
        
        #### Specify setup 
        luminosity = 3000
        setup, selection, length, channels, distance = setup, setup_dict[setup]["selection"], setup_dict[setup]["length"], setup_dict[setup]["channels"], setup_dict[setup]["distance"]
        foresee.set_detector(selection=selection, channels=channels, length=length, distance=distance, luminosity=luminosity)

        
        #### Get reach 
        list_nevents = []    
        for mass in masses:
            #couplings, _, nevents, _, _ , _, _ = foresee.get_events(mass=mass, energy=energy, couplings=np.logspace(-8,-3,12))
            couplings, ctaus, nevents, energies, weights, seps, thetas = foresee.get_events(mass=mass, energy=energy, couplings=np.logspace(-8,-3,12),nsample=100 )

            list_nevents.append(nevents)  
            print(len(couplings),len(thetas))
            
            # Make HepMC files
            for n,coup in enumerate(couplings):
                massname=("%.4f"%mass).replace('.','p')
                coupname='{:5.3e}'.format(coup).replace('.','p')
                
                outname=f'output_{setup}_mass{massname}GeV_coup{coupname}.hepmc'
                print(outname)
                
                weighted_energies, weighted_thetas, weighted_weights = energies[n], thetas[n], weights[n]
                unweighted_energies, unweighted_thetas, unweighted_weights = unweight_events (energies[n], thetas[n], weights[n], 10000)#int(len(energies[n])/10))
                
                print(len(energies),len(unweighted_energies))
        
                datax=[]
                datay=[]
                dataz=[]
                datae=[]
                dataw=[]
                #print "INFO:   - Unweighted energies":
                #print unweighted_energies
                
                # get events
                f= open(outname,"w")
                f.write("\nHepMC::Version 2.06.09\n")
                f.write("HepMC::IO_GenEvent-START_EVENT_LISTING\n")
                for i,(itheta,ienergy,iweight) in enumerate(zip(unweighted_thetas, unweighted_energies, unweighted_weights)):
            
                    #Event Info
                    f.write("E %s -1 -1.0000000000000000e+00 -1.0000000000000000e+00 -1.0000000000000000e+00 20 0 1 1 2 0 0\n"%(i+1))
                    #f.write("N 1 \"0\" \n")
                    f.write("U GEV MM\n")
                    #f.write("C "+str(iweight)+" 0 \n")
                    
                    #vertex
                    posz = random.uniform(0,length)*1000.
                    phi  = random.uniform(-math.pi,math.pi)
                    posy = itheta*480.*1000*np.sin(phi)
                    posx = itheta*480.*1000*np.cos(phi) 
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
                    datae.append(ienergy)
                    dataw.append(iweight)
            
                    #particles
                    particles = foresee.decay_llp(mass=mass, energy=ienergy)
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
 

            
        #### Save results
        np.save(outfile,[masses,couplings,list_nevents])
    
    

##########
# Plot the Results - Configurations
if run_plotreach:
    
    print("INFO   :   - Plotting reach")
    
    #Now let's plot the results. We first specify all detector setups for which we want to show result (filename in model/results directory, label, color, linestyle, opacity alpha for filled contours, required number of events).
    
    setups=[]
    for setup in setup_dict:
        if not scan_search or setup in scan_search[scan_name]:
            setups.append(["14TeV_%s.npy"%setup, setup_dict[setup]["name"],setup_dict[setup]["color"], "solid", 0., 3, setup_dict[setup]["separation"] if "separation" in setup_dict[setup] else None])

    
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
    projections = []
    #projections = [
    #    ["limits_SeaQuest.txt",   "lime",         "SeaQuest", 1.350, 2.1*10**-7, 0  ],
    #    ["limits_NA62.txt",       "limegreen",    "NA62"    , 0.999, 1.1*10**-7, 0  ],
    #    ["limits_SHiP.txt",       "forestgreen",  "SHiP"    , 1.750, 8.2*10**-7, 0  ],
    #    ["limits_HPS.txt",        "deepskyblue",  "HPS"     , 0.050, 1.5*10**-4, 0  ],
    #    ["limits_HPS-1.txt",      "deepskyblue",  None      , 0    , 0         , 0  ],
    #    ["limits_Belle2.txt",     "blue",         "Belle2"  , 0.570, 1.3*10**-4, 0  ],
    #    ["limits_LHCb.txt",       "dodgerblue",   "LHCb"    , 0.135, 2.8*10**-4, 0  ],
    #    ["limits_LHCb-mumu1.txt", "dodgerblue",   None      , 0    , 0         , 0  ],
    #    ["limits_LHCb-mumu2.txt", "dodgerblue",   None      , 0    , 0         , 0  ],
    #]
    
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
        legendloc=(1.02,0.72),
        figsize=(8,6),
        sepfile="/Users/mcfayden/Work/FASER/FASER2/G4_test/FASER2_HepMC_v4_FASER2_Cavern_3rdTrkStation-build/sep_effs.csv"
    )
    
    #data = foresee.readfile("files/models/"+modelname+"/lines/scalar_DM_Oh2_intermediate_eps_vs_mAprime.txt")
    #plot.plot(data.T[0], data.T[1], color="k", lw=2)
    #plot.text(0.010, 3.40*10**-5, "relic target",  fontsize=13,color="k",rotation=25)
    #plot.text(0.011, 2.15*10**-5, r"$m_\chi\!=\!0.6 m_{A'}$",fontsize=13,color="k",rotation=25)
    #plot.text(0.013, 1.20*10**-5, r"$\alpha_D\!=\!0.6$",fontsize=13,color="k",rotation=25)
    
    plot.subplots_adjust(left=0.12, right=0.97, bottom=0.10, top=0.95)
    plot.savefig("NewConfigs_v2-SepScan-DarkPhoton_F2Cavern3rdStation-EPOSLHC-Reach%s.pdf"%("_"+scan_name if scan_search else ""))
    #plot.show()





