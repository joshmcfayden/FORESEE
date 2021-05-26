import numpy as np
from src.foresee import Foresee, Utility, Model
import os



run_plothadrons=False
run_setupmodel=True
run_LLPspectra=False
run_rateexample=False
run_setupscans=True
run_runscans=False
run_plotreach=True

    
#############
# Initialization
print("INFO   : Initialise FORESEE")
foresee = Foresee()


#############
# Plot hadrons
if run_plothadrons:
    print("INFO   : Plot hadrons")
    plot=foresee.get_spectrumplot(pid="511", generator="Pythia8", energy="14")
    plot.savefig("NewConfigs_v1-DarkHiggs-Pythia8_Hadron-Angle_vs_Momentum.pdf")

    
#############
# Specifying the Model: Dark Higgss

if run_setupmodel:    
    print("INFO   : Setting up Dark Higgs model")
    energy = "14"
    modelname = "DarkHiggs"
    model = Model(modelname)

    
    print("INFO   :   - Adding production modes")
    
    ## 2-body decays
    model.add_production_2bodydecay(
    pid0 = "5",
    pid1 = "321",
    br = "5.7 * coupling**2 * pow(1.-pow(mass/5.279,2),2)",
    generator = "Pythia8",
    energy = energy,
    nsample = 10,
    )

    model.add_production_2bodydecay(
        pid0 = "-5",
        pid1 = "321",
        br = "5.7 * coupling**2 * pow(1.-pow(mass/5.279,2),2)",
        generator = "Pythia8",
        energy = energy,
        nsample = 10,
    )
    
    model.add_production_2bodydecay(
        pid0 = "25",
        pid1 = "0",
        br = "2*0.05",
        generator = "Pythia8",
        energy = energy,
        nsample = 100,
        scaling = 0,
    )

    # 3-body

    model.add_production_3bodydecay(
        label= "5_di",
        pid0 = "5",
        pid1 = "321",
        pid2 = "0",
        br = "7.37e-10*np.sqrt(1-4*mass**2/q**2)*(1-q**2/4.5**2)**2",
        generator = "Pythia8",
        energy = energy,
        nsample = 10,
        scaling = 0, 
    )
    
    model.add_production_3bodydecay(
        label= "-5_di",
        pid0 = "-5",
        pid1 = "321",
        pid2 = "0",
        br = "7.37e-10*np.sqrt(1-4*mass**2/q**2)*(1-q**2/4.5**2)**2",
        generator = "Pythia8",
        energy = energy,
        nsample = 10,
        scaling = 0, 
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
    0.1   ,  0.1122,  0.1259,  0.1413,  0.1585,  0.1778,  0.1995,  
    0.2239,  0.2512,  0.2818,  0.3162,  0.3548,  0.3981,  0.4467,  
    0.5012,  0.5623,  0.6026,  0.631 ,  0.6457,  0.6607,  0.6761,  
    0.6918,  0.7079,  0.7244,  0.7413,  0.7586,  0.7762,  0.7943,  
    0.8128,  0.8318,  0.8511,  0.871 ,  0.8913,  0.912 ,  0.9333,  
    0.955 ,  0.9772,  1.    ,  1.122 ,  1.2589,  1.4125,  1.5   ,
    1.5849,  1.7783,  1.9953,  2.2387,  2.5119,  2.8184,  3.1623,  
    3.5   ,  3.7   ,  3.9811,  5.0119,  6.3096,  7.9433,  10.   ,
    11.22 ,  12.589,  14.125,  15.849,  17.783,  19.953,  22.387,  
    25.119,  28.184,  31.623,  39.811,  50.119,  55.000,  60.000,
    63.096,  79.430,  99.9
]
    
##########
# Generate LLP Spectra
if run_LLPspectra:    
    print("INFO   : Generating LLP Spectra")
    

    ## Look at benchmark scenario with m_phi=1.5 GeV and theta=10^-4
    print("INFO   :   - Generating for m_phi=1.5 GeV and theta=10^-4")
    plt = foresee.get_llp_spectrum(mass=1.5, coupling=1e-4, do_plot=True, print_stats=True)
    plt.savefig("NewConfigs_v1-DarkHiggs-Pythia8_LLP_m1500_t10m4-Angle_vs_Momentum.pdf")
    
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
    
    
    print("INFO   :   - For dark higgs (m_phi=1.5 GeV) check #events in decay volume")
    ## For one dark higgs (m_phi=1.5 GeV) look at how many particles decay inside the decay volume.
    output = foresee.get_events(mass=1.5, energy=energy, couplings = np.logspace(-8,-3,6))
    coups, ctaus, nsigs, energies, weights, _ = output
    for coup,ctau,nsig in zip(coups, ctaus, nsigs):
        print ("epsilon =", '{:5.3e}'.format(coup), ": nsig =", '{:5.3e}'.format(nsig))
    
    print("INFO   :   - Plot energy distribution for different couplings")
    ## Look at energy distribution of the dark higgs which decay inside the detector

    fig = plt.figure(figsize=(7,5))
    ax = plt.subplot(1,1,1)
    for coup,en,weight in zip(coups,energies,weights):
        if sum(weight)<10**-5 : continue
        ax.hist(en, weights=weight, bins=np.logspace(2,4, 20), histtype='step', label=r"$\theta=$"+str(coup)) 
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(10**-7,10**5) 
        ax.set_xlabel("E($\phi$) [GeV]") 
        ax.set_ylabel("Number of Events per bin") 
        ax.legend(frameon=False, labelspacing=0)
    plt.tight_layout()
    plt.savefig("NewConfigs_v1-DarkHiggs-Pythia8_LLP_m1500-Energy.pdf")




##########
# Parameter Scan
setup_dict={}
if run_setupscans:
    print("INFO   : Setup parameter scan for different detector configs")
    
    ## Get the LLP sensitivity reach for different detector configuraions. Just need to loop over different masses and use the previously introduced funtion get_events
    
    
    setup_dict={
        "default":{
            "name":"L=5m D=2m (Default)",
            "color":"red",
            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
            "length":5,
            "distance":480,
            "channels": None
        },
        "L5-D1":{
            "name":"L=5m D=1m",
            "color":"firebrick",
            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.5",
            "length":5,
            "distance":480,
            "channels": None
        },
        "L5-D0p5":{
            "name":"L=5m D=0.5m",
            "color":"maroon",
            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.25",
            "length":5,
            "distance":480,
            "channels": None
        },
        "L15-D2":{
            "name":"L=15m D=2m",
            "color":"mediumorchid",
            "selection":"np.sqrt(x.x**2 + x.y**2)< 1",
            "length":15,
            "distance":480,
            "channels": None
        },
        "L15-D1":{
            "name":"L=15m D=1m",
            "color":"darkorchid",
            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.5",
            "length":15,
            "distance":480,
            "channels": None
        },
        "L15-D0p5":{
            "name":"L=15m D=0.5m",
            "color":"rebeccapurple",
            "selection":"np.sqrt(x.x**2 + x.y**2)< 0.25",
            "length":15,
            "distance":480,
            "channels": None
        },
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
            couplings, _, nevents, _, _ , _ = foresee.get_events(mass=mass, energy=energy, modes=["5","-5"], couplings = np.logspace(-9,-2,71))
            list_nevents.append(nevents)  
            
        #### Save results
        np.save(outfile,[masses,couplings,list_nevents])
    
    

##########
# Plot the Results - Configurations
if run_plotreach:
    
    print("INFO   :   - Plotting reach")
    
    #Now let's plot the results. We first specify all detector setups for which we want to show result (filename in model/results directory, label, color, linestyle, opacity alpha for filled contours, required number of events).
    
    setups = [ ["14TeV_%s.npy"%setup, setup_dict[setup]["name"],setup_dict[setup]["color"], "solid", 0., 3] for setup in setup_dict ]
    print(f"INFO   :     - Found {len(setups)} setups")
    #Then we specify all the existing bounds (filename in model/bounds directory, label, label position x, label position y, label rotation)
    
    bounds = [ 
        ["bounds_1508.04094.txt", "LHCb $B^0$"  , 0.430, 2.2*10**-3, 90 ],
        ["bounds_1612.08718.txt", "LHCb $B^+$"  , 0.330, 2.2*10**-3, 90 ],
        ["bounds_1612.08718.txt", "LHCb $B^+$"  , 2.500, 2.2*10**-3, 90 ],
        ["bounds_LSND.txt"      , "LSND"        , 0.250, 9.0*10**-5, 90 ],
        ["bounds_CHARM.txt"     , "CHARM"       , 0.250, 4.0*10**-4, 90 ],
        ["bounds_MicroBoone.txt", "$\mu$BooNE"  , 0.138, 3.4*10**-4, 90 ],
        ["bounds_E949.txt"      , "E949"        , 0.102, 1.5*10**-4, 90 ],
        ["bounds_2011.11329.txt", "NA62 $K^+$"  , 0.170, 6.2*10**-4, 90 ],
        ["bounds_2010.07644.txt", "NA62 $\pi^+$", 0.125, 2.4*10**-3, 90 ],
    ]
    
    
    #We then specify other projected sensitivitities (filename in model/bounds directory, color, label, label position x, label position y, label rotation)
    
    projections = [
        ["limits_SHiP.txt",       "teal",         "SHiP"    , 2.700, 3.2*10**-5, 0  ],
        ["limits_MATHUSLA.txt",   "dodgerblue",   "MATHUSLA", 0.120, 5.0*10**-6, 0  ],
        ["limits_CodexB.txt",     "deepskyblue",  "CodexB"  , 1.700, 2.0*10**-5, 0  ],
        ["limits_LHCb.txt",       "cyan",         "LHCb"    , 3.800, 1.0*10**-4, 0  ],
    ]
    
    # Finally, we can plot everything using foresee.plot_reach(). It returns a matplotlib instance, to which we can add further lines and which we can show or save. Below, we add the dark matter relict target line for a specific benchmark.
    plot = foresee.plot_reach(
        setups=setups,
        bounds=bounds,
        projections=projections,
        title="Dark Higgs", 
        xlims=[0.1,10], 
        ylims=[10**-5.5,10**-2.5],
        xlabel=r"Dark Higgs Mass $m_{\phi}$ [GeV]", 
        ylabel=r"Mixing $\theta$",
        legendloc=(1.00,0.95),
        figsize=(8,6),
    )

    plot.subplots_adjust(left=0.12, right=0.97, bottom=0.10, top=0.95)
    plot.savefig("NewConfigs_v1-DarkHiggs-Pythia8-Reach.pdf")





