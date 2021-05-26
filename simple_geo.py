import math


def get_deflection(p,B,L,verbose=False):
    
    #if verbose:
    #    print(f"Momentum:   {p:.2f} GeV/c")
    #    print(f"B field:    {B:.2f} T")
    #    print(f"Length:     {L:.2f} m")


    # radius of curvature
    R = p / (0.3 * B)
    
    if verbose: print(f"Radius:     {R:.2f} m")
    
    # theta - angle showing how far through circular trajectory particle is
    th = math.asin( L / R )
    #if verbose: print(f"Theta:      {th*1000.:.2f} mrad")

    # displacement of particle in bending direction
    d = R * ( 1. - math.cos(th) )
        
    return th,d

def process_geo(p,m,a,geo,pos_z=None,verbose=False):
    tot_d=0
        
    # iterate through geometry
    for n,g in enumerate(geo):
    
        #if verbose:
        #    print(f"\n\nGeo Element #{n+1}")
        #    print("------------------------")
        #    print(f"alpha:      {a*1000:.2f} mrad")
        
        (B,L)=g
        
        if n==0 and pos_z:
            L=L-pos_z
        
        d=0
        if B:
            #calculate bending
            
            p_prime=p*math.cos(a)
            L_prime=L/math.cos(a)
            B_prime=B*math.cos(a)
            
            #th,d_prime,s_prime=get_deflection(p_prime,B_prime,L_prime)
            th,d_prime=get_deflection(p,B,L_prime)
            #if verbose: print(f"theta:      {th} rad")
            
            d_a=L*math.tan(a)       
            d_B=d_prime/math.cos(a)
    
            #if verbose:
            #    print(f"d_alpha:    {d_a*1e6:.2f} um")
            #    print(f"d_B:        {d_B*1e6:.2f} um")
            
            d=d_a+d_B
            a=a+th
            
        else:
            # just extrapolate straight line
            d=L*math.tan(a)
            
        tot_d=tot_d+d
        
        #if verbose:
        #    print(f"geo{n+1} def:   {d*1e6:.2f} um")
        #    print(f"geo{n+1} sep:   {2*d*1e6:.2f} um")

    return tot_d

def calc_separation(p,m,geo,pos_z=None,th=None):
    #if True or verbose:
    #    print(f"Momentum:   {p} GeV/c")
    #    print(f"Mass:       {m} GeV/c^2")
    #
    #    print("\nGeometry:")
    #    for n,g in enumerate(geo):
    #        print(f"#{n+1}   B={g[0]} T, L={g[1]} m")
        
    tot_d=0

    a=0.5*m/p
    if th:
        a+th
        
    d=process_geo(p,m,a,geo,pos_z=pos_z)

    (B,L)=geo[0]

    if pos_z:
        L=L-pos_z
    
    tot_d=2.*d
    
        
    #if verbose:
    #print("\n-----------------------")
    #print(f"Total deflection: {d*1e6:.2f} um")
    ##print(f"Incoming theta: {th*1e3:.2f} mrad")
    #print(f"Tofal separation: {tot_d*1e6:.2f} um")

    return tot_d
