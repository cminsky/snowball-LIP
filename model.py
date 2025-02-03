from dependencies import *

### the main model

# evolves bedrock and regolith height and weathering rate
def erosion_conditions(dt, # time step [Myr]
        B, # bedrock depth [m]
        H, # regolith depth [m]
        P, # sap. prod. rate [m/Myr]
        E # erosion rate [m/Myr]
        ):
    # chemical weathering
    w_max = P*dt # max sap. produced [m]
    w = np.min([w_max,B]) # can't weather more than B
    H += w
    B -= w

    # physical erosion
    dz = E*dt # max total rock eroded [m]
    H -= dz
    if H < 0: # eroded through regolith
        B += H

    H = np.max([0,H])
    B = np.max([0,B])

    return w,H,B

# carbon cycle model with LIP
def run_model(
        ## Model characteristics
        dt, # time step [Myr]
        t_max,  # time to run until [Myr]
        ## LIP emplacement characteristics
        emp_dur, # length of emplacement [Myr]
        A0, # LIP area [m2]
        B0, # total extrusive LIP height [m]
        erup_freq, # eruption frequency [Myr]
        degass, # [examol CO2]
        ## LIP weathering characteristics
        P0, # saprolite production rate [m/Myr]
        #E0, # erosion rate [m/Myr]
        E_P, # E0/P0 ratio
        d, # weathering attenuation height [m]
        c, # erosion dependence on slope []
        Xm, # Ca Mg fraction [mol/m3]
        ## Background climate
        N0, # initial surficial carbon inventory [examol]
        V, # global volcanic degassing flux [examol/Myr]
        ## Climate sensitivities
        n, # global weathering feedback strength []
        n_p, # LIP weathering feedback strength []
        n_e, # LIP erosion feedback strength []
        # Output
        prognostics = False
        ):

    ## Emplacement characteristic calculations
    erup_num = emp_dur/erup_freq # number of eruptions 
    B_e = B0/erup_num # height extruded each eruption [m]
    degass_e = degass/erup_num # CO2 degassed each eruption [examol CO2]

    ## Erosion
    E0 = E_P*P0

    ## Set up model
    t = [0] # timekeeping array [Myr]
    N = [N0] # surficial carbon array [examol]
    last_erup = 0 # time of last eruption [Myr]
    
    ## Prognostic variables
    B = [B_e] # bedrock [m]
    H = [0] # regolith [m]
    P = [P0]
    E = [E0]
    degass_arr = [degass_e]

    while t[-1] < t_max:
        # global weathering
        with np.errstate(over='ignore', invalid='ignore'):
            y = (N[-1]/N0)**2 # normalized atmospheric pCO2
        W = V*(y**n) # global silicate weathering [examol/Myr]
            
        # if an eruption occurs
        if t[-1]-last_erup > erup_freq and t[-1] < emp_dur:
            B[-1] += B_e # erupt bedrock
            H[-1] = 0 # reset regolith
            degass_arr[-1] = degass_e # degass
            last_erup = t[-1] # keep track of eruption time
        
        # LIP saprolite production
        P_i = P0*np.exp(-H[-1]/d) # soil production function
        with np.errstate(invalid='ignore'):
            P_i *= y**n_p # climate dependence

        # LIP erosion
        E_i = E0*(B[-1]/B_e)**c # dependence on slope
        with np.errstate(invalid='ignore'):
            E_i *= y**n_e # climate dependence

        # LIP weathering
        w_i,H_i,B_i = erosion_conditions(dt=dt,B=B[-1],H=H[-1],P=P_i,E=E_i)
        dC_i = w_i*A0*Xm/1e18 # m m2 mol/m3 --> examol
        
        # increment surficial carbon
        N_i = N[-1] + V*dt - W*dt - dC_i + degass_arr[-1]
        
        # main outputs
        t.append(t[-1]+dt)
        N.append(N_i)

        # prognostic variables
        B.append(B_i)
        H.append(H_i)
        P.append(P_i)
        E.append(E_i)
        degass_arr.append(0) # don't degass if no eruption

    if prognostics:
        return [np.array(t),np.array(N),
                np.array(B),np.array(H),
                np.array(P),np.array(E),
                np.array(degass_arr)]
    else:
        return np.array(t),np.array(N)

# temperature response to CO2
def x(N,N0,b,a):
    with np.errstate(over='ignore'):
        y = (N/N0)**2 # normalized pCO2
    x = b/a*np.log(y) # radiative forcing
    return x
