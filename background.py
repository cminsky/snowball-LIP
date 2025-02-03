from dependencies import *

### generates PDFs from Krissansen-Totton et al. (2018) model runs

# function for interpolating between time steps
def interp(t,t_arr,data):
    tstep = t_arr[1]-t_arr[0]
    i1,i2 = np.where(np.abs(t_arr-t) < tstep)[0]
    times = t_arr[i1:i2+1] # t1,t2 = t_arr[i1],t_arr[i2]
    data_t = data[i1:i2+1]
    linfit = interp1d(times,data_t,axis=0)
    return linfit(t)

def background_fits(t): # t in Myrs from present
    # load data
    try:
        all_output = np.load('data/kt_data.npy')
    except:
        all_output = np.load('kt_data.npy')
    t_arr = all_output[4,:,0]/1e6 # time steps in Myrs

    # interpolate model runs at specified age
    CO2 = interp(t,t_arr,all_output[6,:,:])*1e6 # CO2 (ppm)
    T = interp(t,t_arr,all_output[17,:,:]) # surface temp (K)
    V = interp(t,t_arr,all_output[24,:,:])/1e12 # volcanic degassing (Tmol/yr to examol/Myr)

    # create fits
    CO2_fit = st.lognorm.fit(CO2) # log normal: s, loc, scale
    T_fit = st.norm.fit(T) # normal: loc, scale
    V_fit = st.norm.fit(V) # normal: loc, scale

    return CO2_fit,T_fit,V_fit

def background_ranges(t_Earth,iters,
                      N_pi = 2.83, # pre-industrial surficial carbon [examol]
                      ppm_pi = 280 # pre-industrial pCO2 [ppm]
                      ):
    CO2_fit,T_fit,V_fit = background_fits(t_Earth)
    CO2 = st.lognorm.rvs(*CO2_fit,size=iters) # pCO2 [ppm]
    N0 = N_pi*np.sqrt(CO2/ppm_pi) # surficial carbon [examol]
    V = st.norm.rvs(*V_fit,size=iters) # volcanic degassing [examol/Myr]
    T0 = st.norm.rvs(*T_fit,size=iters) # surface temp [K]
    return N0,V,T0

def y_t(t,t0=4.54e3,b=5.35,alb=0.3,F0=1366):
    L_normed = 1/(1+2/5*(1-(t0-t)/t0)) # Gough et al 1981
    F = F0*L_normed
    delS = (F-F0)*(1-alb)/4
    y = np.exp(-delS/b)
    return y

def N_t(t,b=5.35,alb=0.3,F0=1366,N0=2.83):
    n = np.sqrt(y_t(t,b=b,alb=alb,F0=F0))
    N = N0*n
    return N

# read in Foster et al. 2017 data and estimate pCO2
def Foster_CO2(t_Earth,iters,
               N_pi = 2.83, # pre-industrial surficial carbon [examol]
               ppm_pi = 280 # pre-industrial pCO2 [ppm]
                    ):
    df = pd.read_excel('data/foster-2017.xlsx',skiprows=1)
    row = df.iloc[(df["Age (Ma)"] - t_Earth).abs().idxmin()]
    mean, lw68, up68 = row["pCO2 probability maximum"],row["lw68%"], row["up68%"]
    std = (up68 - lw68) / (2 * st.norm.ppf(0.84))
    CO2_fit = st.norm(loc=mean, scale=std)
    CO2 = CO2_fit.rvs(size=iters) # list of pCO2 values [ppm]
    N0 = N_pi*np.sqrt(CO2/ppm_pi) # surficial carbon [examol]
    return N0

# read in Scotese et al. 2021 data and estimate temperature
def Scotese_T(t_Earth,iters):
    xls = pd.ExcelFile('data/scotese_2021.xlsx')
    df = pd.read_excel(xls, '1my version')
    row = df.iloc[(df["Age"] - t_Earth).abs().idxmin()]
    T = row["GAT"]+273.15 # C to K
    return [T]*iters