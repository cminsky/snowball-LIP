from dependencies import *
from model import *
from background import *
from defaults import *

import streamlit as st
import scipy.stats as sst

st.set_page_config(layout="wide")

"""
###  Climate response to Large Igneous Province weathering

This carbon cycle model simulates the response of the global climate to the emplacement and subsequent silicate weathering of the Franklin Large Igneous Province (LIP).
The model is described in C. Minsky, R. Wordsworth, F. A. Macdonald, and A. H. Knoll (in review) "Neoproterozoic Snowball Earth Initiation from Silicate Weathering of a Large Igneous Province."

This GUI can be used to test which combinations of background climate and LIP characteristics cause initiation of a Snowball Earth transition.
Use the sliders in each tab to modify the model input parameters.
"""

col1,col2 = st.columns(2)

with col1:
    tab1,tab2,tab3 = st.tabs(["Background climate",
                            "LIP emplacement",
                            "Weathering and erosion"])

with tab1:
    t_Earth = st.number_input(label='Time in Earth history (Mya)',value=719,step=1)
    CO2_fit,T_fit,V_fit = background_fits(t_Earth)

    T0 = st.number_input(label='Initial temperature (K)',value=288,step=1)
    #T95 = sst.norm.interval(0.95,*T_fit)
    #T_default = sst.norm.median(*T_fit)
    #T0 = st.slider(label='Initial temperature (K)',value=T_default,
    #               min_value=int(np.max([T95[0],0])),max_value=int(T95[1]))

    CO295 = sst.lognorm.interval(0.95,*CO2_fit)
    CO2_default = sst.lognorm.median(*CO2_fit)
    CO2 = st.slider(label='pCO$_2$ (ppm)',value=int(CO2_default),
                    min_value=np.max([int(CO295[0]),0]),max_value=int(CO295[1]))
    N0 = N_pi*np.sqrt(CO2/ppm_pi)

    V95 = sst.norm.interval(0.95,*V_fit)
    V_default = sst.norm.median(*V_fit)
    V = st.slider(label='Global CO$_2$ degassing (10$^{12}$ mol/yr)',value=V_default,
                    min_value=np.max([V95[0],0]),max_value=V95[1])
    
    b = st.number_input(label='CO2 radiative forcing coefficient (W/m$^2$)',value=5.35)

def make_slider(label,value,key,
                unit_div=1,make_int=False,
                **kwargs):
    min_value = np.min(param_ranges[key])/unit_div
    max_value = np.max(param_ranges[key])/unit_div

    if make_int:
        min_value = int(min_value)
        max_value = int(max_value)

    if not value:
        value = np.mean([min_value,max_value])
    
    return st.slider(label=label,value=value,
                     min_value=min_value,max_value=max_value,
                     **kwargs)

with tab2:
    A0 = make_slider(label='Area (Mkm$^2$)',value=6.8,key="A0",unit_div=1e12)
    A0 *= 1e12
    B0 = make_slider(label='Total height (km)',value=2.,key="B0",unit_div=1e3)
    B0 *= 1e3
    emp_dur = make_slider('Emplacement duration (Myr)',value=0.5,key="emp_dur")
    erup_freq = st.slider(label='Eruption hiatus duration (Myr)',value=0.25,min_value=dt,max_value=emp_dur)
    degass = make_slider(label='Total volatile release (10$^{18}$ mol C)',value=None,key="degass")

with tab3:
    P0 = make_slider(label='Initial regolith production rate (m/Myr)',value=330,key="P0",make_int=True,)
    E_P = make_slider(label='Initial erosion-to-regolith-production ratio',value=1.,key="E_P",
                      help="<1 → erosion-limited weathering, ≥1 → kinetic-limited weathering")
    
    st.write("Feedback sensitivities (power law exponents):")
    n = make_slider('Global silicate weathering feedback',value=0.5,key="n")
    n_p = make_slider('Local LIP silicate weathering feedback',value=0.4,key="n_p")
    n_e = make_slider('Local LIP erosion-climate feedback',value=0.4,key="n_e")
    c = make_slider('Local LIP erosion-relief feedback',value=0.5,key="c")


with col2:
    t_max = st.number_input(label='Run time (Myr)',
                                value=3,step=1,min_value=1)

t,N,B,H,P,E,degass_arr = run_model(dt=dt,t_max=t_max, # model setup
                           emp_dur=emp_dur,A0=A0, # LIP emplacement characteristics
                           B0=B0,erup_freq=erup_freq, # LIP emplacement characteristics
                           degass=degass, # LIP degassing characteristics
                           P0=P0,E_P=E_P,d=d,c=c,Xm=Xm, # Rock weathering characteristics
                           N0=N0,V=V, # Background climate
                           n=n,n_p=n_p,n_e=n_e, # Climate sensitivities
                           prognostics=True)

T = T0+x(N,N0,b,a)


with col2:
    fig,ax = plt.subplots(figsize=(6.4,4))
    ax.plot(t,T,color='k')
    ax.set_ylabel('Temperature (K)')
    ax.set_xlabel('Time (Myr)')
    ax.axhline(280,c='cadetblue',label='Snowball threshold')
    ax.set_ylim(275,290)
    ax.legend(loc='lower left',frameon=False)
    st.write(fig)

    tab0,tab1,tab2,tab3 = st.tabs(["More plots:","Weathering & erosion","LIP height","Regolith"])

    with tab1:
        fig,ax = plt.subplots(figsize=(6.4,4))
        ax.plot(t,P,color='rosybrown',label='Regolith production')
        ax.plot(t,E,color='grey',label='Erosion')
        ax.set_ylabel('Rate (m/Myr)')
        ax.set_xlabel('Time (Myr)')
        ax.legend()
        st.write(fig)

    with tab2:
        fig,ax = plt.subplots(figsize=(6.4,4))
        ax.fill_between(t,0,B,color='grey',label='Unweathered basalt')
        ax.fill_between(t,B,B+H,color='rosybrown',label='Regolith')
        ax.set_ylabel('Height (m)')
        ax.set_xlabel('Time (Myr)')
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncols=2)
        st.write(fig)

    with tab3:
        fig,ax = plt.subplots(figsize=(6.4,4))
        ax.fill_between(t,0,H,color='rosybrown',label='Regolith')
        ax.set_ylabel('Regolith thickness (m)')
        ax.set_xlabel('Time (Myr)')
        st.write(fig)

def footer(text):
    st.markdown(
    f"<p style='text-align: center; color: grey; font-size: 12px;'>{text}</p>",
    unsafe_allow_html=True)


""""""
footer("Contact: Charlotte Minsky, cminsky@g.harvard.edu")
footer("Source code: <a href='https://github.com/cminsky/snowball-LIP' style='color: grey;'>https://github.com/cminsky/snowball-LIP</a>")