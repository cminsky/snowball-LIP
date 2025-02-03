from dependencies import *

### default parameter values/ranges and metadata

dt = 1e3 / 1e6  # Myr
t_max = 3.0     # Myr

N_pi = 2.83 # pre-industrial surficial C (examol)
ppm_pi = 280 # pre-industrial pCO2 (ppm)

# Franklin LIP specific
t_Earth = 719 # Myr ago

param_ranges = {
    # Franklin LIP specific
    "A0": (2.62e12,11e12),        # area [m2]

    # LIP emplacement characteristics
    "emp_dur": (dt, 1),           # length of emplacement [Myr]
    "B0": (1e2, 4e3),            # LIP height [m]
    "degass": (0, 10),           # [examol CO2]

    # LIP weathering characteristics
    "P0": (0, 400),              # optimal saprolite production [m/Myr]
    "E_P": (0, 10),              # E0/P0 ratio []
    "d": 0.5,                    # fixed weathering attenuation height [m]
    "c": (0, 1),                 # erosion dependence on slope []
    "Xm": 10317,                 # fixed Ca Mg fraction [mol/m3]

    # Climate sensitivities
    "n": (0, 1),                 # global weathering feedback strength []
    "n_p": (0, 1),               # LIP weathering feedback strength []
    "n_e": (0, 1),               # LIP erosion feedback strength []

    # Temperature response
    "b": (5.22, 5.76),           # CO2 radiative forcing [W/m2]
    "a": 2                       # fixed temperature response radiative forcing [W/m2/K]
}

# restate constants for easier integration with other parts of code
d = param_ranges["d"]
Xm = param_ranges["Xm"]
a = param_ranges["a"]


# metadata
param_metadata = {
    # Background climate
    "N0": {"long_name": "Initial surficial carbon $N_0$", "units": "examol"},
    "T0": {"long_name": "Initial temperature $T_0$", "units": "K"},
    "b": {"long_name": "CO2 radiative forcing coefficient $b$", "units": "W/m^2"},
    
    # Degassing
    "V": {"long_name": "Global degassing flux $V$", "units": "examol/Myr"},
    "degass": {"long_name": "$CO_2$ release", "units": "examol"},
    
    # Emplacement
    "A0": {"long_name": "LIP area $A_0$", "units": "m^2"},
    "B0": {"long_name": "LIP height $B_0$", "units": "m"},
    "emp_dur": {"long_name": "Emplacement duration", "units": "Myr"},
    "erup_freq": {"long_name": "Hiatus duration", "units": "Myr"},
    
    # Weathering and erosion
    "P0": {"long_name": "Regolith production rate $P_0$", "units": "m/Myr"},
    "E_P": {"long_name": "Ratio of erosion to regolith production rate $E_0/P_0$", "units": ""},
    "c": {"long_name": "Dependence of erosion on relief $c$", "units": ""},
    
    # Climate feedback sensitivities
    "n": {"long_name": "Global weathering feedback $n$", "units": ""},
    "n_p": {"long_name": "LIP weathering feedback $n_p$", "units": ""},
    "n_e": {"long_name": "LIP erosion feedback strength $n_e$", "units": ""},
}