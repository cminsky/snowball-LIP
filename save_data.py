from dependencies import *
from defaults import *

# save data from an MC run into a NetCDF file
def save_data(filename,t,T,params,description=None,print_progress=False,param_metadata=param_metadata):
    with Dataset(filename, "w", format="NETCDF4") as ncfile:
        
        ## basic setup  and saving results
        if print_progress:
            print("Saving results...")
        
        # add comments
        ncfile.description = description

        # iteration dimension/variable
        iters_dim = ncfile.createDimension("iterations", len(T))
        #iter_var = ncfile.createVariable("iteration", "i4", ("iterations",))
        #iter_var[:] = np.arange(len(T))

        # time dimension/variable
        time_dim = ncfile.createDimension("time_steps", len(t))
        time_var = ncfile.createVariable("time", "f8", ("time_steps",))        
        time_var[:] = t
        time_var.units = "Myr"
        
        # temperature results variable
        temp_var = ncfile.createVariable("temperature", "f8", ("iterations", "time_steps"))
        temp_var[:, :] = T
        temp_var.units = "K"

        ## saving parameters
        if print_progress:
            print("Saving parameters...")
            
        params_dim = ncfile.createDimension("parameters", len(params.keys()))
        params_group = ncfile.createGroup("params")

        # initial record
        param_vars = {}
        for key, val in params.items():
            param_vars[key] = params_group.createVariable(key, "f8", ("iterations",))
            param_vars[key][:] = np.array(val)
            if key in param_metadata:
                param_vars[key].long_name = param_metadata[key]["long_name"]
                param_vars[key].units = param_metadata[key]["units"]
                
        # change units of area and height
        params_group.variables["A0"][:] /= 1e12 # m to Mkm2
        params_group.variables["A0"].units = "Mkm^2"
        params_group.variables["B0"][:] /= 1e3 # m to km
        params_group.variables["B0"].units = "km"
        
        # new parameter: erosion rate (not as ratio)
        E0_var = params_group.createVariable("E0", "f8", ("iterations",))
        E0_var[:] = params_group.variables["P0"][:]*params_group.variables["E_P"][:]
        E0_var.units = params_group.variables["P0"].units
        E0_var.long_name = "Erosion rate"
        
        # new parameter: number of eruptions
        erup_var = params_group.createVariable("erup_num", "f8", ("iterations",))
        erup_var[:] = params_group.variables["emp_dur"][:]/params_group.variables["erup_freq"][:]
        erup_var.units = "dimensionless"
        erup_var.long_name = "Number of eruptions"
        
        print(f"Data saved to {filename}")
        
# analyze for number of Snowballs and other info
def summary_stats(filename,
                  threshold_temp=280,  # temperature threshold [K]
                  min_time=0.9,  # min time to snowball [Myr]
                  max_time=2.15,  # max time to snowball [Myr]
                  print_progress=True):
    
    with Dataset(filename, "a") as ncfile:
        #if print_progress:
            #print("Analyzing, creating summary statistics...")
        stats_group = ncfile.createGroup("stats")

        # read in data
        T = ncfile.variables["temperature"][:]
        t = ncfile.variables["time"][:]
        
        # normalized temperature
        norm_temp_var = stats_group.createVariable("normed_temp", "f8", ("iterations", "time_steps"))
        for i in range(len(T)):
            norm_temp_var[i, :] = T[i, :] - T[i, 0]
        norm_temp_var.units = "K"
        norm_temp_var.long_name = "Normalized temperature (relative to initial value)"

        # boolean flags for conditions
        early_index = (t >= min_time).argmax()
        late_index = (t >= max_time).argmax()

        late_condition = (T[:, :late_index] < threshold_temp).any(axis=1)
        early_condition = (T[:, :early_index] < threshold_temp).any(axis=1)

        snowball_condition = late_condition & ~early_condition

        # min temperature before max time
        min_temp_var = stats_group.createVariable("min_temp", "f8", ("iterations",))
        min_temp_var[:] = T[:, :late_index].min(axis=1)
        min_temp_var.units = "K"
        min_temp_var.long_name = f"Min temp during model runs before {max_time} Myr"

        # min normalized temperature before max time
        min_normed_temp_var = stats_group.createVariable("min_normed_temp", "f8", ("iterations",))
        min_normed_temp_var[:] = norm_temp_var[:, :late_index].min(axis=1)
        min_normed_temp_var.units = "K"
        min_normed_temp_var.long_name = f"Min normalized temp during model runs before {max_time} Myr"
        
        if print_progress:
            print(f"Average cooling: {np.nanmean(min_normed_temp_var[:]):0.2f} K")

        # save flags
        late_flag_var = stats_group.createVariable("late_flag", "i1", ("iterations",))
        late_flag_var[:] = late_condition.astype("i1")
        late_flag_var.long_name = f"Flag if temp dropped below {threshold_temp}K before {max_time} Myr"
        late_flag_var.units = "boolean"

        early_flag_var = stats_group.createVariable("early_flag", "i1", ("iterations",))
        early_flag_var[:] = early_condition.astype("i1")
        early_flag_var.long_name = f"Flag if temp dropped below {threshold_temp}K before {min_time} Myr"
        early_flag_var.units = "boolean"

        snowball_flag_var = stats_group.createVariable("snowball_flag", "i1", ("iterations",))
        snowball_flag_var[:] = snowball_condition.astype("i1")
        snowball_flag_var.long_name = f"Flag if temp dropped below {threshold_temp}K between {min_time} and {max_time} Myr"
        snowball_flag_var.units = "boolean"

        # Save summary percentages
        iters = len(T)
        percentages = {
            "late_flag_percentage": (late_condition.sum() / iters * 100),
            "early_flag_percentage": (early_condition.sum() / iters * 100),
            "snowball_flag_percentage": (snowball_condition.sum() / iters * 100),
        }

        for name, value in percentages.items():
            var = stats_group.createVariable(name, "f8")
            var[:] = value
            var.long_name = f"Percentage of iterations for {name.replace('_', ' ')}"
            var.units = "%"

        if print_progress:
            for name, value in percentages.items():
                print(f"{name.replace('_', ' ').capitalize()}: {value:.2f}%")

            #print(f"Summary statistics added to {filename}")
            
# read in summary stats
def read_summary_stats(filename):
    with Dataset(filename, "r") as ncfile:
        print(f"Reading summary statistics from {filename}")
        
        stats_group = ncfile.groups["stats"]
        
        min_normed_temp = stats_group.variables["min_normed_temp"][...]
        late_flag_percentage = stats_group.variables["late_flag_percentage"][...]
        early_flag_percentage = stats_group.variables["early_flag_percentage"][...]
        snowball_flag_percentage = stats_group.variables["snowball_flag_percentage"][...]
        
        print(f"Average cooling: {np.nanmean(min_normed_temp):0.2f} K")
        print(f"Late flag percentage: {late_flag_percentage:.2f}%")
        print(f"Early flag percentage: {early_flag_percentage:.2f}%")
        print(f"Snowball flag percentage: {snowball_flag_percentage:.2f}%")