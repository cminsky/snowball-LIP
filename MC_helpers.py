from dependencies import *
from model import *
from background import *

### helper functions to perform Monte Carlo sampling

# generate parameter values to iterate through
def sample_params(iters,param_ranges,dt,t_max):
    params = {}
    
    # fixed timing parameters
    params['dt'] = [dt]*iters
    params['t_max'] = [t_max]*iters
    
    # randomly sample from ranges
    for key,val in param_ranges.items():
        if isinstance(val,tuple):
            params[key] = np.random.uniform(*val,iters)
        else:
            params[key] = [val]*iters
            
    # eruption frequency depends on emplacement duration
    params["erup_freq"] = np.random.uniform(dt,params["emp_dur"])
    
    return params

# split the parameter dictionary into parameters for the model vs. x function
def slice_params(params):
    model_params = {key: params[key] for key in params if key not in ['b','a','T0']} # arguments for model function
    x_params = {key: params[key] for key in ['N0','b','a']} # arguments for x function
    return model_params, x_params

# turn a dictionary of arrays into an array of dictionaries
def dict_array(d):
    length = len(list(d.values())[0])
    return[{key:d[key][i] for key in d} for i in range(length)]

# print time in seconds or minutes
def sec_min_str(time_sec):
    return f"{int(time_sec // 60)} minutes" if time_sec >= 60 else f"{int(time_sec)} seconds"

# run the model (not using this function - replaced by parallelized version)
def run_model_iters(params):
    model_params, x_params = slice_params(params)
    model_array, x_array = dict_array(model_params),dict_array(x_params)
    results = []
    
    iters = len(model_array)
    for i in range(iters):
        t,N = run_model(**model_array[i])
        T = params['T0'][i] + x(N,**x_array[i])
        results.append(T)
    
    return np.array(t),np.array(results)

### MAIN MC FUNCTIONS ###
def run_single_iteration(i,model_array,x_array,T0):
    t, N = run_model(**model_array[i])
    T = T0[i] + x(N, **x_array[i])
    return t, T

def run_model_iters_parallel(params, print_progress=True):
    # slice the dictionary into model function vs. x function and rearrange
    model_params, x_params = slice_params(params)
    model_array, x_array = dict_array(model_params), dict_array(x_params)

    # setup
    iters = len(model_array)
    T0 = params['T0']

    if print_progress:
        # set up progress tracking
        if iters < 500:
            print_interval = 10
        elif iters < 100000:
            print_interval = 30
        else:
            print_interval = 120

        start_time = time.time()
        next_print_time = start_time + print_interval

        est_sec = iters * 10 / 250  # observed about 250 iters every 10 seconds
        start_time_str = datetime.fromtimestamp(start_time).strftime('%I:%M %p')
        est_time_str = datetime.fromtimestamp(start_time + est_sec).strftime('%I:%M %p')
        print(f"Starting {iters:,} iterations at {start_time_str}")
        print(f"Estimating completion in {sec_min_str(est_sec)} at {est_time_str}")

    # parallelization
    results = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {
            executor.submit(run_single_iteration, i, model_array, x_array, T0): i
            for i in range(iters)
        }
        for future in concurrent.futures.as_completed(futures):
            iteration_index = futures[future]
            try:
                result = future.result()  # get (t, T)
                results.append((iteration_index, *result))  # include iteration index to sort later

                # print progress every time interval
                if print_progress:
                    current_time = time.time()
                    if current_time >= next_print_time:
                        elapsed = current_time - start_time
                        completed = len(results)
                        print(f"    Completed {completed:,}/{iters:,} iterations after {sec_min_str(elapsed)}")
                        next_print_time += print_interval
                        
            except Exception as e:
                print(f"Iteration {iteration_index} generated an exception: {e}")

    # ending print statements
    if print_progress:
        current_time = time.time()
        elapsed = current_time - start_time
        end_time_str = datetime.fromtimestamp(current_time).strftime('%I:%M %p')
        print(f"Completed {iters:,} iterations after {sec_min_str(elapsed)} at {end_time_str}")

    # combine results
    results.sort(key=lambda x: x[0])  # sort by iteration index to preserve order alignment with parameters
    _, t_all, T_all = zip(*results)

    t = np.array(t_all[0])
    T = np.array(T_all)

    return t, T