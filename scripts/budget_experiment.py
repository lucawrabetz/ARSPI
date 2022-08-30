import subprocess
import pandas as pd
import seaborn as sns
from datetime import date
import os

COLUMNS = ["setname", "n", "p", "k", "r_0", "objective_value", "runtime"]
EXPERIMENTS = "experiments"
# ./bin/sp tests-06_20_22-0 2 3 32 3 1 5 500 50 200 0.1 0 3

def check_make_dir(path, i):
    """
    Recursively check if an experiment directory exists, or create one with the highest number
        - example - if "path" string is "/dat/experiments/test-01_29_22", and there already exist:
            - "/dat/experiments/test-01_29_22-0"
            - "/dat/experiments/test-01_29_22-1"
            - "/dat/experiments/test-01_29_22-2"
        we have to create the dir "/dat/experiments/test-01_29_22-3"
    """
    isdir = os.path.isdir(path + "-" + str(i))
    # if the directory exists, call on the next i
    if isdir:
        return check_make_dir(path, i + 1)
    # base case - create directory for given i (and return final path)
    else:
        os.mkdir(path + "-" + str(i))
        return path + "-" + str(i)

def single_run(setname, p_0, k_0, n, p, k, r_0, M=500, cost_a=20, cost_b=200, fraction=0.1, costs=1, dist=3):
    sp_call = ["./bin/sp", setname, str(p_0), str(k_0), str(n), str(p), str(k), str(r_0), str(M), str(cost_a), str(cost_b), str(fraction), str(costs), str(dist)]
    sp_result = subprocess.run(sp_call, stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')[-2]
    objective_value = sp_result.split()[0]
    runtime = int(sp_result.split()[1]) * 0.001 # (in seconds)
    result = [setname, str(n), str(p), str(k), str(r_0), str(objective_value), str(runtime)]
    return result

def budget_experiment(setname, p_0, k_0, n, p, r_0_list):
    today = date.today()
    date_str = today.strftime("%m_%d_%y")
    exp_path = os.path.join(EXPERIMENTS, "budget-same", setname+"-"+date_str+"-density"+str(p_0)+"-n"+str(n))
    exp_path_finished = check_make_dir(exp_path, 0)
    final_path = os.path.join(exp_path_finished, "results.csv")
    fig_path = os.path.join(exp_path_finished, "figure.png")
    results = []
    first = True
    for r_0 in r_0_list:
        for k in [1, 2, 3]:
            if (first):
                single_result = single_run(setname, p_0, k_0, n, p, k, r_0, costs=1)
                first = False
            else:
                single_result = single_run(setname, p_0, k_0, n, p, k, r_0, costs=0)
            results.append(single_result)
    results_df = pd.DataFrame(results, columns=COLUMNS)
    results_df.to_csv(final_path)

    # plotting
    sns.set_context("paper")
    sns.set(style="darkgrid")
    print(results_df)
    return results_df

def main():
    setname = "tests-06_20_22-0"
    p_0_list = [1, 2, 3, 4]
    k_0 = 3
    n_list = [122]
    p = 3
    k = 1
    for n in n_list:
        r_0_list = range(1, int(n / 2))
        for p_0 in p_0_list:
            for i in range(3):
                results_df = budget_experiment(setname, p_0, k_0, n, p, r_0_list)

if __name__=="__main__":
    main()
