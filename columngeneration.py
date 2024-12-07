from masterproblem import *
from setup import *
from subproblem import *
import pandas as pd
import numpy as np
from compactmodel import *

random.seed(22236)
R_p, Entry_p, Max_t, P, D, T, w, W, W_min, E_app, E_min = generate_dicts(26, 14, 9, 4, 4)
N_c, R_p_c, Entry_p_c = get_unique_combinations_and_list_with_dicts(R_p, Entry_p)

# **** Prerequisites ****
# Create Dataframes
data = pd.DataFrame({
    'P': P + [np.nan] * (max(len(P), len(T), len(D)) - len(P)),
    'T': T + [np.nan] * (max(len(P), len(T), len(D)) - len(T)),
    'D': D + [np.nan] * (max(len(P), len(T), len(D)) - len(D))
})

## Dataframe
results = pd.DataFrame(columns=['P', 'D', 'T', 'objective_value', 'time', 'gap', 'mip_gap'])
results_cg = pd.DataFrame(columns=['it', 'P', 'D', 'T', 'objective_value', 'time', 'lagrange', 'lp-bound'])

# **** Column Generation ****
# Prerequisites
modelImprovable = True
reached_max_itr = False
time_Limit = 1800
max_itr = 200
threshold = 1e-7

# Get Starting Solutions
problem_start = Problem(data, W, W_min, R_p, Entry_p, Max_t, E_app, E_min)
problem_start.buildModel()
problem_start.solveStart()

# Create
start_e = {(p, d): problem_start.e[p, d].x for p in P for d in D}
start_l = {(p, d): problem_start.l[p, d].x for p in P for d in D}
start_x = {(p, t, d): problem_start.x[p, t, d].x for p in P for t in T for d in D}
start_y = {(p, d): problem_start.y[p, d].x for p in P for d in D}
start_z = {(p, t): problem_start.z[p, t].x for p in P for t in T}
start_w = {(p, d): problem_start.w[p, d].x for p in P for d in D}
start_LOS = {(p): problem_start.LOS[p].x for p in P}

while True:
    # Initialize iterations
    itr = 0
    total_subproblem_count = 0
    t0 = time.time()
    last_itr = 0

    # Create empty results lists
    histories = ['objValHistSP', 'timeHist', 'objValHistRMP', 'avg_rc_hist', 'lagrange_hist', 'sum_rc_hist',
                 'avg_sp_time', 'gap_rc_hist']
    histories_dict = {}
    for history in histories:
        histories_dict[history] = []
    objValHistSP, timeHist, objValHistRMP, avg_rc_hist, lagrange_hist, sum_rc_hist, avg_sp_time, gap_rc_hist = histories_dict.values()

    # Create Dicts
    Y_Schedules = plan_dict(start_y, P, None, D)
    Z_Schedules = plan_dict(start_z, P, T, None)
    W_Schedules = plan_dict(start_w, P, None, D)
    LOS_Schedules = plan_dict(start_LOS, P, None, None)
    L_Schedules = plan_dict(start_l, P, None, D)
    E_Schedules = plan_dict(start_e, P, None, D)
    X_Schedules = plan_dict(start_x, P, T, D)

    master = MasterProblem(data, max_itr, itr, Max_t)
    master.buildModel()
    master.initCoeffs()
    master.startSol(start_x, start_LOS)

    # Initialize and solve relaxed model
    master.solRelModel()

    # Retrieve dual values
    duals_td0, duals_p0 = master.getDuals()

    while (modelImprovable) and itr < max_itr:
        itr += 1

        # Solve RMP
        master.current_iteration = itr + 1
        master.solRelModel()

        # Get and Print Duals
        duals_td, duals_p = master.getDuals()

        # Solve SPs
        modelImprovable = False
        for index in P:
            # Build SP
            subproblem = Subproblem(duals_p, duals_td, data, index, itr, R_p, Entry_p, W, W_min, E_app, E_min)
            subproblem.buildModel()

            # Save time to solve SP
            subproblem.solModel()

            reducedCost = subproblem.Model.objval

            # Increase latest used iteration
            last_itr = itr + 1

            # Generate and add columns with reduced cost
            if reducedCost < -threshold:
                Schedules_x = subproblem.getOptVals('x')
                Schedules_LOS = subproblem.getOptVals('LOS')
                master.addCol(index, itr, Schedules_x, Schedules_LOS)
                modelImprovable = True

            master.Model.update()

    if not modelImprovable:
        break

    if modelImprovable and itr == max_itr:
        break

        # Update Model
        master.Model.update()

        if not modelImprovable:
            break

    if modelImprovable and itr == max_itr:
        break

master.finSol()