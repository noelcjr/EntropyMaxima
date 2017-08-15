# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 20:25:03 2016

@author: noel
"""
import pandas as pd
import numpy as np
from collections import Counter
from random import randint
file_path = '/home/noel/Projects/Job_Interviews/2ndGenome/scop2_processed_data.csv'
df_pdb = pd.read_csv(file_path)

df_pdb = df_pdb.set_index(['Unnamed: 0'])
df_pdb.index.names = ['indx']
len(df_pdb.Structural_Class[df_pdb.Structural_Class == 'alpha'])          # 156
len(df_pdb.Structural_Class[df_pdb.Structural_Class == 'beta'])           # 52
len(df_pdb.Structural_Class[df_pdb.Structural_Class == 'alpha_and_beta']) # 224
len(df_pdb.Structural_Class[df_pdb.Structural_Class == 'alpha_plus_beta'])# 236
len(df_pdb.Structural_Class[df_pdb.Structural_Class == 'small'])          # 7
# I will get rid of small proteins because they do not have secondary structure info
t = []
drop = []
for i in range(df_pdb.shape[0]):
    if df_pdb.ix[i]['Structural_Class'] == 'alpha':
        t.append(1)
    elif df_pdb.ix[i]['Structural_Class'] == 'beta':
        t.append(2)
    elif df_pdb.ix[i]['Structural_Class'] == 'alpha_and_beta':
        t.append(3)
    elif df_pdb.ix[i]['Structural_Class'] == 'alpha_plus_beta':
        t.append(4)
    elif df_pdb.ix[i]['Structural_Class'] == 'small':
        drop.append(i)
    else:
        print(df_pdb.ix[i]['Structural_Class'], 'Unknown Structural Class')
#  Drop small proteins and renumber index
df_pdb = df_pdb.drop(list(set(drop)))
df_pdb = df_pdb.reset_index(drop=True)
df_pdb['x'] = df_pdb.index
prot_struct_aa = df_pdb[['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE', \
                         'LYS','LEU','MET','ASN','PRO','GLN','ARG','SER', \
                         'THR','VAL','TRP','TYR','x','Structural_Class']]
prot_struct_tp1 = df_pdb[['Hydro','Polar','Charged','Structural_Class']]
###############################################################################
####################   KNN Euclidian Distance #################################
def euclidean_distance(np1, np2):
    return np.linalg.norm(np1-np2)
#euclidean_distance(np.array([1,2,3,4]),np.array([1,3,4,8]))
def predict(current_df,unknown, k = 3):
    '''
    Input:
        unknown  == four attributes of an unknown flower
        k        == the number of neighbors used
    Output:
        A prediction of the species of flower (str)
    '''
    distances = [(euclidean_distance(unknown, row[:-2]),row[-1]) for row in current_df]
    nearest = sorted(distances)[:k]

    return Counter([n[1] for n in nearest]).most_common(1)[0][0]
###############################################################################
'''Putting it all together KNN only'''
n = 20
total_cv_results = []
All_results = pd.DataFrame(columns=['KNN_cv','KNN_t','KNN_HMM_C_vc','KNN_HMM_C_t','KNN_HMM_O_cv','KNN_HMM_O_t'])
for k in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]:
    results = []
    np.random.seed(10)
    for i in range(10):
        l = np.random.choice(range(550), 550, replace=False)
        current_df = prot_struct_aa.values[l]
        predictions = np.array([predict(current_df[0:412],row[:n],k) for row in current_df[412:]])
        actual = np.array([row[-1] for row in current_df[412:]])
        results.append(np.mean(predictions == actual))
    total_cv_results.append((k,sum(results)/len(results)))
    All_results.loc[k,'KNN_cv'] = sum(results)/len(results)
    print(k,sum(results)/len(results))
  # best result (10, 0.61884057971014494)
total_cv_results = []
np.random.seed(10)
k10 = []
for k in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]:
    results = []
    np.random.seed(10)
    for i in range(10):
        l = np.random.choice(range(len(prot_struct_aa)), len(prot_struct_aa), replace=False)
        current_df = prot_struct_aa.values[l]
        predictions = np.array([predict(current_df[0:550],row[:n],k) for row in current_df[550:]])
        actual = np.array([row[-1] for row in current_df[550:]])
        if k == 10:
            k10.append(np.mean(predictions == actual))
        results.append(np.mean(predictions == actual))
    total_cv_results.append((k,sum(results)/len(results)))
    All_results.loc[k,'KNN_t'] = sum(results)/len(results)
    #np.mean(predictions == actual)   # k 10 0.65254237288135597
    # k10 =
    #[0.65254237288135597,
    # 0.61864406779661019,
    # 0.5847457627118644,
    # 0.55932203389830504,
    # 0.61016949152542377,
    # 0.65254237288135597,
    # 0.53389830508474578,
    # 0.66101694915254239,
    # 0.69491525423728817,
    # 0.56779661016949157]
###############################################################################
def viterbi_display(obs, states, start_p, trans_p, emit_p):
    V = [{}]
    for st in states:
        V[0][st] = {"prob": start_p[st] * emit_p[st][obs[0]], "prev": None}
    # Run Viterbi when t > 0
    for t in range(1, len(obs)):
        V.append({})
        for st in states:
            max_tr_prob = max(V[t-1][prev_st]["prob"]*trans_p[prev_st][st] for prev_st in states)
            for prev_st in states:
                if V[t-1][prev_st]["prob"] * trans_p[prev_st][st] == max_tr_prob:
                    max_prob = max_tr_prob * emit_p[st][obs[t]]
                    V[t][st] = {"prob": max_prob, "prev": prev_st}
                    break
    for line in dptable(V):
        print line
    opt = []
    # The highest probability
    max_prob = max(value["prob"] for value in V[-1].values())
    previous = None
    # Get most probable state and its backtrack
    for st, data in V[-1].items():
        if data["prob"] == max_prob:
            opt.append(st)
            previous = st
            break
    # Follow the backtrack till the first observation
    for t in range(len(V) - 2, -1, -1):
        opt.insert(0, V[t + 1][previous]["prev"])
        previous = V[t + 1][previous]["prev"]

    print 'The steps of states are ' + ' '.join(opt) + ' with highest probability of %s' % max_prob
def viterbi(obs, states, start_p, trans_p, emit_p):
    V = [{}]
    for st in states:
        V[0][st] = {"prob": start_p[st] * emit_p[st][obs[0]], "prev": None}
    # Run Viterbi when t > 0
    for t in range(1, len(obs)):
        V.append({})
        for st in states:
            max_tr_prob = max(V[t-1][prev_st]["prob"]*trans_p[prev_st][st] for prev_st in states)
            for prev_st in states:
                if V[t-1][prev_st]["prob"] * trans_p[prev_st][st] == max_tr_prob:
                    max_prob = max_tr_prob * emit_p[st][obs[t]]
                    V[t][st] = {"prob": max_prob, "prev": prev_st}
                    break
    opt = []
    # The highest probability
    max_prob = max(value["prob"] for value in V[-1].values())
    previous = None
    # Get most probable state and its backtrack
    for st, data in V[-1].items():
        if data["prob"] == max_prob:
            opt.append(st)
            previous = st
            break
    # Follow the backtrack till the first observation
    for t in range(len(V) - 2, -1, -1):
        opt.insert(0, V[t + 1][previous]["prev"])
        previous = V[t + 1][previous]["prev"]
    return opt

def viterbi2(obs, states, start_p, trans_p, emit_p):
    V = [{}]
    for st in states:
        V[0][st] = {"prob": start_p[st] * emit_p[st][obs[0]]['1'], "prev": None}
    # Run Viterbi when t > 0
    for t in range(1, len(obs)):
        V.append({})
        for st in states:
            max_tr_prob = max(V[t-1][prev_st]["prob"]*trans_p[prev_st][st] for prev_st in states)
            for prev_st in states:
                if V[t-1][prev_st]["prob"] * trans_p[prev_st][st] == max_tr_prob:
                    max_prob = max_tr_prob * emit_p[st][obs[t]][str(t+1)]
                    V[t][st] = {"prob": max_prob, "prev": prev_st}
                    break
    opt = []
    # The highest probability
    max_prob = max(value["prob"] for value in V[-1].values())
    previous = None
    # Get most probable state and its backtrack
    for st, data in V[-1].items():
        if data["prob"] == max_prob:
            opt.append(st)
            previous = st
            break
    # Follow the backtrack till the first observation
    for t in range(len(V) - 2, -1, -1):
        opt.insert(0, V[t + 1][previous]["prev"])
        previous = V[t + 1][previous]["prev"]
    return opt

def dptable(V):
    # Print a table of steps from dictionary
    yield " ".join(("%10d" % i) for i in range(len(V)))
    for state in V[0]:
        yield "%.7s: " % state + " ".join("%.7s" % ("%f" % v[state]["prob"]) for v in V)
#############################################################################################
''' KNN with HMM'''
def train_HMM(current_df,unknown, k = 3):
    '''This checks neighbors sequences and transitions winth a sample for traing (%70)
    The main difference is that it avoids calculating the distance of some observations to itsel, which would be 0'''
    distances = [(euclidean_distance(row[:-2],unknown[:-1]),row[-1],row[-2],unknown[-1]) for row in current_df]
    # sorting by first field in tuple
    nearest = sorted(distances, key=lambda tup: tup[0])[:k+1]
    return nearest
def observe_HMM(current_df,unknown, k = 3):
    '''Similar to predict'''
    distances = [(euclidean_distance(row[:-2],unknown),row[-1]) for row in current_df]
    nearest = sorted(distances, key=lambda tup: tup[0])[:k]
    return nearest
########################################################################################################################
# The following is a markov model using a count of states observed in the readings plus a count of transition of states
# because of the nature of the model design. It turns out that the emision_probabilities and transition_probabilities
# are identical. This is add, but it gace results in the 50% ranges. So I kept the results but I tried anothe HMM model
# after this.
n = 20
total_cv_results = []
k = 10
t_set = 550
full_set = len(prot_struct_aa)
best_from_CV = 4
emmision_prop_calc = ['count','distance']
epc = emmision_prop_calc[0]
for k in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]:
    results = []
    np.random.seed(10)
    for ite in range(10):
        states = ('alpha','beta','alpha_and_beta','alpha_plus_beta')
        # We will do it for k - 10 to test the code, then we will do cross validation
        observations = ('1','2','3','4','5','6','7','8','9','10')
        # we know how many of each stae are in the test set, this knowledge could bias the HMM,
        # but we can say that there are more alpha beta combinations than alpha and betas alone
        start_probability = {'alpha': 0.2, 'beta': 0.1, 'alpha_and_beta':0.3, 'alpha_plus_beta': 0.4}
        # Transition and emission probabilities need to be calculated from 412 structures first, then 550
        # They are initialized to 0
        transition_probability = {'alpha' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0},\
                                  'beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0},\
                                  'alpha_and_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0},\
                                  'alpha_plus_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0}}
        emission_probability = {'alpha' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'alpha_and_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'alpha_plus_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0}}
        emission_counter = {'alpha' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'alpha_and_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'alpha_plus_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0}}
        emission_max = {'alpha' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'alpha_and_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'alpha_plus_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0}}    
        l = np.random.choice(range(full_set), full_set, replace=False)
        current_df = prot_struct_aa.values[l]
        traing_set_for_HMMmodel = np.array([train_HMM(current_df[0:t_set],row[:n+1],k) for row in current_df[t_set:]])
        for i  in traing_set_for_HMMmodel:
            state = i[0][1]
            for j in i[1:]:
                transition_probability[state][j[1]] += 1
                if epc == 'count':
                    emission_probability[state][j[1]] += 1
                elif epc == 'distance':
                    emission_probability[state][j[1]] += float(j[0])
                    emission_counter[state][j[1]] += 1
                    if float(j[0]) > emission_max[state][j[1]]:
                        emission_max[state][j[1]] = float(j[0])
                state = j[1]        
        # Normalization emission_probabilities must add up to observations for training 412x10 = 4120
        for i in states:
            if epc == 'count':
                suma = 0.0
                for j in states:
                    suma += emission_probability[i][j]
                for j in states:
                    emission_probability[i][j] = emission_probability[i][j]/suma
            elif epc == 'distance':
                for j in states:
                    emission_probability[i][j] = emission_probability[i][j]/emission_counter[i][j]
                m = max(emission_probability[i].values())
                for j in states:
                    emission_probability[i][j] = 1 - emission_probability[i][j]/m
                s = sum(emission_probability[i].values())
                for j in states:
                    emission_probability[i][j] = emission_probability[i][j]/s                           
        for i in states:
            suma = 0.0
            for j in states:
                suma += transition_probability[i][j]
            for j in states:
                transition_probability[i][j] = transition_probability[i][j]/suma
        count = 0
        matched = 0
        for row in current_df[t_set:]:
            tup = tuple([tt[1] for tt in observe_HMM(current_df[0:t_set],row[:n],k)])
            results2 = viterbi(tup, states, start_probability, transition_probability, emission_probability)
            count += 1
            if results2[-1] == row[-1]:
                matched += 1
            succes_rate = float(matched)/float(count)
        results.append(succes_rate)
        if ite == 0 and k == 1:
            t_p = transition_probability
            e_p = emission_probability
            s_r = succes_rate
            b_k = k
        elif ite == 0 and k == best_from_CV:
            t_p2 = transition_probability
            e_p2 = emission_probability
            s_r2 = succes_rate
            b_k2 = k            
        else:
            if succes_rate > s_r:
                t_p = transition_probability
                e_p = emission_probability
                s_r = succes_rate
                b_k = k
            if k == best_from_CV and succes_rate > s_r2:
                t_p2 = transition_probability
                e_p2 = emission_probability
                s_r2 = succes_rate
                b_k2 = k
    print(k,results)
    All_results.loc[k,'KNN_HMM_C_t'] = sum(results)/10
########################################################################################################################
# After Cross validation. Test models on data not used for training tests.
n = 20
total_cv_results = []
k = 10
t_set = 550
full_set = len(prot_struct_aa)
best_from_CV = 4
emmision_prop_calc = ['count','distance']
epc = emmision_prop_calc[0]
for k in [3,4]:
    results = []
    np.random.seed(10)
    states = ('alpha','beta','alpha_and_beta','alpha_plus_beta')
    observations = ('1','2','3','4','5','6','7','8','9','10')
    start_probability = {'alpha': 0.2, 'beta': 0.1, 'alpha_and_beta':0.3, 'alpha_plus_beta': 0.4}
    if k == 3:
        transition_probability = t_p
        emission_probability = e_p
    elif k == 4:
        transition_probability = t_p2
        emission_probability = e_p2
    current_df = prot_struct_aa.values[0:full_set]
    count = 0
    matched = 0
    for row in current_df[t_set:]:
        tup = tuple([tt[1] for tt in observe_HMM(current_df[0:t_set],row[:n],k)])
        results2 = viterbi(tup, states, start_probability, transition_probability, emission_probability)
        count += 1
        if results2[-1] == row[-1]:
            matched += 1
        succes_rate = float(matched)/float(count)
    results.append(succes_rate)
    print(k,results)
# Test set gave best HMM model to K=3 in 10 different cross validations, and this model is used un data not used in training.
# Test set gave the K = 4 value to be four from 10 different crossvalidations, from which the Best HMM model is taken 
# and used on untrain data and the result is below
#(k=3, [0.5338983050847458])
#(k=4, [0.4491525423728814])
########################################################################################################################
# k = 1 0.5677966 NN_HMM_C_t            k = 6   KNN_HMM_C_t
# [0.576271186440678,                  [0.5254237288135594,
#  0.5932203389830508,                  0.5338983050847458,
#  0.5338983050847458,                  0.5508474576271186, 
#  0.5508474576271186,                  0.5,
#  0.5254237288135594,                  0.5847457627118644,
#  0.5508474576271186,                  0.559322033898305,
#  0.5508474576271186,                  0.5169491525423728, 
#  0.6016949152542372,                  0.5932203389830508,
#  0.652542372881356,                   0.635593220338983,
#  0.5423728813559322]                  0.5169491525423728
########################################################################################################################
# This HMM model is more sophisticated than the previous one. The trasition_probability matrix is identical as the previous
# method, but the emision_probabilites are calculated differently. Statistics on the positions of states observed are
# kept for every state. When alpha is observe, the probabilities that alpha is 1st, 2nd, 3rd, 4th.. etc is calculated
# relative to that observation in other states. Probabilities add up to one for each position 1st, 2nd etc.. that is
# that posibilities depend on the probability in which observed states are found. Only order2 worked. The rest, gave 
# results close to random guessing.
n = 20
k = 10
t_set = 550
full_set = len(prot_struct_aa)
emmision_prop_calc = ['count','distance','order','order2','consensus']
epc = emmision_prop_calc[3]
for k in range(1,21):
    results = []
    np.random.seed(10)
    for ite in range(10):
        states = ('alpha','beta','alpha_and_beta','alpha_plus_beta')
        # We will do it for k - 10 to test the code, then we will do cross validation
        #observations = ('1','2','3','4','5','6','7','8','9','10')
        observations = tuple([str(i) for i in range(1,k+1)])
        # we know how many of each stae are in the test set, this knowledge could bias the HMM,
        # but we can say that there are more alpha beta combinations than alpha and betas alone
        start_probability = {'alpha': 0.2, 'beta': 0.1, 'alpha_and_beta':0.3, 'alpha_plus_beta': 0.4}
        # Transition and emission probabilities need to be calculated from 412 structures first, then 550
        # They are initialized to 0
        transition_probability = {'alpha' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0},\
                                  'beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0},\
                                  'alpha_and_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0},\
                                  'alpha_plus_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0}}
        emission_probability = {'alpha':{'alpha':{},'beta':{},'alpha_and_beta':{},'alpha_plus_beta':{}},\
                                'beta':{'alpha':{},'beta':{},'alpha_and_beta':{},'alpha_plus_beta':{}},\
                                'alpha_and_beta':{'alpha':{},'beta':{},'alpha_and_beta':{},'alpha_plus_beta':{}},\
                                'alpha_plus_beta':{'alpha':{},'beta':{},'alpha_and_beta':{},'alpha_plus_beta':{}}}
        emission_counter = {'alpha':{},'beta':{},'alpha_and_beta':{},'alpha_plus_beta':{}}
        emission_max = {'alpha':{},'beta':{},'alpha_and_beta':{},'alpha_plus_beta':{}}
        for i in states:
            for j in states:
                for ll in tuple([str(jj) for jj in range(1,k+1)]):
                    emission_probability[i][j][ll] = 0.0
            for j in tuple([str(jj) for jj in range(1,k+1)]):
                emission_counter[i][j] = 0.0
                emission_max[i][j] = 0.0
        l = np.random.choice(range(full_set), full_set, replace=False)
        current_df = prot_struct_aa.values[l]
        traing_set_for_HMMmodel = np.array([train_HMM(current_df[0:t_set],row[:n+1],k) for row in current_df[t_set:]])
        cnt = 0
        for i  in traing_set_for_HMMmodel:
            #if cnt == 549:
            state = i[0][1]
            count_e = 1
            for j in i[1:]:
                transition_probability[state][j[1]] += 1
                if epc == 'count':
                    emission_probability[state][j[1]] += 1
                elif epc == 'distance':
                    emission_probability[state][j[1]] += float(j[0])
                    emission_counter[state][j[1]] += 1
                    if float(j[0]) > emission_max[state][j[1]]:
                        emission_max[state][j[1]] = float(j[0])
                elif epc == 'order':
                    #print(cnt,'order',count_e,i[0][1],j[1],emission_probability[i[0][1]][j[1]][str(count_e)])
                    emission_counter[i[0][1]][str(count_e)] += 1
                    emission_probability[i[0][1]][j[1]][str(count_e)] += 1
                    #print('    ',cnt,'order',count_e,i[0][1],j[1],emission_probability[i[0][1]][j[1]][str(count_e)])
                    count_e += 1
                elif epc == 'order2':
                    #print(cnt,'order',count_e,i[0][1],j[1],emission_probability[i[0][1]][j[1]][str(count_e)])
                    emission_counter[i[0][1]][str(count_e)] += 1
                    emission_probability[i[0][1]][j[1]][str(count_e)] += 1
                    #print('    ',cnt,'order',count_e,i[0][1],j[1],emission_probability[i[0][1]][j[1]][str(count_e)])
                    count_e += 1
                state = j[1]
            #cnt += 1
        # Normalization emission_probabilities must add up to observations for training 412x10 = 4120
        for i in states:
            if epc == 'count':
                suma = 0.0
                for j in states:
                    suma += emission_probability[i][j]
                for j in states:
                    emission_probability[i][j] = emission_probability[i][j]/suma
            elif epc == 'distance':
                for j in states:
                    emission_probability[i][j] = emission_probability[i][j]/emission_counter[i][j]
                m = max(emission_probability[i].values())
                for j in states:
                    emission_probability[i][j] = 1 - emission_probability[i][j]/m
                s = sum(emission_probability[i].values())
                for j in states:
                    emission_probability[i][j] = emission_probability[i][j]/s
            elif epc == 'order':
                for ii in states:
                    suma = 0.0
                    for j in tuple([str(jj) for jj in range(1,k+1)]):
                        suma += emission_probability[i][ii][j]
                    for j in tuple([str(jj) for jj in range(1,k+1)]):
                        emission_probability[i][ii][j] = emission_probability[i][ii][j]/suma
            elif epc == 'order2':                
                for j in tuple([str(jj) for jj in range(1,k+1)]):
                    suma = 0.0
                    for ii in states:
                        suma += emission_probability[i][ii][j]
                    for ii in states:
                        emission_probability[i][ii][j] = emission_probability[i][ii][j]/suma                         
        for i in states:
            suma = 0.0
            for j in states:
                suma += transition_probability[i][j]
            for j in states:
                transition_probability[i][j] = transition_probability[i][j]/suma
        
        count = 0
        matched = 0
        for row in current_df[t_set:]:
            tup = tuple([tt[1] for tt in observe_HMM(current_df[0:t_set],row[:n],k)])
            #results = viterbi(tup, states, start_probability, transition_probability, emission_probability)
            #def viterbi_display(obs, states, start_p, trans_p, emit_p):
            V = [{}]
            for st in states:
                V[0][st] = {"prob": start_probability[st] * emission_probability[st][tup[0]]['1'], "prev": None}
            # Run Viterbi when t > 0
            for t in range(1, len(tup)):
                V.append({})
                for st in states:
                    max_tr_prob = max(V[t-1][prev_st]["prob"]*transition_probability[prev_st][st] for prev_st in states)
                    for prev_st in states:
                        if V[t-1][prev_st]["prob"] * transition_probability[prev_st][st] == max_tr_prob:
                            max_prob = max_tr_prob * emission_probability[st][tup[t]][str(t+1)]
                            V[t][st] = {"prob": max_prob, "prev": prev_st}
                            break
            opt = []
            # The highest probability
            max_prob = max(value["prob"] for value in V[-1].values())
            previous = None
            # Get most probable state and its backtrack
            for st, data in V[-1].items():
                if data["prob"] == max_prob:
                    opt.append(st)
                    previous = st
                    break
            # Follow the backtrack till the first observation
            for t in range(len(V) - 2, -1, -1):
                opt.insert(0, V[t + 1][previous]["prev"])
                previous = V[t + 1][previous]["prev"]
        
            #print 'The steps of states are ' + ' '.join(opt) + ' with highest probability of %s' % max_prob
            count += 1
            if opt[-1] == row[-1]:
                matched += 1
            succes_rate = float(matched)/float(count)
        results.append(succes_rate)
        if ite == 0 and k == 1:
            t_p = transition_probability
            e_p = emission_probability
            s_r = succes_rate
            b_k = k
        elif ite == 0 and k == best_from_CV:
            t_p2 = transition_probability
            e_p2 = emission_probability
            s_r2 = succes_rate
            b_k2 = k            
        else:
            if succes_rate > s_r:
                t_p = transition_probability
                e_p = emission_probability
                s_r = succes_rate
                b_k = k
            if k == best_from_CV and succes_rate > s_r2:
                t_p2 = transition_probability
                e_p2 = emission_probability
                s_r2 = succes_rate
                b_k2 = k
    print(k,results)
    All_results.loc[k,'KNN_HMM_O_t'] = sum(results)/10
########################################################################################################################
# After Cross validation. Test models on data not used for training tests.Both the best HMM model and best average model
# for a given k are 4
n = 20
total_cv_results = []
k = 10
t_set = 550
full_set = len(prot_struct_aa)
best_from_CV = 4
emmision_prop_calc = ['count','distance']
epc = emmision_prop_calc[0]
for k in [6,4]:
    results = []
    np.random.seed(10)
    states = ('alpha','beta','alpha_and_beta','alpha_plus_beta')
    observations = ('1','2','3','4','5','6','7','8','9','10')
    start_probability = {'alpha': 0.2, 'beta': 0.1, 'alpha_and_beta':0.3, 'alpha_plus_beta': 0.4}
    if k == 6:
        transition_probability = t_p
        emission_probability = e_p
    elif k == 4:
        transition_probability = t_p2
        emission_probability = e_p2
    current_df = prot_struct_aa.values[0:full_set]
    count = 0
    matched = 0
    for row in current_df[t_set:]:
        tup = tuple([tt[1] for tt in observe_HMM(current_df[0:t_set],row[:n],k)])
        #results2 = viterbi(tup, states, start_probability, transition_probability, emission_probability)
        #def viterbi_display(obs, states, start_p, trans_p, emit_p):
        V = [{}]
        for st in states:
            V[0][st] = {"prob": start_probability[st] * emission_probability[st][tup[0]]['1'], "prev": None}
        # Run Viterbi when t > 0
        for t in range(1, len(tup)):
            V.append({})
            for st in states:
                max_tr_prob = max(V[t-1][prev_st]["prob"]*transition_probability[prev_st][st] for prev_st in states)
                for prev_st in states:
                    if V[t-1][prev_st]["prob"] * transition_probability[prev_st][st] == max_tr_prob:
                        max_prob = max_tr_prob * emission_probability[st][tup[t]][str(t+1)]
                        V[t][st] = {"prob": max_prob, "prev": prev_st}
                        break
        opt = []
        # The highest probability
        max_prob = max(value["prob"] for value in V[-1].values())
        previous = None
        # Get most probable state and its backtrack
        for st, data in V[-1].items():
            if data["prob"] == max_prob:
                opt.append(st)
                previous = st
                break
        # Follow the backtrack till the first observation
        for t in range(len(V) - 2, -1, -1):
            opt.insert(0, V[t + 1][previous]["prev"])
            previous = V[t + 1][previous]["prev"]
        #print 'The steps of states are ' + ' '.join(opt) + ' with highest probability of %s' % max_prob
        count += 1
        if opt[-1] == row[-1]:
            matched += 1
        succes_rate = float(matched)/float(count)
    results.append(succes_rate)
    print(k,results)
########################################################################################################################
#
#(k=6, [0.5423728813559322])
#(k=4, [0.5169491525423728])
########################################################################################################################
# Random Test
count = 0
matched = 0
for i in current_df[t_set:]:
    count += 1
    if states[randint(0,3)] == row[-1]:
        matched += 1
succes_rate = float(matched)/float(count)
print(succes_rate)
########################################################################################################################
########################################################################################################################
# HMM No cross validation. Just train on all test set and check on testing set
n = 20
total_cv_results = []
k = 10
t_set = 550
full_set = len(prot_struct_aa)
best_from_CV = 4
emmision_prop_calc = ['count','distance']
epc = emmision_prop_calc[0]
for k in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]:
    results = []
    np.random.seed(10)
    for ite in range(10):
        states = ('alpha','beta','alpha_and_beta','alpha_plus_beta')
        # We will do it for k - 10 to test the code, then we will do cross validation
        observations = ('1','2','3','4','5','6','7','8','9','10')
        # we know how many of each stae are in the test set, this knowledge could bias the HMM,
        # but we can say that there are more alpha beta combinations than alpha and betas alone
        start_probability = {'alpha': 0.2, 'beta': 0.1, 'alpha_and_beta':0.3, 'alpha_plus_beta': 0.4}
        # Transition and emission probabilities need to be calculated from 412 structures first, then 550
        # They are initialized to 0
        transition_probability = {'alpha' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0},\
                                  'beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0},\
                                  'alpha_and_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0},\
                                  'alpha_plus_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0}}
        emission_probability = {'alpha' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'alpha_and_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'alpha_plus_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0}}
        emission_counter = {'alpha' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'alpha_and_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'alpha_plus_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0}}
        emission_max = {'alpha' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'alpha_and_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0},\
                                'alpha_plus_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta': 0.0, 'alpha_plus_beta': 0.0}}    
        l = np.random.choice(range(full_set), full_set, replace=False)
        current_df = prot_struct_aa.values[l]
        traing_set_for_HMMmodel = np.array([train_HMM(current_df[0:t_set],row[:n+1],k) for row in current_df[0:t_set]])
        for i  in traing_set_for_HMMmodel:
            state = i[0][1]
            for j in i[1:]:
                transition_probability[state][j[1]] += 1
                if epc == 'count':
                    emission_probability[state][j[1]] += 1
                elif epc == 'distance':
                    emission_probability[state][j[1]] += float(j[0])
                    emission_counter[state][j[1]] += 1
                    if float(j[0]) > emission_max[state][j[1]]:
                        emission_max[state][j[1]] = float(j[0])
                state = j[1]        
        # Normalization emission_probabilities must add up to observations for training 412x10 = 4120
        for i in states:
            if epc == 'count':
                suma = 0.0
                for j in states:
                    suma += emission_probability[i][j]
                for j in states:
                    emission_probability[i][j] = emission_probability[i][j]/suma
            elif epc == 'distance':
                for j in states:
                    emission_probability[i][j] = emission_probability[i][j]/emission_counter[i][j]
                m = max(emission_probability[i].values())
                for j in states:
                    emission_probability[i][j] = 1 - emission_probability[i][j]/m
                s = sum(emission_probability[i].values())
                for j in states:
                    emission_probability[i][j] = emission_probability[i][j]/s                           
        for i in states:
            suma = 0.0
            for j in states:
                suma += transition_probability[i][j]
            for j in states:
                transition_probability[i][j] = transition_probability[i][j]/suma
        count = 0
        matched = 0
        for row in current_df[t_set:]:
            tup = tuple([tt[1] for tt in observe_HMM(current_df[0:t_set],row[:n],k)])
            results2 = viterbi(tup, states, start_probability, transition_probability, emission_probability)
            count += 1
            if results2[-1] == row[-1]:
                matched += 1
            succes_rate = float(matched)/float(count)
        results.append(succes_rate)
        if ite == 0 and k == 1:
            t_p = transition_probability
            e_p = emission_probability
            s_r = succes_rate
            b_k = k
        elif ite == 0 and k == best_from_CV:
            t_p2 = transition_probability
            e_p2 = emission_probability
            s_r2 = succes_rate
            b_k2 = k            
        else:
            if succes_rate > s_r:
                t_p = transition_probability
                e_p = emission_probability
                s_r = succes_rate
                b_k = k
            if k == best_from_CV and succes_rate > s_r2:
                t_p2 = transition_probability
                e_p2 = emission_probability
                s_r2 = succes_rate
                b_k2 = k
    print(k,results)
    All_results.loc[k,'KNN_HMM_C2_t'] = sum(results)/10
########################################################################################################################
# order2 HMM model
n = 20
k = 10
t_set = 550
full_set = len(prot_struct_aa)
emmision_prop_calc = ['count','distance','order','order2','consensus']
epc = emmision_prop_calc[3]
for k in range(1,21):
    results = []
    np.random.seed(10)
    for ite in range(10):    
        states = ('alpha','beta','alpha_and_beta','alpha_plus_beta')
        # We will do it for k - 10 to test the code, then we will do cross validation
        #observations = ('1','2','3','4','5','6','7','8','9','10')
        observations = tuple([str(i) for i in range(1,k+1)])
        # we know how many of each stae are in the test set, this knowledge could bias the HMM,
        # but we can say that there are more alpha beta combinations than alpha and betas alone
        start_probability = {'alpha': 0.2, 'beta': 0.1, 'alpha_and_beta':0.3, 'alpha_plus_beta': 0.4}
        # Transition and emission probabilities need to be calculated from 412 structures first, then 550
        # They are initialized to 0
        transition_probability = {'alpha' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0},\
                                  'beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0},\
                                  'alpha_and_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0},\
                                  'alpha_plus_beta' : {'alpha': 0.0, 'beta': 0.0, 'alpha_and_beta':0.0,'alpha_plus_beta':0.0}}
        emission_probability = {'alpha':{'alpha':{},'beta':{},'alpha_and_beta':{},'alpha_plus_beta':{}},\
                                'beta':{'alpha':{},'beta':{},'alpha_and_beta':{},'alpha_plus_beta':{}},\
                                'alpha_and_beta':{'alpha':{},'beta':{},'alpha_and_beta':{},'alpha_plus_beta':{}},\
                                'alpha_plus_beta':{'alpha':{},'beta':{},'alpha_and_beta':{},'alpha_plus_beta':{}}}
        emission_counter = {'alpha':{},'beta':{},'alpha_and_beta':{},'alpha_plus_beta':{}}
        emission_max = {'alpha':{},'beta':{},'alpha_and_beta':{},'alpha_plus_beta':{}}
        for i in states:
            for j in states:
                for ll in tuple([str(jj) for jj in range(1,k+1)]):
                    emission_probability[i][j][ll] = 0.0
            for j in tuple([str(jj) for jj in range(1,k+1)]):
                emission_counter[i][j] = 0.0
                emission_max[i][j] = 0.0
        l = np.random.choice(range(full_set), full_set, replace=False)
        current_df = prot_struct_aa.values[l]
        traing_set_for_HMMmodel = np.array([train_HMM(current_df[0:t_set],row[:n+1],k) for row in current_df[0:t_set]])
        cnt = 0
        for i  in traing_set_for_HMMmodel:
            #if cnt == 549:
            state = i[0][1]
            count_e = 1
            for j in i[1:]:
                transition_probability[state][j[1]] += 1
                if epc == 'count':
                    emission_probability[state][j[1]] += 1
                elif epc == 'distance':
                    emission_probability[state][j[1]] += float(j[0])
                    emission_counter[state][j[1]] += 1
                    if float(j[0]) > emission_max[state][j[1]]:
                        emission_max[state][j[1]] = float(j[0])
                elif epc == 'order':
                    #print(cnt,'order',count_e,i[0][1],j[1],emission_probability[i[0][1]][j[1]][str(count_e)])
                    emission_counter[i[0][1]][str(count_e)] += 1
                    emission_probability[i[0][1]][j[1]][str(count_e)] += 1
                    #print('    ',cnt,'order',count_e,i[0][1],j[1],emission_probability[i[0][1]][j[1]][str(count_e)])
                    count_e += 1
                elif epc == 'order2':
                    #print(cnt,'order',count_e,i[0][1],j[1],emission_probability[i[0][1]][j[1]][str(count_e)])
                    emission_counter[i[0][1]][str(count_e)] += 1
                    emission_probability[i[0][1]][j[1]][str(count_e)] += 1
                    #print('    ',cnt,'order',count_e,i[0][1],j[1],emission_probability[i[0][1]][j[1]][str(count_e)])
                    count_e += 1
                state = j[1]
            #cnt += 1
        # Normalization emission_probabilities must add up to observations for training 412x10 = 4120
        for i in states:
            if epc == 'count':
                suma = 0.0
                for j in states:
                    suma += emission_probability[i][j]
                for j in states:
                    emission_probability[i][j] = emission_probability[i][j]/suma
            elif epc == 'distance':
                for j in states:
                    emission_probability[i][j] = emission_probability[i][j]/emission_counter[i][j]
                m = max(emission_probability[i].values())
                for j in states:
                    emission_probability[i][j] = 1 - emission_probability[i][j]/m
                s = sum(emission_probability[i].values())
                for j in states:
                    emission_probability[i][j] = emission_probability[i][j]/s
            elif epc == 'order':
                for ii in states:
                    suma = 0.0
                    for j in tuple([str(jj) for jj in range(1,k+1)]):
                        suma += emission_probability[i][ii][j]
                    for j in tuple([str(jj) for jj in range(1,k+1)]):
                        emission_probability[i][ii][j] = emission_probability[i][ii][j]/suma
            elif epc == 'order2':                
                for j in tuple([str(jj) for jj in range(1,k+1)]):
                    suma = 0.0
                    for ii in states:
                        suma += emission_probability[i][ii][j]
                    for ii in states:
                        emission_probability[i][ii][j] = emission_probability[i][ii][j]/suma                         
        for i in states:
            suma = 0.0
            for j in states:
                suma += transition_probability[i][j]
            for j in states:
                transition_probability[i][j] = transition_probability[i][j]/suma
        
        count = 0
        matched = 0
        for row in current_df[t_set:]:
            tup = tuple([tt[1] for tt in observe_HMM(current_df[0:t_set],row[:n],k)])
            #results = viterbi(tup, states, start_probability, transition_probability, emission_probability)
            #def viterbi_display(obs, states, start_p, trans_p, emit_p):
            V = [{}]
            for st in states:
                V[0][st] = {"prob": start_probability[st] * emission_probability[st][tup[0]]['1'], "prev": None}
            # Run Viterbi when t > 0
            for t in range(1, len(tup)):
                V.append({})
                for st in states:
                    max_tr_prob = max(V[t-1][prev_st]["prob"]*transition_probability[prev_st][st] for prev_st in states)
                    for prev_st in states:
                        if V[t-1][prev_st]["prob"] * transition_probability[prev_st][st] == max_tr_prob:
                            max_prob = max_tr_prob * emission_probability[st][tup[t]][str(t+1)]
                            V[t][st] = {"prob": max_prob, "prev": prev_st}
                            break
            opt = []
            # The highest probability
            max_prob = max(value["prob"] for value in V[-1].values())
            previous = None
            # Get most probable state and its backtrack
            for st, data in V[-1].items():
                if data["prob"] == max_prob:
                    opt.append(st)
                    previous = st
                    break
            # Follow the backtrack till the first observation
            for t in range(len(V) - 2, -1, -1):
                opt.insert(0, V[t + 1][previous]["prev"])
                previous = V[t + 1][previous]["prev"]
        
            #print 'The steps of states are ' + ' '.join(opt) + ' with highest probability of %s' % max_prob
            count += 1
            if opt[-1] == row[-1]:
                matched += 1
            succes_rate = float(matched)/float(count)
        results.append(succes_rate)
    print(k,results)
    All_results.loc[k,'KNN_HMM_O2_t'] = sum(results)/10
'''
       KNN_cv      KNN_t KNN_HMM_C_vc KNN_HMM_C_t KNN_HMM_O_cv KNN_HMM_O_t  1   0.5282609  0.5694915    0.5202899   0.5677966    0.5202899   0.5677966   
2   0.5123188  0.5466102    0.5028986   0.5432203    0.5036232   0.5449153   
3    0.523913  0.5686441    0.5275362   0.5347458    0.5086957   0.5211864   
4   0.5688406  0.5881356    0.5384058   0.5423729    0.5195652   0.5152542   
5   0.5702899  0.5915254    0.5318841   0.5381356    0.5369565   0.5381356   
6    0.592029  0.6067797    0.5224638   0.5516949    0.5202899   0.5381356   
7   0.5855072        0.6    0.4949275   0.5313559    0.5057971    0.520339   
8   0.5884058  0.6110169    0.5094203   0.5271186    0.4833333   0.5152542   
9   0.5963768  0.6059322    0.5137681   0.5228814    0.5057971   0.4966102   
10  0.6188406  0.6135593    0.5137681   0.5169492    0.5130435    0.490678   
11  0.6123188  0.6127119    0.5108696    0.529661    0.4949275   0.5084746   
12  0.6094203  0.6110169    0.5086957    0.529661    0.4978261   0.5177966   
13  0.6115942  0.6050847    0.5101449   0.5076271    0.4949275   0.4898305   
14  0.6101449  0.6050847    0.4992754   0.5220339     0.473913   0.4983051   
15  0.6101449  0.6016949    0.5014493   0.5076271    0.4797101   0.4855932   
16  0.6094203  0.6084746     0.484058   0.4983051    0.4775362   0.4711864   
17  0.6065217  0.6084746    0.4963768   0.5169492    0.4891304   0.4737288   
18  0.6094203  0.6033898    0.4652174   0.5076271    0.4804348   0.4677966   
19  0.6036232        0.6    0.4789855         0.5    0.4847826   0.4567797   
20  0.6036232  0.5940678    0.4782609   0.5016949    0.4746377   0.4838983   

    KNN_HMM_C2_t  KNN_HMM_O2_t  
1       0.564407      0.564407  
2       0.507627      0.555085  
3       0.539831      0.551695  
4       0.527119      0.556780  
5       0.548305      0.541525  
6       0.553390      0.543220  
7       0.538983      0.554237  
8       0.536441      0.552542  
9       0.530508      0.538136  
10      0.521186      0.513559  
11      0.530508      0.522881  
12      0.538136      0.521186  
13      0.507627      0.492373  
14      0.520339      0.505932  
15      0.511864      0.500000  
16      0.502542      0.500000  
17      0.509322      0.501695  
18      0.511864      0.495763  
19      0.489831      0.470339  
20      0.500000      0.489831 
'''
