from urn_packages import *
from binary_tree_functions import *

'''
Functions below run balls in bins process or non-linear urn model, with
feedback function f(w)=aw^g.

Parameters:
w0: vector of initial wealth,
w_add: wealth packet to be added at each simulation
a: fitness value
g: gamma value
iters: number of iterations
'''

#No binary tree - slower
def nl_Polya_urn_stand(w0,w_add,g,iters):
    N=len(w0)
    individuals=np.arange(0,N,1)
    w=w0.copy()
    for t in range(iters):
        powers=w**g
        probs=powers/sum(powers)
        i=np.random.choice(individuals, 1, p=probs)[0]
        w[i]+=w_add
    return w

#With binary tree
def nl_Polya_urn_stand_bt(w0,w_add,g,iters):
    N=len(w0)
    powers = w0**g
    R_tree = r_tree(powers)
    w=w0.copy()
    for t in range(iters):
        s=R_tree[1]
        powers=w**g
        p=np.random.uniform(0,s)
        tree_pos,rates_pos=find_position(R_tree,p)
        w[rates_pos]+=w_add
        delta=w[rates_pos]**g-(w[rates_pos]-w_add)**g
        R_tree = r_tree_update(R_tree, tree_pos, delta)
    return w

#Urn simulation with binary search but with fitness as part of feedback function
def nl_Polya_urn_stand_bt_fit(w0,w_add,a,g,iters):
    N=len(w0)
    powers = a*w0**g
    R_tree = r_tree(powers)
    w=w0.copy()
    for t in range(iters):
        s=R_tree[1]
        powers=a*w**g
        p=np.random.uniform(0,s)
        tree_pos,rates_pos=find_position(R_tree,p)
        w[rates_pos]+=w_add
        delta=a[rates_pos]*(w[rates_pos]**g-(w[rates_pos]-w_add)**g)
        R_tree = r_tree_update(R_tree, tree_pos, delta)
    return w

#### Pure birth process that through exponential embedding produces the same model

#Single run of birth process
def pure_birth_process_single(w,w_add,a,g,jt_max,w_max):
    jump_time=0
    while jump_time<jt_max and w<w_max:
        f=a*w**(-g) #feedback fn
        holding_time=expon.rvs(scale=f)
        jump_time+=holding_time
        w+=w_add
    return w

#Jump time of single run of birth process
def pure_birth_process_single_t(w,w_add,a,g,w_max):
    jump_time=0
    while w<w_max:
        f=a*w**(-g) #feedback fn
        holding_time=expon.rvs(scale=f)
        jump_time+=holding_time
        w+=w_add
    return jump_time

#Mean of jump times
def jump_time_runs(runs,w_add,a,g,w_max):
    ts=[]
    for i in range(runs):
        w=1
        t=pure_birth_process_single_t(w,w_add,a,g,w_max)
        ts.append(t)
    return np.mean(ts)

#Lower bound for expected explosion time g>1
def exp_expl_time(g):
    return 1/(g-1)

#####################################################################################################

#Run urn model ns times and append each run to an array
def nl_Polya_urn_mult_times(w0,w_add,g,T,ns):
    Ws=[]
    for i in range(ns):
        W=nl_Polya_urn_stand_bt(w0,w_add,g,T)
        Ws.append(W)
    return Ws

#Produces new vector with largest value removed
def get_losing(w):
    w1=np.sort(w)
    return w1[:-1]


#################################################################################################

#Predict wealth

#g=1 ########################################
#Function for expected wealth when g=1
def expected_wealth_fit_pred_g1(t,w0,a):
    return w0*np.exp(a*t)

#Function for difference in sum of expected wealth and sum of simulated values
def dif_sum_func_g1(t,total,w0,a):
    total1=np.sum(expected_wealth_fit_pred_g1(t,w0,a))
    return total-total1

#Solve for t to minimise above sum
def get_time_g1(init_t,total,w0,a):
    func = lambda t : total-np.sum(expected_wealth_fit_pred_g1(t,w0,a))
    return(fsolve(func,init_t)[0])
#####################################################################

#g>1 ###################################################################
#Function for expected wealth when g>1
def expected_wealth_fit_pred(t,w0,g,a):
    k=1/w0**(g-1)-(g-1)*a*t
    k1=k**(-1/(g-1)) #g>1
    return k1

#Function for difference in sum of expected wealth and sum of simulated values
def dif_sum_func(t,total,w0,g,a):
    total1=np.sum(expected_wealth_fit_pred(t,w0,g,a))
    return total-total1

#Solve for t to minimise above sum
def get_time(init_t,total,w0,g,a):
    func = lambda t : total-np.sum(expected_wealth_fit_pred(t,w0,g,a))
    return(fsolve(func,init_t)[0])
##############################################################################

########################################

#Find time to approximate expected wealth
#g=1
def find_t_g1(t_tests, w_fit, w0, fitness):
    t=0
    i=0
    d=dif_sum_func_g1(t_tests[i],np.sum(w_fit),w0,fitness)
    while (d>0):
        i+=1
        d=dif_sum_func_g1(t_tests[i],np.sum(w_fit),w0,fitness)
    return t_tests[i-1]

#g>1
def find_t_g(t_tests, w_fit, w0, g, fitness):
    t=0
    i=0
    d=dif_sum_func(t_tests[i],np.sum(w_fit),w0,g,fitness)
    while (d>0):
        i+=1
        d=dif_sum_func(t_tests[i],np.sum(w_fit),w0,g,fitness)
    return t_tests[i-1]

###########################################################################################

'''
Sum formula function for mass function p_t(\omega)

Inputs:
g: gamma
t: time t
w0: initial value of w
w_s: max w to go up to

Outputs:
w_v: vector of w
A: array of vector of coefficients for each w
p: prediction for p_t(\omega)
'''

def get_sum_p(g,t,w0,w_s):
    A=[[1]]
    p=[]
    s=0
    w_v=np.arange(w0,w_s+1,1)
    for w in range(w0+1,w_s+1):
        A_w=[]
        n=len(A[s])
        for i in range(1,n+1):
            a=(w-1)**g/(w**g-i**g)*A[s][i-1]
            A_w.append(a)
        A_w.append(-sum(A_w))
        A.append(A_w)
        s+=1
    for i in range(len(A)):
            j=np.arange(w0,w0+len(A[i]),1)
            p_s=np.sum(A[i]*np.exp(-(j**g)*t))
            p.append(p_s)
    return w_v,A,p

#Functions to approximate sum formula
def master_approx(w_0,w,g,t):
    return np.exp(-w_0**g*t)*(w_0/w)**g

def master_approx1(w_0,w,g,t):
    k=1
    ks=[k]
    for i in range(len(w)-1):
        k=k*1/(1-(w_0/w[i+1])**g)
        ks.append(k)
    ks=np.asarray(ks,float)
    return ks*np.exp(-w_0**g*t)*(w_0/w)**g
