from urn_packages import *
from binary_tree_functions import *

'''
Functions below run balls in bins process or non-linear urn model, with
feedback function f(w)=aw^g.

Parameters:
w0: vector of initial wealth,
w_add: wealth packet to be added at each simulation
g: gamma value in urn model
T: number of iterations
'''

def nl_Polya_urn_stand(w0,w_add,g,T):
    N=len(w0)
    individuals=np.arange(0,N,1)
    w=w0.copy()
    for t in range(T):
        powers=w**g
        probs=powers/sum(powers)
        i=np.random.choice(individuals, 1, p=probs)[0]
        w[i]+=w_add
    return w


def nl_Polya_urn_stand_bt(w0,w_add,g,T):
    N=len(w0)
    powers = w0**g
    R_tree = r_tree(powers)
    w=w0.copy()
    for t in range(T):
        s=R_tree[1]
        powers=w**g
        p=np.random.uniform(0,s)
        tree_pos,rates_pos=find_position(R_tree,p)
        w[rates_pos]+=w_add
        delta=w[rates_pos]**g-(w[rates_pos]-w_add)**g
        R_tree = r_tree_update(R_tree, tree_pos, delta)
    return w

#Urn simulation with binary search but with fitness as part of feedback function
def nl_Polya_urn_stand_bt_fit(w0,w_add,a,g,T):
    N=len(w0)
    powers = a*w0**g
    R_tree = r_tree(powers)
    w=w0.copy()
    for t in range(T):
        s=R_tree[1]
        powers=a*w**g
        p=np.random.uniform(0,s)
        tree_pos,rates_pos=find_position(R_tree,p)
        w[rates_pos]+=w_add
        delta=a[rates_pos]*(w[rates_pos]**g-(w[rates_pos]-w_add)**g)
        R_tree = r_tree_update(R_tree, tree_pos, delta)
    return w

#####################################################################################################

#Run urn model ns times and append each run to an array

def nl_Polya_urn_mult_times(w0,w_add,g,T,ns):
    Ws=[]
    for i in range(ns):
        W=nl_Polya_urn_stand_bt(w0,w_add,g,T)
        Ws.append(W)
        print(i)
    return Ws

#Produces new vector with largest value removed

def get_losing(w):
    w1=np.sort(w)
    return w1[:-1]


#################################################################################################

#Predict wealth

#g=1
#Function for expected wealth when g=1
def expected_wealth_fit_pred_g1(t,w0,a):
    return w0*np.exp(a*t)

#Function for different
def dif_sum_func_g1(t,total,w0,a):
    total1=np.sum(expected_wealth_fit_pred_g1(t,w0,a))
    return total-total1

def get_time_g1(init_t,total,w0,a):
    func = lambda t : total-np.sum(expected_wealth_fit_pred_g1(t,w0,a))
    return(fsolve(func,init_t)[0])

#g>1

def expected_wealth_fit_pred(t,w0,g,a):
    k=1/w0**(g-1)-(g-1)*a*t
    k1=k**(-1/(g-1)) #g>1
    return k1

def dif_sum_func(t,total,w0,g,a):
    total1=np.sum(expected_wealth_fit_pred(t,w0,g,a))
    return total-total1

def get_time(init_t,total,w0,g,a):
    func = lambda t : total-np.sum(expected_wealth_fit_pred(t,w0,g,a))
    return(fsolve(func,init_t)[0])

########################################

def find_t_g1(t_tests, w_fit, w0, fitness):
    t=0
    i=0
    d=dif_sum_func_g1(t_tests[i],np.sum(w_fit),w0,fitness)
    while (d>0):
        i+=1
        d=dif_sum_func_g1(t_tests[i],np.sum(w_fit),w0,fitness)
    return t_tests[i-1]

def find_t_g(t_tests, w_fit, w0, g, fitness):
    t=0
    i=0
    d=dif_sum_func(t_tests[i],np.sum(w_fit),w0,g,fitness)
    while (d>0):
        i+=1
        d=dif_sum_func(t_tests[i],np.sum(w_fit),w0,g,fitness)
    return t_tests[i-1]

###########################################################################################

#Sum formula

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
