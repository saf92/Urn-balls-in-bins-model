
from urn_packages import *

'''
If length of array is to a power of 2 returns the array, if not appends zeros to the array so it's length is
to the nearest power of 2 greater than the length of the array
'''

def get_new_div_two_array(x):
    x1=len(x)
    x2=np.log2(x1)
    if x2%1==0.0:
        return x
    x3=np.int(x2)+1
    x4=2**x3-x1
    return np.concatenate((x,np.zeros(x4)))

'''
Creates binary tree of rates array
'''

def r_tree(rates): 
    rates=get_new_div_two_array(rates)
    N=len(rates) #len(rates)=N=2^M
    M=int(np.log2(N))
    T=2**(M+1)
    R_tree = np.zeros(T)
    for i in range(N):
        R_tree[N+i]=rates[i];
        
    for i in range(N-1,0,-1):
        R_tree[i]= R_tree[i<<1]+R_tree[(i<<1)+1]
    
    return R_tree

'''
Finds the tree and rates position of the random value p (such that 0<=p<=sum(rates)) 
'''

def find_position(R_tree,p):
    l=len(R_tree)
    M=int(np.log2(l)-1)
    N=2**M
    tree_pos=1
    check=R_tree[tree_pos<<1]
    
    for i in range(M):
        if p<check:
            tree_pos=tree_pos<<1 #go left
            if tree_pos<N:
                check=R_tree[tree_pos<<1]
        else:
            p=p-check
            tree_pos=(tree_pos<<1)+1 #go right
            if tree_pos<N:
                check=R_tree[tree_pos<<1]
    rates_pos=tree_pos-N
    return tree_pos,rates_pos

'''
Updates the binary tree for new rates vector with delta added to one of the rates
'''

def r_tree_update(R_tree, tree_pos, delta):
    l=len(R_tree)
    k=int(np.log2(l))
    for i in range(k):
        R_tree[tree_pos]=R_tree[tree_pos]+delta
        tree_pos=tree_pos>>1
    return R_tree
