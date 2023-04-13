
from urn_packages import *

#Finds emperical distribution of data
def ecdf_points(data):
    L=len(data)
    data_array=[]
    for i in range(L):
        data_array.append(data[i])
    np.asarray(data_array,float)
    x=np.sort(data_array)
    y=np.arange(0,1,1/L)
    X=[x[0]]
    Y=[y[0]]
    for i in range(0,L-1):
        if(x[i]!=x[i+1]):
            X.append(x[i+1])
            Y.append(y[i+1])
    X=np.asarray(X,float)
    Y=np.asarray(Y,float)
    return X,Y

#Finds tail of data
def tail(data):
    x,y=ecdf_points(data)
    return x, 1-y


#######################################################################

#Finds index of array with value of that index closest to the value imputted
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#Power law function
def power_law(x,a,b):
    return a*x**b
###########################################################################

#Fit exponential tail with MLE
def exp_tail(x,l):
    return np.exp(-l*x)

def mle_exp_lambda(x_sample):
    n=len(x_sample)
    return n/(np.sum(x_sample))

###########################################################################

'''Functions from Clauset et al paper'''

#Get's the MLE power beta parameter
def mle_pl_beta(x,x_min):
    n=len(x)
    return n/np.sum(np.log(x/x_min))

#Fit's the power law with MLE outputting parameters, predictions and index for inputted
#minimum x value for sample fit
def get_power_law_fit(data,x_min):
    data_s=np.sort(data)
    x,y=tail(data_s)
    ind=find_nearest(data_s,x_min)
    ind1=find_nearest(x,x_min)
    x_sample=data_s[ind:]
    x1=x[ind1:]
    b=mle_pl_beta(x_sample,x_min)
    a=y[ind1]*x_min**b
    y_pred=power_law(x1,a,-b)
    return x1,y_pred,a,b


#KS statistic
def KS_stat(y,y_pred):
    return np.max(np.abs(y-y_pred))

#Get's estimate for x_min. Outputs this estimate which is the minimum of the
#KS-statistics also outputted
def x_min_pred(data):
    es1=[]
    x,y=tail(data)
    for i in range(len(x)-2):
        xs,ys,a,b=get_power_law_fit(data,x[i])
        y_norm=y[i]
        yn=y[i:]/y_norm
        ysn=ys/y_norm
        e1=KS_stat(yn,ysn)
        es1.append(e1)
    ind1=np.argmin(es1)
    x_min_pred=x[ind1]
    return es1,x_min_pred

############################################################################
