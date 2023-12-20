#ultimate a value combiner and what not.
#rn - reads in from files e1, e2, e3 etc. howver will eventually make it grasp.outs
import numpy as np
num_levels = 40
#print(e1_raw) 
transition_matrix = np.zeros([num_levels,num_levels])

def readin_and_add(transition_matrix,filepath):
    raw = np.loadtxt(filepath,usecols=[0,1,2,3])
    lower = raw[:,0].astype(int) - 1 #subtracting as python is 0 indexed, god save us
    upper = raw[:,1].astype(int) - 1
    avalues = raw[:,3]

    for kk in range(0,len(avalues)):
        i = lower[kk]
        j = upper[kk] 
        avalue = avalues[kk]
        transition_matrix[i,j]+= avalue
    return transition_matrix 

def finalise(transition_matrix):
    si = np.shape(transition_matrix)
    for jj in range(0,si[0]):
        for ii in range(0,si[1]):
            if transition_matrix[jj,ii] == 0.0:
                transition_matrix[jj,ii] = 1.0e-30
    return transition_matrix

def output(transition_matrix,output_file,num_levels):
    #f = open(output_file,'w')

    for ii in range(0,num_levels):
        for jj in range(ii+1,num_levels):
            lower = ii+1#python is 0 indexed. god save us.
            upper = jj+1
            avalue = transition_matrix[ii,jj]
            #print(upper,lower,avalue)
            string_to_be_written = '{:5}'.format(upper)
            string_to_be_written += '{:5}'.format(lower)
            avalue_to_be_written = round(avalue,ndigits=2)
            string_to_be_written += '{:}'.format(" ")
            avalue_string = '{:.2E}'.format(avalue)
            avalue_string = avalue_string.replace("E","")
            string_to_be_written += avalue_string
            print(string_to_be_written)


    return 0

transition_matrix = readin_and_add(transition_matrix,'e1') 
#output(transition_matrix,'we',num_levels=10)
transition_matrix = readin_and_add(transition_matrix,'e2') 
#output(transition_matrix,'we',num_levels=10)
transition_matrix = readin_and_add(transition_matrix,'m1') 
transition_matrix = readin_and_add(transition_matrix,'m2') 
transition_matrix = finalise(transition_matrix)
output(transition_matrix,'we',num_levels=10)