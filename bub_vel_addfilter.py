# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 11:14:40 2017

@author: Xin
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 15:45:35 2017
you can use higher central difference method
the default is 2nd central scheme by numpy.gradient
@author: Xin
"""



import h5py
import numpy as np
#import numdifftools as nd
import pylab
import matplotlib.pyplot as plt

lz=3.2
variables = ['Prho', 'PVz']
h5file = h5py.File('tests_single_new.h5', 'r')

#h5file = h5py.File('temp.h5', 'r+')
#read dataset dimensions
mylist = ['Fields/','Prho','/','002000']
delimiter=''
filepath = delimiter.join(mylist)
databk = h5file.get(filepath)
m1 = np.array(databk)
nz=m1.shape[0]
ny=m1.shape[1]
nx=m1.shape[2]

dx=lz/nz
def high_order_gradient(fx,dx,order):
    length=len(fx)
    fxgrad = np.zeros(length)
    if order==4: 
        for i in range(2):
            fxgrad[i]=(-25*fx[i]+48*fx[i+1]-36*fx[i+2]+16*fx[i+3]-3*fx[i+4])/(12*dx)
        for i in range(2,length-2):
            fxgrad[i]=(-fx[i+2]+fx[i+1]*8-fx[i-1]*8+fx[i-2])/(12*dx)
        for i in range(length-2,length):
            fxgrad[i]=(25*fx[i]-48*fx[i-1]+36*fx[i-2]-16*fx[i-3]+3*fx[i-4])/(12*dx)
    if order==6:
        for i in range(3):
            fxgrad[i]=(-49/20*fx[i]+6*fx[i+1]-15/2*fx[i+2]+20/3*fx[i+3]-15/4*fx[i+4]+6/5*fx[i+5]-1/6*fx[i+6])/(dx)
        for i in range(3,length-3):
            fxgrad[i]=(fx[i+3]-9*fx[i+2]+45*fx[i+1]-45*fx[i-1]+9*fx[i-2]-fx[i-3])/(60*dx)
        for i in range(length-3,length):
            fxgrad[i]=(49/20*fx[i]-6*fx[i-1]+15/2*fx[i-2]-20/3*fx[i-3]+15/4*fx[i-4]-6/5*fx[i-5]+1/6*fx[i-6])/(dx)

        
        
        
        
    return fxgrad

step = []
for i in range(721):
    step.append(str((i+1)*1000).zfill(6))



bub_loc_all = np.zeros(len(step))
bub_loc_all_ori = np.zeros(len(step))
sp_loc_all = np.zeros(len(step))
bub_velo_all = np.zeros(len(step))
bub_velo_all_aver = np.zeros(len(step))
sp_velo_all = np.zeros(len(step))
bub_velo_all_ori = np.zeros(len(step))

ii=0
for istep in step:
    delimiter = ''
    mylist = ['Fields/', variables[0], '/', istep]
    filepath = delimiter.join(mylist)
    databk = h5file.get(filepath)
    np_data = np.array(databk)
    m1 = np_data[:, ny/2-1, 0]
    m2 = np_data[:, 0, 0]
    m1_filter=m1.copy();    
    m2_filter=m2.copy();    
    for jstep in range(2,nz-3):
        m1_filter[jstep]=(m1[jstep-2]+m1[jstep-1]+m1[jstep]+m1[jstep+1]+m1[jstep+2])/5;
        m2_filter[jstep]=(m2[jstep-2]+m2[jstep-1]+m2[jstep]+m2[jstep+1]+m2[jstep+2])/5;
    #    m2_filter[jstep]=(m2[jstep-1]+m2[jstep]+m2[jstep+1])/3;
    m1_grad = np.gradient(m1)
    m2_grad = np.gradient(m2)

    m1_grad = high_order_gradient(m1_filter,dx,6)
    m2_grad = high_order_gradient(m2_filter,dx,6)
#    m2_grad_ori = high_order_gradient(m2,dx,6)

    sp_loc = np.argmax(m1_grad)
    bub_loc = np.argmax(m2_grad)
    bub_loc_ori = np.argmax(m2_grad_ori)
  #  print 'at step ', ii, ' sp_loc is ',sp_loc, ' bub loc is ',bub_loc
    sp_loc_all[ii] = sp_loc
    bub_loc_all[ii] = bub_loc
 #   bub_loc_all_ori[ii] = bub_loc_ori


    mylist = ['Fields/', variables[1], '/', istep]
    filepath = delimiter.join(mylist)
    databk = h5file.get(filepath)
    np_data = np.array(databk)
    m1 = np_data[:, ny/2-1, 0]
    m2 = np_data[:, 0, 0]
    sp_velo = m1[sp_loc]
    bub_velo = m2[bub_loc]
 #   bub_velo_ori=m2[bub_loc_ori]
  #  print 'at step ', ii, ' sp velo is ',sp_velo ,' bub velo is ',bub_velo
    bub_velo_all[ii] = bub_velo
 #   bub_velo_all_ori[ii] = bub_velo_ori
 #   bub_velo_all_aver=bub_velo_all.copy()
    #calculate avearge buble velocity
#    for jstep in range(10,585-11):
        #bub_velo_all_aver[jstep]=(bub_velo_all[jstep-10]+bub_velo_all[jstep-9]+bub_velo_all[jstep-8]+bub_velo_all[jstep-7]+bub_velo_all[jstep-6]+bub_velo_all[jstep-5]+bub_velo_all[jstep-4]+bub_velo_all[jstep-3]+bub_velo_all[jstep-2]+bub_velo_all[jstep-1]+bub_velo_all[jstep]+bub_velo_all[jstep+1]+bub_velo_all[jstep+2]+bub_velo_all[jstep+3]+bub_velo_all[jstep+4]+bub_velo_all[jstep+5]+bub_velo_all[jstep+6]+bub_velo_all[jstep+7]+bub_velo_all[jstep+8]+bub_velo_all[jstep+9]+bub_velo_all[jstep+10])/21
#        bub_velo_all_aver[jstep]=(bub_velo_all_ori[jstep-2]+bub_velo_all_ori[jstep-1]+bub_velo_all_ori[jstep]+bub_velo_all_ori[jstep+1]+bub_velo_all_ori[jstep+2])/5
    sp_velo_all[ii] = sp_velo
    ii = ii + 1

all_data = np.column_stack((bub_loc_all,bub_velo_all,sp_loc_all,sp_velo_all))
np.savetxt('saved_bub_velo',all_data,delimiter='\t',fmt='%s')

#print('The spike locations are')
#print(sp_loc_all)
#print('The spike velo are')
#print(sp_velo_all)
#print('The bubble locations are')
#print(bub_loc_all)
#print('The bubble velo are')
#print(bub_velo_all)
#fig, ax = plt.subplots( nrows=1, ncols=1 )
#plt.plot(bub_velo_all_ori, label='original')
#plt.plot(bub_velo_all, label='average on density')
#plt.plot(bub_velo_all_aver, label='average on velocity')
#pylab.legend(loc='best')
#plt.plot(m1_filter)
#plt.plot(m1)
#plt.title('2nd order')

#plt.savefig('velo.eps', format='eps', dpi=1000)


