#to plot results
from numpy import*
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter


#plt.figure([2,2])
#plt.figure(figsize=(10,8),facecolor='w') 
#filename1='output.dat'
filename1='Res_test_case5.txt'
data1=loadtxt(filename1)
if (filename1=='Res_test_case1.txt'):
    filename2='shock_tube.dat' #problem 1
elif (filename1=='Res_test_case2.txt'):
    filename2='over_heating.dat' #problem 2    
elif (filename1=='Res_test_case3.txt'):
    filename2='woodward_collela.dat' #problem 3
elif (filename1=='Res_test_case4.txt'):    
    filename2='shock_collision.dat' #problem 4    
elif (filename1=='Res_test_case5.txt'):    
    filename2='prob5.dat' #problem 5
elif (filename1=='Res_test_case6.txt'):    
    filename2='steadycontact.dat' #problem 6 steady contact
elif (filename1=='Res_test_case7.txt'):    
    filename2='ExactProb7.dat' #problem 7 steady shock
elif (filename1=='Res_test_case8.txt'):    
    filename2='slow_movin_shock_t_0.15.dat' #slowly moving shock at t=0.15 prob 8
elif (filename1=='Res_test_case9.txt'):    
    filename2='prob9.out' #problem 9
elif (filename1=='Res_test_case10.txt'):    
    filename2='rand_choice_prob10_038.out' #problem 10 woodward and collela problem for t=0.038    
    #filename2='rand_choice_prob10_026.out' #problem 10 woodward and collela problem for t=0.026

data=loadtxt(filename2)
legend1='KFDS+'
###############################################################################
###############################################################################
#replace '.txt' from filename1 to empty string. This string will be used to save the 
#results
FileToSave=filename1.replace(".txt", "")
###############################################################################
############################################################################### 
#to plot the figures
plt.figure(1)
plt.plot(data1[:,0], data1[:,1],'--', color='k', markersize=10, linewidth=3)
plt.plot(data[:,0], data[:,1], color='k', markersize=10, linewidth=3)
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
plt.xlabel(r'x', fontsize=18)
plt.ylabel(r'$\rho$', fontsize=18)
#plt.legend(['KFDS_MOVERS','Exact'], fontsize=18, loc='0', frameon=False)
plt.legend([legend1,'Exact'], fontsize=18, loc='0', frameon=False)
plt.savefig(FileToSave+'_den'+'.pdf', bbox_inches='tight')
###############################################################################
###############################################################################
plt.figure(2)
plt.plot(data1[:,0], data1[:,2],'--', color='k', markersize=10, linewidth=3)
plt.plot(data[:,0], data[:,2], color='k', markersize=10, linewidth=3)
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
plt.xlabel(r'x', fontsize=18)
plt.ylabel(r'$u$', fontsize=18)
#plt.legend(['KFDS_MOVERS','Exact'], fontsize=18, loc='0', frameon=False)
plt.legend([legend1,'Exact'], fontsize=18, loc='0', frameon=False)
plt.savefig(FileToSave+'_u'+'.pdf',bbox_inches='tight')
###############################################################################
###############################################################################
plt.figure(3)
plt.plot(data1[:,0], data1[:,3],'--', color='k', markersize=10, linewidth=3)
plt.plot(data[:,0], data[:,3], color='k', markersize=10, linewidth=3)
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
plt.xlabel(r'x', fontsize=18)
plt.ylabel(r'$p$', fontsize=18)
#plt.legend(['KFDS_MOVERS','Exact'], fontsize=18, loc='0', frameon=False)
plt.legend([legend1,'Exact'], fontsize=18, loc='0', frameon=False)
plt.savefig(FileToSave+'_p'+'.pdf',bbox_inches='tight')
###############################################################################
###############################################################################
plt.figure(4)
plt.plot(data1[:,0], data1[:,4],'--', color='k', markersize=10, linewidth=3)
plt.plot(data[:,0], data[:,4], color='k', markersize=10, linewidth=3)
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
plt.xlabel(r'x', fontsize=18)
plt.ylabel(r'$e$', fontsize=18)
#plt.legend(['KFDS_MOVERS','Exact'], fontsize=18, loc='0', frameon=False)
plt.legend([legend1,'Exact'], fontsize=18, loc='0', frameon=False)
plt.savefig(FileToSave+'_e'+'.pdf',bbox_inches='tight')
###############################################################################
###############################################################################
#plt.axis([0, 1, data[0,1]-1, max(data[:,1])+1.0])
#plt.axis([0,1.1,0.9,1.5])
##to plot entropy
#Gamma=1.4
#s1=data[0:,3]/data[0:,1]**Gamma
#s=numpy.zeros(len(s1))
#for i in range(len(s1)):
#    s[i]=math.log(s1[i])
#plt.plot(x1,s)
#plt.show()


###############################################################################
###################################END#########################################
###############################################################################
