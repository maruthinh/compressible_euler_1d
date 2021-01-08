#to plot results
from numpy import*
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import pylab

#filename1='output.dat'
filename1='ResultTest1.dat'
data1=loadtxt(filename1)
#filename2='prob0.out' #problem 0., shock tube without sonic point
#filename2='woodward_collela.dat' #problem 3
#filename2='shock_collision.dat' #problem 4
#filename2='over_heating.dat' #problem 2
filename2='shock_tube.dat' #problem 1
#filename2='prob5.dat' #problem 5
#filename2='steadycontact.dat' #problem 6 steady contact
#filename2='rand_choice_prob10_026.out' #problem 10 woodward and collela problem for t=0.026
#filename2='rand_choice_prob10_038.out' #problem 10 woodward and collela problem for t=0.038
#filename2='slow_movin_shock_t_0.15.dat' #slowly moving shock at t=0.15 prob 9
data=loadtxt(filename2)

#to plot the figures
x=data1[:,0]
y=data1[:,1]
plt.figure(1)
plt.subplot(221)
plt.plot(x,y,'--')
plt.hold
x1=data[:,0]
y1=data[:,1]
plt.plot(x1,y1)
plt.xlabel('x')
plt.ylabel(r'$\rho$', fontsize=18)

#to plot second subplot
x=data1[:,0]
y=data1[:,2]
plt.subplot(222)
plt.plot(x,y,'--')
plt.hold
data=loadtxt(filename2)
x1=data[:,0]
y1=data[:,2]
plt.plot(x1,y1)
plt.xlabel('x')
plt.ylabel(r'u', fontsize=18)

#to plot third subplot
x=data1[:,0]
y=data1[:,3]
plt.subplot(223)
plt.plot(x,y,'--')
plt.hold
data=loadtxt(filename2)
x1=data[:,0]
y1=data[:,3]
plt.plot(x1,y1)
plt.xlabel('x')
plt.ylabel(r'p', fontsize=18)

#to plot fourth subplot
x=data1[:,0]
y=data1[:,4]
plt.subplot(224)
plt.plot(x,y,'--')
plt.hold
data=loadtxt(filename2)
x1=data[:,0]
y1=data[:,4]
plt.plot(x1,y1)
plt.xlabel('x')
plt.ylabel(r'e', fontsize=18)


#
#f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
#ax1.plot(data[:,0], data[:,1])
#ax1.plot(data1[:,0], data1[:,1])
#ax1.set_title('Sharing x per column, y per row')
#ax2.scatter(x, y)
#ax3.scatter(x, 2 * y ** 2 - 1, color='r')
#ax4.plot(x, 2 * y ** 2 - 1, color='r')

#plt.show()
plt.savefig("Test1.pdf",bbox_inches='tight')
#save("signal", ext="svg", close=True, verbose=True)
