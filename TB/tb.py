################################################
##### Last update 01/2017  Felipe Crasto de Lima
################################################
import numpy as np
from pythtb import *
import pylab as plt
from matplotlib.ticker import AutoMinorLocator
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['axes.linewidth']=2

############# MODEL GENERATION PARAMETERS #######################
a=             # lattice parameter
lat=[[],[]]    # 2D lattice vectors
norb=          # number of orbitals
orb=[[],[],[]] # orbitals coordinates in terms of lattice vectors

# list of k-points nodes for band calculations in
# term of reciprocal vectors
path=[[],[],[],[]]
# labels of the nodes
label=(r'', r'', r'', r'')
# total number of interpolated k-points along the path
nk=

#####  TB parameters  #####
t=1.0    #normalization of NN-hopping   
alpha=   #parameter controling further neighbors hopping strength
###########################



kg=tb_model(2,2,lat,orb,per=None)

y=[0]*norb

kg.set_onsite(y)

lat=kg.get_lat()
norb=kg.get_num_orbitals()
xyzf=kg.get_orb()

print "number of orbitals:"
print norb
print 'lattice vectors:'
print lat
print 'orbitals fractional coordinates:'
print xyzf

dd = 0.01*a
dnn=7*a 
dnnn=10*a
print '...setting NN and NNN distance...'
#xyz = np.dot(xyzf, lat)
for nada in range(3):
  for k in [-2,-1,0,1,2]:
    for l in [-2,-1,0,1,2]:
      for i in range(norb):
        for j in range(norb):
          fdif = xyzf[i] - (xyzf[j] +[k, 0] + [0, l])
          dif = np.dot(fdif, lat)
          dist = np.sqrt(np.dot(dif, dif))
          if dist != 0 and dist > 0.000001 and dist < dnn+dd :
            dnn=dist
          elif dist != 0 and dist > 0.000001 and dist < dnnn+dd :
            dnnn=dist
print 'NN distance:', dnn
print 'NNN distance', dnnn
alpha=falpha/dnn
alps=alpha

nor=1./abs(t*np.exp(-alpha*dnn))



print '.....setting hopping.....'
for k in [-2,-1,0,1,2]:
  for l in [-2,-1,0,1,2]:
    if k==0 and l==0 :
      for i in range(norb):
        for j in range(norb):
          if i != j:
            fdif = xyzf[i] - xyzf[j]
            dif = np.dot(fdif, lat)
            dist = np.sqrt(np.dot(dif, dif)) 
            kg.set_hop(nor*t*np.exp(-alpha*dist), i, j, [ k, l], mode="add", allow_conjugate_pair=True)
    else:
      for i in range(norb):
        for j in range(norb):
          fdif = xyzf[i] - (xyzf[j] + [k, 0] + [0, l])
          dif = np.dot(fdif, lat)
          dist = np.sqrt(np.dot(dif, dif))
          kg.set_hop(nor*t*np.exp(-alpha*dist), i, j, [ k, l], mode="add", allow_conjugate_pair=True)


############################
# print tight-binding model
kg.display()

# call function k_path to construct the actual path
(k_vec,k_dist,k_node)=kg.k_path(path,nk)

################### BAND
print '---------------------------------------'
print 'starting calculation'
print '---------------------------------------'
print 'Calculating bands...'

# obtain eigenvalues to be plotted
evals=kg.solve_all(k_vec)

# figure for bandstructure
fig, ax = plt.subplots()
# set range of horizontal axis
ax.set_xlim([0,k_node[-1]])
# put tickmarks and labels at node positions
ax.set_xticks(k_node)
ax.set_xticklabels(label)
# add vertical lines at node positions
for n in range(len(k_node)):
  ax.axvline(x=k_node[n],linewidth=2, color='k')
# put title
ax.set_ylabel("Energy (t)", fontsize=25)
for tick in ax.yaxis.get_major_ticks():
  tick.label1.set_fontsize(25)
for tick in ax.xaxis.get_major_ticks():
  tick.label1.set_fontsize(25)
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params(which='both', width=2, direction='in')
plt.tick_params(which='major', length=7,left=True, right=True)
plt.tick_params(which='minor', length=4,left=True, right=True)

for n in range(norb):
   ax.plot(k_dist,evals[n], 'r', linewidth=1)


# make an PDF figure of a plot
fig.tight_layout()
fig.savefig("band.ps")
###############END BAND

print 'Done.\n'
