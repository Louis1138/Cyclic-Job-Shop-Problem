# -*- coding: utf-8 -*-

#import sys
from docplex.cp.model import CpoModel
import time
#from myfct.py import *
#import plotly.figure_factory as ff
#from plotly.offline import plot


TimeLimit=60  # secondes
   
################### Lecture fichier ###################

fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_3_15_6_#1.dat", "r")
#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_50_5_#1.dat", "r")
#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_50_10_#1.dat", "r")
#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_50_10_#10.dat", "r")
#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_100_10_#9.dat", "r")


ligne=fichier.readline()
data=ligne.rstrip('\n\r').split(" ")
nbTaches=int(data[0])
nbArcs=int(data[1])
print("nbtaches :",nbTaches)
print("nbArcs :",nbArcs)

ligne=fichier.readline()
data=ligne.rstrip('\n\r').split(" ")
data.remove('')
p=[int(x) for x in data]

ligne=fichier.readline()
data=ligne.rstrip('\n\r').split(" ")
data.remove('')
m=[int(x) for x in data]
nbMachines=max(m)
print("nbMachines :",nbMachines)

ConstraintsUniform=[]
for  ligne  in fichier :
    data = ligne.rstrip('\n\r').split(" ")
    data.remove('')
    constraint=[int(x) for x in data]
    ConstraintsUniform.append(constraint)


#################### Data processing ########################

nbJobs=0
a=0
LBj=[]
for  i  in range(nbArcs):
    ConstraintsUniform[i][0]-=1
    ConstraintsUniform[i][1]-=1
    a+=ConstraintsUniform[i][2]
    if ConstraintsUniform[i][1]==nbTaches-1:
        LBj.append(a)
        a=0
        nbJobs+=1

Machine = [[] for _ in range(nbMachines)]
nbTachesMachine=[0]*(nbMachines)

for i in range(nbTaches):
    Machine[m[i]-1].append(i)
    nbTachesMachine[m[i]-1]=nbTachesMachine[m[i]-1]+1;
# Machine[n] contains list of tasks processed on machine n+1

LBm=[]
for i in range(nbMachines):
    LBm.append(sum (p[p_i] for p_i in Machine[i]))


#################### CP model ###################################

# Create a new model
m=CpoModel()

# Create variables
tau=m.integer_var(min=0,max=int(10000/max(LBm)),name='tau')
u=[m.integer_var(name='u'+str(i),min=0,max=int(10000/max(LBm))) for i in range(nbTaches)]

# Add constraint: conjunctive contraints
for i in range(nbArcs):
    m.add(u[ConstraintsUniform[i][1]]-u[ConstraintsUniform[i][0]]>=ConstraintsUniform[i][2]*tau-ConstraintsUniform[i][3]*10000)

K=[]
c=0;
for l in range(nbMachines):
    for k in range(nbTachesMachine[l]):
        for j in range(k+1,nbTachesMachine[l]) :
            K.append(m.integer_var(min=-int(-10000),max=int(10000),name='K'+str(c)))
            K.append(m.integer_var(min=-int(10000),max=int(10000),name='K'+str(c+1)))
            m.add(u[Machine[l][j]]+K[c]*10000>=u[Machine[l][k]]+tau*p[Machine[l][k]])
            c+=1
            m.add(u[Machine[l][k]]+K[c]*10000>=u[Machine[l][j]]+tau*p[Machine[l][j]])
            m.add(K[c]+K[c-1]==1)
            c+=1

# Set objective
m.add(m.maximize(tau))

# Solve problem
print("solving...")
start_time=time.time()
solution=m.solve(TimeLimit=60)



run_time=time.time()-start_time
alpha=10000/float(solution.get_value(tau))

print('\n------ cycle time ------')
print('alpha =',alpha)
print('LBm = ', max(LBm),'; LBj =',max(LBj)*0.5)
"""print(opt)
"""
print('runtime =',run_time)
print('-------------------------------\n')

# Display variables
for i in range(nbTaches):
    print("t[",i,"]", alpha*float(solution.get_value('u'+str(i))))



"""
c=0
for l in range(nbMachines):
    for k in range(nbTachesMachine[l]):
        for j in range(k+1,nbTachesMachine[l]) :
            print("K",Machine[l][k],"->",Machine[l][j],":",solution.get_value('K'+str(c)))
            c+=1
            print("K",Machine[l][j],"->",Machine[l][k],solution.get_value('K'+str(c)))
            c+=1
"""


