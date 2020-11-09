# -*- coding: utf-8 -*-

#import sys
import time
from docplex.mp.model import Model
from cplex.exceptions import CplexError
#from myfct.py import *
#import plotly.figure_factory as ff
#from plotly.offline import plot

# Que faut-il afficher ?
print_constraints=False
print_cycle_results=True
print_ti=False
print_Kij=False
print_gantt_chart=False

TimeLimit=10000  # secondes

################### Lecture fichier ###################

#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_3_15_6_#1.dat", "r")
#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_50_5_#1.dat", "r")
#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_50_10_#1.dat", "r")
#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_50_10_#10.dat", "r")
fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_100_10_#9.dat", "r")

"""
f=sys.argv[1]
fichier = open(f, "r")
"""

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
if print_constraints:
    print("Constraints :")
for  ligne  in fichier :
    data = ligne.rstrip('\n\r').split(" ")
    data.remove('')
    constraint=[int(x) for x in data]
    if print_constraints:
        print(constraint)
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

#################### MIP model ###################################


try:
#if True:
    # Create a new model
    m=Model("CJSP")
    m.set_time_limit(TimeLimit)
    #m.setParam('OutputFlag', False) # turns off solver chatter
    
    # Create variables
    tau=m.continuous_var(lb=0,ub=1,name='tau')
    
    #m.addConstr(tau <= 1.0/(max(max(LBm),max(LBj)*0.5)), "UB")
    u=[m.continuous_var(name='u'+str(i)) for i in range(nbTaches)]

    # Add constraint: non-reentrance contraints
    m.add_constraint_(1>=max(p)*tau,"non reentrance")

    # Add constraint: conjunctive contraints
    for i in range(nbArcs):
        m.add_constraint_(u[ConstraintsUniform[i][1]]-u[ConstraintsUniform[i][0]]>=ConstraintsUniform[i][2]*tau-ConstraintsUniform[i][3],"precedence"+str(i))

    K=[]
    c=0;
    for l in range(nbMachines):
        for k in range(nbTachesMachine[l]):
            for j in range(k+1,nbTachesMachine[l]) :
                K.append(m.integer_var(lb=-m.infinity,ub=m.infinity,name='K'+str(c)))
                K.append(m.integer_var(lb=-m.infinity,ub=m.infinity,name='K'+str(c+1)))
                m.add_constraint_(u[Machine[l][j]]+K[c]>=u[Machine[l][k]]+tau*p[Machine[l][k]])
                c+=1
                m.add_constraint_(u[Machine[l][k]]+K[c]>=u[Machine[l][j]]+tau*p[Machine[l][j]])
                m.add_constraint_(K[c]+K[c-1]==1)
                c+=1

    # Set objective
    m.maximize(tau)

    # Solve problem
    print("solving...")
    start_time=time.time()
    solution=m.solve()
    run_time=time.time()-start_time
    alpha=1/solution.get_value('tau')
    ## status test 2:opt, 9:TimeLimit
    if str(m.get_solve_status())=="JobSolveStatus.OPTIMAL_SOLUTION":
        opt="optimal solution"
    elif str(m.get_solve_status())=="JobSolveStatus.FEASIBLE_SOLUTION":
        opt="feasible solution (TimeLimit reached)"
    else :
        opt="!!!!!!!!!!!!!!!!!!!new_status!!!!!!!!!!!!!!!!!!"
    if print_cycle_results:
        print('\n------ cycle time ------')
        print('alpha =',alpha)
        print('LBm = ', max(LBm),'; LBj =',max(LBj)*0.5)
        print(opt)
        print('runtime =',run_time)
        print('-------------------------------\n')

    # Display variables
    if print_ti:
        for i in range(nbTaches):
            print("t[",i,"]", alpha*solution.get_value('u'+str(i)))
    if print_Kij:
        c=0
        for l in range(nbMachines):
            for k in range(nbTachesMachine[l]):
                for j in range(k+1,nbTachesMachine[l]) :
                    print("K",Machine[l][k],"->",Machine[l][j],":",solution.get_value('K'+str(c)))
                    c+=1
                    print("K",Machine[l][j],"->",Machine[l][k],solution.get_value('K'+str(c)))
                    c+=1

 
#################### storing results ###################################

#    fichier = open("CJSP_resume_Python.csv", "a")
#    fichier.write(sys.argv[1] + ";" + str(opt) +";" + str(alpha) +  ";" + str(max(LBj)*0.5) +  ";" + str(max(LBm)) +  ";" + str(m.Runtime) + ";" + str(m.NodeCount) + ";" + str(m.SolCount))
#    fichier.write("\n")
#    fichier.close()

except CplexError:
    print('Encountered a Cplex error')

except AttributeError:
    print('Encountered an attribute error')


######################################################################

fichier.close()

if print_ti:
    for i in range(nbTaches):
        print("t[",i,"]",alpha*solution.get_value('u'+str(i)))

