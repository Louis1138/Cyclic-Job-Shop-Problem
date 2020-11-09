# -*- coding: utf-8 -*-


import time
from docplex.mp.model import Model
import matplotlib.pyplot as plt

# Que faut-il afficher ?
print_constraints=False
print_cycle_results=True
print_ti=False
print_Kij=False
print_gantt_chart=False

TimeLimit=3600  # secondes
eps=0.5
coeff=0.0001


################### Lecture fichier ###################

#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_3_15_6_#1.dat", "r")
#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_50_5_#1.dat", "r")
#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_50_10_#1.dat", "r")
#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_50_10_#10.dat", "r")
fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_100_10_#9.dat", "r")


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
    data=ligne.rstrip('\n\r').split(" ")
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

print('LBm = ', max(LBm),'; LBj =',max(LBj)*0.5)


#################### MIP model ###################################


# Create a new model
m=Model("CJSP")
m.set_time_limit(TimeLimit)

# Create variables
tau=m.continuous_var(name='tau')
u=[m.continuous_var(name='u'+str(i)) for i in range(nbTaches)]

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
            # constraint dualized
            #m.add_constraint_(K[c]+K[c-1]==1)
            c+=1

penalties=[m.continuous_var(-m.infinity,m.infinity,name='p'+str(i)) for i in range(c//2)]

# Dualized constraints: sum j: K[c]+K[c-1]==1
c=0;
for p in penalties:
    m.add_constraint_(p==K[c]+K[c+1]-1)
    c=c+2
l=[0]*len(penalties)
s=[0]*len(penalties)

######################### Run program ############################
sum_p=[]
start_time=time.time()
for k in range(2):
    m.maximize(tau
               +sum(
                    # Penalties for dualized constraints
                    l_j*p_j for l_j,p_j in zip(l, penalties)))
    solution=m.solve()
    if solution.get_value('tau')==0:
        alpha=m.infinity
    else: alpha=1/solution.get_value('tau')
    obj=solution.get_objective_value()
    if obj==0:
        inv_obj=m.infinity
    else: inv_obj=1/obj
    #print('penalties =',[solution.get_value('p'+str(i)) for i in range(len(penalties))])
    sum_p.append(sum([abs(solution.get_value('p'+str(i))) for i in range(len(penalties))]))
    print(sum_p[-1],"     ",1/solution.get_value('tau'),"         ",time.time()-start_time)


# Test for complementary slackness
    stop=True
    for i in range(len(penalties)):
        if abs(solution.get_value('p'+str(i)))>eps:
            stop=False
            break
        
    if stop:
        print('primal feasible & optimal')
        opt=1
        break
    
    else:
        LB=1.0/(max((max(LBj)*0.5),max(LBm)))
        for i in range(len(penalties)):
            if solution.get_value('p'+str(i))!=0:
                s[i]=coeff/solution.get_value('p'+str(i))**2
            l[i]+=-s[i]*(solution.get_value('p'+str(i)))

if not stop:
    print('feasible solution not found')
    print('lower bound = ',1/solution.get_value('tau'))


run_time=time.time()-start_time
print('runtime =',run_time)
plt.clf()
plt.plot(sum_p)
plt.show()
fichier.close()
