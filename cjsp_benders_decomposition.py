# -*- coding: utf-8 -*-

from docplex.mp.model import Model
import sys
import time


epsilon=0.01
TimeLimit=20.0  # secondes

print_constraints=False
print_cycle_results=True
print_ti=True
print_Kij=False
print_gantt_chart=False

################### Lecture fichier ###################

#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_2_4_0_#1.dat", "r")
fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_3_15_6_#1.dat", "r")
#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_50_5_#1.dat", "r")
#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_50_10_#1.dat", "r")
#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_50_10_#10.dat", "r")
#fichier=open("C:\\Users\\Louis\\Documents\\2A\\Stage\\Python\\cjsp_5_100_10_#9.dat", "r")

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



#################### MIP model ###################################

#try:
if True:
    mdl=Model('CJSP')
    mdl.clear()
    mdl.set_time_limit(TimeLimit)
    tau=mdl.continuous_var(lb=0,ub=1,name='tau')
    u=[mdl.continuous_var(name='u'+str(i)) for i in range(nbTaches)]
    mdl.add_constraint_(1>=max(p)*tau,"non reentrance")
    for i in range(nbArcs):
        mdl.add_constraint_(u[ConstraintsUniform[i][1]]-u[ConstraintsUniform[i][0]]>=ConstraintsUniform[i][2]*tau-ConstraintsUniform[i][3],"precedence"+str(i))
    K=[]
    c=0;
    for l in range(nbMachines):
        for k in range(nbTachesMachine[l]):
            for j in range(k+1,nbTachesMachine[l]) :
                K.append(mdl.integer_var(-mdl.infinity,mdl.infinity,name='K'+str(c)))
                K.append(mdl.integer_var(-mdl.infinity,mdl.infinity,name='K'+str(c+1)))
                mdl.add_constraint_(u[Machine[l][j]]+K[c]>=u[Machine[l][k]]+tau*p[Machine[l][k]])
                c+=1
                mdl.add_constraint_(u[Machine[l][k]]+K[c]>=u[Machine[l][j]]+tau*p[Machine[l][j]])
                mdl.add_constraint_(K[c]+K[c-1]==1)
                c+=1
    mdl.maximize(tau)

    # Solve problem
    print("solving...")
    start_time=time.time()
    solution_mdl=mdl.solve()
    run_time=time.time()-start_time
    alpha=1/solution_mdl.get_value('tau')
    ## status test 2:opt, 9:TimeLimit
    if str(mdl.get_solve_status())=="JobSolveStatus.OPTIMAL_SOLUTION":
        opt="optimal solution"
    elif str(mdl.get_solve_status())=="JobSolveStatus.FEASIBLE_SOLUTION":
        opt="feasible solution (TimeLimit reached)"
    else :
        opt="!!!!!!!!!!!!!!!!!!!new_status!!!!!!!!!!!!!!!!!!"
    if print_cycle_results:
        print('\n_________CLASSIC METHOD________')
        print('Solved in',run_time,'seconds')
        print('alpha =',alpha)

    # Display variables
    if print_ti:
        for i in range(nbTaches):
            print("t[",i,"]", alpha*solution_mdl.get_value('u'+str(i)))
    if print_Kij:
        c=0
        for l in range(nbMachines):
            for k in range(nbTachesMachine[l]):
                for j in range(k+1,nbTachesMachine[l]):
                    print("K",Machine[l][k],"->",Machine[l][j],":",solution_mdl.get_value('K'+str(c)))
                    c+=1
                    print("K",Machine[l][j],"->",Machine[l][k],solution_mdl.get_value('K'+str(c)))
                    c+=1
    mdl.clear()

#################### storing results ###################################

#    fichier = open("CJSP_resume_Python.csv", "a")
#    fichier.write(sys.argv[1] + ";" + str(opt) +";" + str(alpha) +  ";" + str(max(LBj)*0.5) +  ";" + str(max(LBm)) +  ";" + str(m.Runtime) + ";" + str(m.NodeCount) + ";" + str(m.SolCount))
#    fichier.write("\n")
#    fichier.close()

#except CplexError:
#    print('Encountered a Cplex error')

#except AttributeError:
#    print('Encountered an attribute error')


######################################################################

# A=Matrix of constraints of the sub pb
A=[]
for i in range(nbArcs):
    L=[0 for i in range(nbTaches)]+[ConstraintsUniform[i][2]]
    L[ConstraintsUniform[i][1]],L[ConstraintsUniform[i][0]]=-1,1
    A.append(L)

for l in range(nbMachines):
    for k in range(nbTachesMachine[l]):
        for j in range(k+1,nbTachesMachine[l]):
            L=[0 for i in range(nbTaches)]+[p[Machine[l][k]]]
            L[Machine[l][j]],L[Machine[l][k]]=-1,1
            A.append(L)
            L=[0 for i in range(nbTaches)]+[p[Machine[l][j]]]
            L[Machine[l][k]],L[Machine[l][j]]=-1,1
            A.append(L)
# Transpose of the matrix of constraints of the sub pb (matrix of contraints of the dual sub pb)
# used in function dual_SP in order to set the constraints conviniently
B=[[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]

# Subproblem Solver Function
def dual_SP(K_):
    # K_ = list of [Kij,i,j]
    # Solve dual sub pb    
    dsp=Model('dual_SP')
    dsp.clear()
    dsp.set_time_limit(TimeLimit)
    pi=[dsp.continuous_var(name='pi'+str(i)) for i in range(len(A))]
    for i in range(len(B)-1):
        dsp.add_constraint_(sum(B[i][j]*pi[j] for j in range(len(pi)))>=0,ctname='ct'+str(i))
    dsp.add_constraint_(sum(B[-1][j]*pi[j] for j in range(len(pi)))>=1,ctname='ct'+str(len(B)-1))
    dsp.minimize(sum(ConstraintsUniform[i][3]*pi[i] for i in range(nbArcs))
              +sum(K_[i][0]*pi[i+nbArcs] for i in range(len(K_))))
    solution_dsp=dsp.solve()

    # Case 1: unbounded dual sub pb -> infeasible sub pb
    if str(dsp.get_solve_status()) in('JobSolveStatus.UNBOUNDED_SOLUTION','JobSolveStatus.INFEASIBLE_OR_UNBOUNDED_SOLUTION'):
        #print("unbounded dual sub pb -> infeasible sub pb")
        dsp.add_constraint_(sum(ConstraintsUniform[i][3]*pi[i] for i in range(nbArcs))
              +sum(K_[i][0]*pi[i+nbArcs] for i in range(len(K_)))>=-1e10)
        #dsp.add_constraints_(pi[i]<=9e19 for i in range(len(pi)))
        #dsp.add_constraints_(pi[i]>=-9e19 for i in range(len(pi)))
        solution_dsp=dsp.solve()
        if str(dsp.get_solve_status()) in('JobSolveStatus.UNBOUNDED_SOLUTION','JobSolveStatus.INFEASIBLE_OR_UNBOUNDED_SOLUTION'):
            print('pb!!!!!!!!!!!!!!!!!!!!!! borne artificielle défaillante dans dsp')
            sys.exit()
        #extreme_ray=[]
        #for i in range(len(pi)):
        #    if solution_dsp.get_value('pi'+str(i))>=8e19: extreme_ray.append(1)
        #    elif solution_dsp.get_value('pi'+str(i))<=-8e19: extreme_ray.append(-1)
        #    else: extreme_ray.append(0)
        extreme_ray=[solution_dsp.get_value('pi'+str(i)) for i in range(len(pi))]
        #print('extreme_ray =',[i for i in extreme_ray])
        return('infeasible',extreme_ray)

    # Case 2: feasible dual sub pb -> feasible sub pb
    elif str(dsp.get_solve_status())=="JobSolveStatus.OPTIMAL_SOLUTION":
        #print('feasible dual sub pb -> feasible sub pb : alpha =',solution_dsp.get_objective_value())
        dual_vals=[dsp.dual_values(dsp.find_matching_linear_constraints('ct'+str(i)))[0] for i in range(len(B))]
        return('feasible',[solution_dsp.get_value('pi'+str(i)) for i in range(len(pi))],dual_vals,solution_dsp.get_objective_value())

    # Case 3: Other... Error
    else:
        print('PROBLEM !!!!!!!!!!!!!!!!! dual_SP case 3')
        print('case not taken into account yet')
        sys.exit()

# Solve the problem with the benders decomposition
# Define master pb
mp=Model('Master_Pb')
mp.clear()
mp.set_time_limit(TimeLimit)
# K_ and K : lists of [Kij,i,j]
K_,K=[],[]
for l in Machine:
    for i in l[:-1]:
        for j in l[l.index(i)+1:]:
            K.append([mp.integer_var(lb=-mp.infinity,name='K'+str(i)+'.'+str(j)),i,j])
            K.append([mp.integer_var(lb=-mp.infinity,name='K'+str(j)+'.'+str(i)),j,i])
            K_.append([1,i,j])
            K_.append([0,j,i])
            mp.add_constraint_(K[-1][0]+K[-2][0]==1)

UB=mp.infinity
LB=-mp.infinity
z=mp.continuous_var(lb=-mp.infinity,name='z')
mp.maximize(z)
nb_iterations=0
start_time=time.time()
print()
solution_dual_SP=[1,1,1,1]
while UB-LB>=epsilon or solution_dual_SP[0]!='feasible':
    nb_iterations+=1
    solution_dual_SP=dual_SP(K_)
    if solution_dual_SP[0]=='infeasible':
        pi_=solution_dual_SP[1]
        #print('add feasibility cut')
        mp.add_constraint_(sum(ConstraintsUniform[i][3]*pi_[i] for i in range(nbArcs))
              +sum(K[i][0]*pi_[i+nbArcs] for i in range(len(K)))>=0)
    elif solution_dual_SP[0]=='feasible':
        LB=max(LB,solution_dual_SP[3])
        pi_=solution_dual_SP[1]
        #print('add optimality cut')
        mp.add_constraint_(sum(ConstraintsUniform[i][3]*pi_[i] for i in range(nbArcs))
              +sum(K[i][0]*pi_[i+nbArcs] for i in range(len(K)))>=z)
    solution_mp=mp.solve()
    if str(mp.get_solve_status()) in('JobSolveStatus.UNBOUNDED_SOLUTION','JobSolveStatus.INFEASIBLE_OR_UNBOUNDED_SOLUTION'):
        mp.add_constraint(z<=1000)
        solution_mp=mp.solve()
    if str(mp.get_solve_status()) in('JobSolveStatus.UNBOUNDED_SOLUTION','JobSolveStatus.INFEASIBLE_OR_UNBOUNDED_SOLUTION'):
        print('pb!!!!!!!!!!!!!!!!!!!!!! borne artificielle défaillante dans le master')
    """
    for l in Machine:
        for i in l[:-1]:
            for j in l[l.index(i)+1:]:
                print('K',i,j,solution_mp.get_value('K'+str(i)+'.'+str(j)))
                print('K',j,i,solution_mp.get_value('K'+str(j)+'.'+str(i)))
    """
    UB=solution_mp.get_objective_value()
    if nb_iterations%20==0:
        print('UB-LB=',UB-LB)
    for i in range(len(K_)):
        K_[i][0]=solution_mp.get_value('K'+str(K_[i][1])+'.'+str(K_[i][2]))

run_time=time.time()-start_time

fichier.close()

dual_variables=solution_dual_SP[2]
alpha=1/solution_dual_SP[3]
if print_cycle_results:
    print('\n_________BENDERS METHOD________')
print('Solved in',nb_iterations,'iterations and in',run_time,'seconds')
print('alpha =',alpha)

# Display variables
if print_ti:
    for i in range(nbTaches):
        print("t[",i,"]",alpha*dual_variables[i])
if print_Kij:
    for l in Machine:
        for i in l[:-1]:
            for j in l[l.index(i)+1:]:
                print('K',i,j,solution_mp.get_value('K'+str(i)+'.'+str(j)))
                print('K',j,i,solution_mp.get_value('K'+str(j)+'.'+str(i)))


