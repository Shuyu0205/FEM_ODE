#Shuyu's work
#This is a program for using finite elements method to solve 1D helmhotz equation with Sommerfeld radiation neumann bounadary!
from numpy import *
from sympy.abc import s
from sympy import *
import math
import cmath
from time import *
import matplotlib.pyplot as plt
###########################################################################
#Specify plotting styles
line_type=['o:','o--','.-']
line_colour_r=['#228B22','#4169E1','#FF1493']#Forestgreen+Royalblue+Deeppink!
line_colour_i=['lime','blue','yellow']
elements_solve=[3,15,120]#Maybe overfit?
Domain=[0,2*pi]#Tell the domain to work on, but not functional yet!
Domain_end=Domain[1]
Coefficient_type='real'#Try the complex coefficient just for fun
###########################################################################
#Build local coordinates
phi_local=[0,(0.5*(1-s)),(0.5*(1+s))]#Local functions
dphi_local=[0,-0.5,0.5]#Store the derivative by hand, also can be done by sy.diff
W_i=[0,5/9,8/9,5/9]#Gauss rule's weight, '0' for indices consistent
S_i=[0,-math.sqrt(3/5),0,math.sqrt(3/5)]#Gauss rule's point
gausspoint=[1,2,3]#Gauss rule's point easy reader for later use
###########################################################################
#Set up
#The thought of this part is to use array to store all the elements and their nodes
#So in the algorithm, we can easily recall them by their indices and names
for N in elements_solve:#Do with different numbers of elements
    begin_time=time()
    Element=[]#To generate N Elements and store them as 1,2,3,...,N for later use
    for n in range(N):
            Element.append(n+1)           
    Globalnodes=[]#To store Global nodes of the whole domian
    for n in range(N+1):
            Globalnodes.append(n+1)           
    Globalsplit=[]#To split the whole domain into elements, with their global nodes
    for n in range(N):
            Globalsplit.append([n+1,n+2])           
    Coordinate=[0]#To split the whole domain and store their actual value on each node
    #All the beginning '0's in lists are just used for making the indices consistent
    #Some very straight forward things don't include this tricky '0'
    #Anyway, go away and check the data example in the beginning if you are confused with my data structure!
    for n in range(N):
            Coordinate.append([(n*Domain_end)/N,((n+1)*Domain_end)/N])           
    Localnodes=[0]#Give each element local nodes
    for n in range(N):
            Localnodes.append([1,2])            
    Globalequation=[]#To give each nodes their global equation number
    for n in Globalnodes:
            if (n==1):
                    Globalequation.append(-1)
            else:
                    Globalequation.append(n-1)                   
    def globalequ(jglobalnode):#Function to look up global equation number by global node number
        return Globalequation[jglobalnode-1]
    def globalnode(jlocal,e):#Function to look up global node number by element number and it's local node
        return Globalsplit[e-1][jlocal-1]
    Localequation=[]#Here gives each node a local equation number
    dof=[]#To store each node's degree of freedom
    localdof=[0]#To store each local node's local dof
    for e in Element:
            count=1
            lst=[]
            jdof=0
            lst2=[]
            for jlocal in Localnodes[e]:
                    jglobal=globalnode(jlocal,e)
                    E_j_global=globalequ(jglobal)
                    if E_j_global==-1:
                            lst.append(-1)
                            lst2.append(0)
                    else:
                            jdof=jdof+1
                            lst.append(count)
                            count=count+1
                            lst2.append(jdof)
            dof.append(jdof)
            Localequation.append(lst)
            localdof.append(lst2)           
    def E_hat_j_e(jdof,e): #A function to look up the global equation number by dof and element number
                    index=localdof[e].index(jdof)
                    j_global=Globalsplit[e-1][index]
                    return globalequ(j_global)           
    def Ndof(e):#A function to read Ndof from one element.
        return dof[e-1]
    def localequ(jlocal,e):#A function to look up the local equation number by local nodes and element number
        return Localequation[e-1][jlocal-1]
    print('It is working on '+str(N)+' 2-node elements')
###########################################################################
#Boundary specify
    U_0=[0,1] #The first zero makes indices consistent, the second one is left Dirichlet boundary
    for i in range(N-1):
            U_0.append(0)#Starting guess
    U_0.append(0)#Should exist here but doing!
###########################################################################
#Main algorithm function
#prepare the global J and r matrix with all 0
#But one thing tricky is make one more row and column
#So that the indices will be consistent with everything starts from 1, not 0
#The extra column and row can easily be deleted after whole calculation
    def solve(N):
            J = []
            lst = []
            r = []
            for i in range(N+1):
                for k in range(N+1):
                    lst.append(0)
                J.append(lst)
                lst = []
                r.append(0)
            for e in Element:#Start the whole algorithm
                    a = Coordinate[e][0]
                    b = Coordinate[e][1]
                    ndof=Ndof(e)
                    r_elements=[]
                    J_elements=[]
                    lst=[]
                    for i in range(ndof+1):
                        for k in range(ndof+1):
                            lst.append(0)
                        J_elements.append(lst)
                        lst = []
                        r_elements.append(0)
                    for i_int in gausspoint:
                        #Strat with setting everything needed.
                            s_int=S_i[i_int]
                            w_int=W_i[i_int]
                            x_s=((a*phi_local[1]+b*phi_local[2]))#The transformation relation
                            Jacob=((a*dphi_local[1]+b*dphi_local[2]))#For integral in local coordinates
                            #f=(30*(sin(math.sqrt(30)*x_s)))  #take sign when it's on right hand
                            f=0
                            k_sqr=16# take sign when it's on left hand. Of couse!
                            k=math.sqrt(k_sqr)
                            u=(U_0[globalnode(1,e)]*phi_local[1])
                            du=(U_0[globalnode(1,e)]*dphi_local[1])#Just because I know the diff of linear function
                            Neumann=k*1j*u
                            Modifier=Neumann/3
                            for j_local in Localnodes[e]:
                                    j_global=globalnode(j_local,e)
                                    E_j_global=globalequ(j_global)
                                    if E_j_global!=-1:
                                            idof=localequ(j_local,e)
                                            #The residual vector should be given here. In galerking form
                                            r_elements[idof]=r_elements[idof]+(((du*dphi_local[j_local]*(Jacob**(-2)))-((k_sqr)*phi_local[j_local]*u)+(f*phi_local[j_local]))*Jacob*w_int).subs(s,s_int)
                                            if (e==N) and (j_local==2):
                                            	r_elements[j_local]=r_elements[j_local]-Modifier
                                            	J_elements[idof][jdof]=J_elements[idof][jdof]-(k*1j)/3
                                            for k_local in Localnodes[e]:
                                                    k_global=globalnode(k_local,e)
                                                    E_k_global=globalequ(k_global)
                                                    if E_k_global!=-1:
                                                       jdof=localequ(k_local,e)
                                                       #The coefficient matrix should be given here.
                                                       J_elements[idof][jdof] = J_elements[idof][jdof]+(((dphi_local[k_local]*dphi_local[j_local]*(Jacob**(-2)))-k_sqr*(phi_local[j_local]*phi_local[k_local]))*Jacob*w_int).subs(s,s_int)
                                                       #if (e==N) and (j_local==2):
                                                       	   
                            for idof in range(ndof):#J and r assembly
                                idof=idof+1
                                E_i_global=E_hat_j_e(idof,e)
                                r[E_i_global]=r[E_i_global]+r_elements[idof]
                                for jdof in range(ndof):
                                    jdof=jdof+1
                                    E_j_global=E_hat_j_e(jdof,e)
                                    J[E_i_global][E_j_global]=J[E_i_global][E_j_global]+J_elements[idof][jdof]
            #print(J)
            J=list(map(lambda x:x[1:],J))#delete first trick column
            J=J[1:]#delete first trick row
            r=r[1:]#delete first row
            print(r)
            r=matrix(r,'complex').T
            J=matrix(J,'complex')        
###########################################################################             
#Solve JU=-r， May apply numerical method here.
            #print(J)
            J_i=J.I
            U=(J_i*(-r)).T
            U=U.tolist()
            U=U[0]#need to apply boundary later, so need to change it to lists to make life easier.
            return U
###########################################################################
#Apply boundary
    U=solve(N)
    U.insert(0,U_0[1])#We store our boundary in the beginning
    U_r=[] #save the real part
    U_i=[] #save the imag part
    #Can do this because of the linearity of a linear system
    #Ax=b <=> ARe(x)+AIm(x)=b
    for u in U:
        U_r.append(u.real)
        U_i.append(u.imag)
    end_time=time()
    run_time=end_time-begin_time
    print('Solving with '+str(N)+' elements takes: ',run_time,'s')
###########################################################################
#Plot by points
    x_set=[]
    y_set=[]
    y_set_r=[]
    y_set_i=[]
    y_set_neg=[]
    plt.plot(x_set,y_set_r,line_type[elements_solve.index(N)],color=line_colour_r[elements_solve.index(N)],label='Re '+'N= '+str(N))#Plot out the legend first
    plt.plot(x_set,y_set_i,line_type[elements_solve.index(N)],color=line_colour_i[elements_solve.index(N)],label='Im '+'N= '+str(N))#Plot out the legend first
    U_r.insert(0,'Allons y!')
    U_i.insert(0,'Bugs') #make index consistent,so drop anything funny on the 0-th member, nobody will find out!
    for e in Element:
        a=Coordinate[e][0]
        b=Coordinate[e][1]
        x=a*phi_local[1]+b*phi_local[2]
        y_r=0
        y_i=0
        for j_local in Localnodes[e]:
            j_global=globalnode(j_local,e)
            U_node_r=U_r[j_global]
            U_node_i=U_i[j_global]
            y_r=y_r+U_node_r*phi_local[j_local]
            y_i=y_i+U_node_i*phi_local[j_local]
        s_split=[-1,0,1]  #how many point to be calculated in one local element, doesn't affect accuracy of the method, can be changed
        for s_p in s_split:
            x_set.append(x.subs(s,s_p))#store points to make plotting easier
            y_set_r.append(y_r.subs(s,s_p))
            y_set_i.append(y_i.subs(s,s_p))
            #y_set_neg.append((-y).subs(s,s_p)) #For poisson problem
            plt.plot(x_set, y_set_r, line_type[elements_solve.index(N)],color=line_colour_r[elements_solve.index(N)])#plot real
            plt.plot(x_set, y_set_i, line_type[elements_solve.index(N)],color=line_colour_i[elements_solve.index(N)])#plot imag
            #plt.plot(x_set,(y_set_neg), line_type[elements_solve.index(N)], color=line_colour[elements_solve.index(N)])
###########################################################################
#Nicely show the figure
plt.title("A FEM solution to "+Coefficient_type+" coefficient helmholtz equ.")
plt.text(x_set[-8],y_set_r[-8]+0.05,"Re",size=15,color='#FF0000')
plt.text(x_set[-8],y_set_i[-8]+0.05,"Im",size=15,color='#FF0000')
plt.xlabel("x")
plt.ylabel("u(x)")
plt.legend(loc=2, bbox_to_anchor=(1.0,1.0),borderaxespad = 0)
plt.show()




