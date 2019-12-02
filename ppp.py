import numpy as np
import matplotlib.pyplot as plt
import math as m


# material properties
poission_ratio = 0.33
youngs_modulus = 210


# geometrical parameters
length = 70
radius = 50

f_ext = 60

# shape function for isoparametric element
E1 = 1
E2 = 1
N = ([[(1-E1)/4,(1-E2)/4],
          [(1+E1)/4,(1-E2)/4],
          [(1+E1)/4,(1+E2)/4],
          [(1-E1)/4,(1+E2)/4]])



# defining elements
nelm = 4
isoparametric_edge = 4
d_o_f = 2
no_nodes = int((np.sqrt(nelm)+1)*(np.sqrt(nelm)+1))

# graphical representation of elements
nelm_length = np.linspace(0,length,np.sqrt(nelm)+1)
nelm_radius = np.linspace(0,radius,np.sqrt(nelm)+1)
yy,xx = np.meshgrid(nelm_length,nelm_radius)

# fig, ax = plt.subplots()
# ax.scatter (xx,yy)
# plt.show()

# deriving all elements coordinates for the calculating jacobian
all_ele_coord = np.zeros((nelm,4,2))
loop = 0
for j in range(int(np.sqrt(nelm))):
    for i in range(int(np.sqrt(nelm))):

        ele_coord = np.matrix(([[nelm_radius[i],nelm_length[j]],
                      [nelm_radius[i+1],nelm_length[j]],
                      [nelm_radius[i+1],nelm_length[j+1]],
                      [nelm_radius[i],nelm_length[j+1]]]))
        ele_coord = ele_coord.reshape((4,2))
        all_ele_coord_zeros = all_ele_coord[loop]+ele_coord  
        all_ele_coord[loop] = all_ele_coord_zeros
        loop += 1


area = all_ele_coord[0,1,0]*all_ele_coord[0,2,1]    

# material rotuine

def material_rotuine(poission_ratio, youngs_modulus):
    c_1 = (youngs_modulus)/((1+poission_ratio)*(1-(2*poission_ratio)))
    c = c_1*np.matrix([[1-poission_ratio,poission_ratio,poission_ratio,0],
                       [poission_ratio,1-poission_ratio,poission_ratio,0],
                       [poission_ratio,poission_ratio,1-poission_ratio,0],
                       [0,0,0,((1-2*poission_ratio)/2)]])
                                                                         
    return c


# element rotuine
def element_rotuine(E1,E2,all_ele_coord,N,area):
  k_all_ele = np.zeros((nelm,isoparametric_edge*2,isoparametric_edge*2))
  for i in range(nelm):
    derivative_N = 1/4*np.matrix([[-(1-E2),(1-E2),(1+E2),-(1+E2)],
                                [-(1-E1),-(1+E1),(1+E1),(1-E1)]])
    x_y_ele = all_ele_coord[i]
    jacobi_1 = derivative_N*x_y_ele
    jacobi_inverse = np.linalg.inv(jacobi_1)
    B_1_matrix = jacobi_inverse*derivative_N
    radius = (np.sum(x_y_ele[:,0]))/4
    
    N_1 = N_2 = N_3 = N_4 = 1/4
    B_matrix = ([[B_1_matrix[0,0],0,B_1_matrix[0,1],0,B_1_matrix[0,2],0,B_1_matrix[0,3],0],
                  [0,B_1_matrix[1,0],0,B_1_matrix[1,1],0,B_1_matrix[1,2],0,B_1_matrix[1,3]],
                  [N_1/radius,0,N_2/radius,0,N_3/radius,0,N_4/radius,0],
                  [B_1_matrix[1,0],B_1_matrix[0,0],B_1_matrix[1,1],B_1_matrix[0,1],B_1_matrix[1,2],B_1_matrix[0,2],B_1_matrix[1,3],B_1_matrix[0,3]]])
    k = 2*3.14*radius*area*np.transpose(B_matrix)*material_rotuine(poission_ratio, youngs_modulus)*B_matrix
    k = k_all_ele[i] + k
    k_all_ele[i] = k

  return k_all_ele

all_k_ele = element_rotuine(E1,E2,all_ele_coord,N,area)



# nodes of the elements
x_cells = int(np.sqrt(nelm))
y_cells = int(np.sqrt(nelm))
elements = np.arange((x_cells+1)*(y_cells+1)).reshape(y_cells+1,x_cells+1)
ID = ([])
for i in range(y_cells):
    for j in range(x_cells):
        id = np.matrix([elements[i,j],elements[i,j+1],elements[i+1,j],elements[i+1,j+1]])
        ID = np.append(ID,id)
 
sum = ([])  
for i in ID:
    multi = np.array([(i*2),((i*2)+1)])
    sum = np.append(sum,multi)

sum = sum.astype(int)
a_column = sum.flatten().reshape(nelm,isoparametric_edge*d_o_f)



#  assignment matrix for all the elements
a_rows = (np.arange(isoparametric_edge*d_o_f)).astype(int)
all_a = np.zeros((nelm,isoparametric_edge*d_o_f , no_nodes*d_o_f))
for k in range(nelm):
  a = np.zeros((isoparametric_edge*d_o_f , no_nodes*d_o_f))
  for i,j in zip(a_rows,a_column[k]):
    a[i,j] = 1
  all_a[k] = a
print(all_a)
# calculating global stiffness matrix
k_global = 0
for i in range(nelm):
  assembly_1 = np.dot(np.transpose(all_a[i]),all_k_ele[i])
  assembly = np.dot(assembly_1,all_a[i])
  k_global = k_global + assembly

print(k_global)


# f = np.zeros((18,1))
# x = (np.arange(no_nodes))
# y = x[-(int(np.sqrt(nelm)+1)):]
# for i in y:
#   a = (i*2)+1
#   f[a] = f_ext

# v = np.linalg.inv(k_global)
# print(v)

 
print()

































