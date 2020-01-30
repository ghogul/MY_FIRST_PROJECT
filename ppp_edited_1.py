import numpy as np
import matplotlib.pyplot as plt



'''
=======================================================
                INITIAL PARAMETERS
=======================================================
'''
'''***MATERIAL PARAMETER***'''
poission_ratio = 0.3
youngs_modulus = 200000 #Mpa


'''***GEOMETRICAL PARAMETER***'''
length = 70  #mm
radius = 5   #mm


'''***TIME LIKE PARAMETERS***'''
tau = 0
time_step = 0.1
total_time = 1

'''***EXTERNAL LOADING***'''
f_ext = 400  #Mpa
yield_stress = 350 #Mpa


'''***lame's constant***'''
mu = youngs_modulus/(2*(1+poission_ratio))
lamda = (poission_ratio*youngs_modulus)/((1+poission_ratio)*(1-2*poission_ratio))
k = youngs_modulus/(3*(1-(2*poission_ratio)))


'''***DISCRETIZATION PARAMETERS***'''
nelm = 4
isoparametric_edge = 4
d_o_f = 2
no_nodes = int((np.sqrt(nelm)+1)*(np.sqrt(nelm)+1))
length_element = radius/np.sqrt(nelm)




# graphical representation of elements
nelm_length = np.linspace(0,length,np.sqrt(nelm)+1)
nelm_radius = np.linspace(0,radius,np.sqrt(nelm)+1)
yy,xx = np.meshgrid(nelm_length,nelm_radius)

# fig, ax = plt.subplots()
# ax.scatter (xx,yy)
# plt.show()


'''***shape function for isoparametric element***'''
E1 = 0
E2 = 0
N = ([[(1-E1)/4,(1-E2)/4],
          [(1+E1)/4,(1-E2)/4],
          [(1+E1)/4,(1+E2)/4],
          [(1-E1)/4,(1+E2)/4]])
weight = 2
'''***deriving all elements coordinates for the calculating jacobian***'''
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


# def external_force(f_ext,length_element,nelm,all_ele_coord,i,all_external_force):
#     if i < 2:
#         all_external_force[i] = (np.zeros((isoparametric_edge*2,1)))
#         return (np.zeros((isoparametric_edge*2,1)))

#     else:
#         x_y_ele = all_ele_coord[i]
#         radius = x_y_ele[1,0]
#         external_force = (2*3.14*radius*(length_element/2)*f_ext*(1/4))*np.transpose(np.matrix([0,0,0,0,0,4,0,4]))
#         all_external_force[i] = external_force
#         return (external_force)   
    
    
# all_external_force = np.zeros((nelm,8,1))
# for i in range(nelm):
#     external_force(f_ext,length_element,nelm,all_ele_coord,i,all_external_force)


# print(all_external_force)





# for i in range(nelm):
#         if i < 2:
#             return (np.zeros((isoparametric_edge*2,1)))

#         else:
#             x_y_ele = all_ele_coord[i]
#             radius = x_y_ele[1,0]
#             external_force = (2*3.14*radius*(length_element/2)*f_ext*(1/4))*np.transpose(np.matrix([0,0,0,0,0,4,0,4]))
#             return (external_force)



# for i in range(nelm):
#     if i < 2:
#         print("if",i)

#     else:
#         print("else",i)
    




# print(length_element)


# *(np.matrix([0,0,0,0,0,4,0,4]).transpose)





# print(2*np.matrix([1,1]))



# nodes of the elements
# x_cells = int(np.sqrt(nelm))
# y_cells = int(np.sqrt(nelm))
# elements = np.arange((x_cells+1)*(y_cells+1)).reshape(y_cells+1,x_cells+1)
# for i in range(0,(x_cells+1),1):
#     if i % 2 != 0:
#         elements_1 = elements[i]
#         elements[i] = elements_1[::-1]
# ID = ([])
# for i in range(y_cells):
#     for j in range(x_cells):
#         id = np.matrix([elements[i,j],elements[i,j+1],elements[i+1,j+1],elements[i+1,j]])
#         ID = np.append(ID,id)
 
# # print(ID)
ID = ([0, 1, 4, 5, 1, 2, 3, 4, 4, 3, 8, 7,5, 4, 7, 6,])
sum = ([])  
for i in ID:
    multi = np.array([(i*2),((i*2)+1)])
    sum = np.append(sum,multi)

sum = sum.astype(int)
a_column = sum.flatten().reshape(nelm,isoparametric_edge*d_o_f)
# print(sum)


#  assignment matrix for all the elements
a_rows = (np.arange(isoparametric_edge*d_o_f)).astype(int)
all_a = np.zeros((nelm,isoparametric_edge*d_o_f , no_nodes*d_o_f))
for k in range(nelm):
  a = np.zeros((isoparametric_edge*d_o_f , no_nodes*d_o_f))
  for i,j in zip(a_rows,a_column[k]):
    a[i,j] = 1
  all_a[k] = a




# print(all_a)



# k = np.array([[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#               [0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#               [0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#               [0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#               [0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0],
#               [0,0,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0],
#               [0,0,0,0,0,0,7,0,0,0,0,0,0,0,0,0,0,0],
#               [0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0],
#               [0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0],
#               [0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0],
#               [0,0,0,0,0,0,0,0,0,0,11,0,0,0,0,0,0,0],
#               [0,0,0,0,0,0,0,0,0,0,0,12,0,0,0,0,0,0],
#               [0,0,0,0,0,0,0,0,0,0,0,0,13,0,0,0,0,0],
#               [0,0,0,0,0,0,0,0,0,0,0,0,0,14,0,0,0,0],
#               [0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,0,0,0],
#               [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0],
#               [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,17,0],
#               [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,18]])
# print(k)

# global_stiffness_matrix = 0
# for i in range(nelm):
#     assembly_1 = np.dot((np.transpose(all_a[i])),k)
#     assembly = np.dot(assembly_1,all_a[i])
#     global_stiffness_matrix = global_stiffness_matrix + assembly

# print(global_stiffness_matrix)


# reduced_1 = k[2:,2:]
# reduced_2 = np.delete(reduced_1,1,1)
# reduced_3 = np.delete(reduced_2,1,0)
# reduced_4 = np.delete(reduced_3,2,1)
# reduced_5 = np.delete(reduced_4,2,0)
# reduced_6 = np.delete(reduced_5,2,1)
# reduced_7 = np.delete(reduced_6,2,0)
# reduced_8 = np.delete(reduced_7,7,1)
# reduced_9 = np.delete(reduced_8,7,0)
# print(reduced_9.shape)




# a = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18])
# b = a.reshape(18,1)
# # print(b)

# c = np.delete(b,0,0)
# d = np.delete(c,0,0)
# e = np.delete(d,1,0)
# f = np.delete(e,2,0)
# g = np.delete(f,2,0)
# h = np.delete(g,7,0)
# print(h.shape)




'''
=========================================================================================================================================================
                                                          material rotuine
=========================================================================================================================================================
'''

def material_rotuine(poission_ratio, youngs_modulus,initial_displacement,global_plastic_strain,yield_stress,mu,lamda):
    # epsilon = np.dot(B_matrix , initial_displacement)
    # strain = epsilon
    
    c_1 = (youngs_modulus)/((1+poission_ratio)*(1-(2*poission_ratio)))
    C = c_1*np.matrix([[1-poission_ratio,poission_ratio,poission_ratio,0],
                       [poission_ratio,1-poission_ratio,poission_ratio,0],
                       [poission_ratio,poission_ratio,1-poission_ratio,0],
                       [0,0,0,((1-2*poission_ratio)/2)]])
    # trial_stress = np.dot(C , (epsilon - global_plastic_strain))
    # trial_stress_deviatoric = trial_stress - (1/3) * np.sum(trial_stress)
    # trial_stress_equivalent = np.sqrt((3/2) * (np.sum(trial_stress_deviatoric)**2))
    # print(strain)
    # if trial_stress_equivalent-yield_stress < 0:
        # print("elastic")
    # stress = trial_stress
    return C

    # else:
    #     # print("plastic")
    #     delta_lamda = (trial_stress_equivalent-yield_stress)/(3*mu)
    #     current_stress = ((1/3)*np.sum(trial_stress)*np.ones(4).reshape(4,1)) + (((trial_stress_equivalent-(3*mu*delta_lamda)))/trial_stress_equivalent)*trial_stress_deviatoric
    #     current_stress_deviatoric = current_stress - (1/3)*np.sum(current_stress)
    #     current_stress_equivalent = np.sqrt ((3/2)*(np.sum(current_stress_deviatoric)**2))
    #     plastic_strain = global_plastic_strain + ((delta_lamda * 3/2) / current_stress_equivalent) * current_stress_deviatoric
    #     stress = current_stress
    #     c_t_1st = (1/3) * (3*lamda + 2*mu)*np.ones((4,4))
    #     identity_deviatoric = ((1/2)*(np.eye(4)+np.eye(4)))-((1/3)*np.ones((4,4)))
    #     c_t_2nd = (2*mu)*((trial_stress_equivalent-(3*mu*delta_lamda))/trial_stress_equivalent)*identity_deviatoric
    #     c_t_3rd = (3*mu / (trial_stress_equivalent)**2)*(trial_stress_deviatoric*trial_stress_deviatoric.reshape(1,4))
    #     C_t = c_t_1st + c_t_2nd - c_t_3rd
    #     return C_t,stress

'''
=========================================================================================================================================================
                                                      element rotuine
=========================================================================================================================================================
'''
# element rotuine
def element_rotuine(E1,E2,all_ele_coord,N,area,weight):
  k_all_ele = np.zeros((nelm,isoparametric_edge*2,isoparametric_edge*2))
  internal_force_matrix_all_ele = np.zeros((nelm,8,1))
  for i in range(nelm):
    derivative_N = 1/4*np.matrix([[-(1-E2),(1-E2),-(1+E2),(1+E2)],
                                [-(1-E1),-(1+E1),(1-E1),(1+E1)]])
    x_y_ele = all_ele_coord[i]
    jacobi_1 = derivative_N*x_y_ele
    jacobi_inverse = np.linalg.inv(jacobi_1)
    B_1_matrix = jacobi_inverse*derivative_N
    radius = x_y_ele[1,0]
    
    N_1 = N_2 = N_3 = N_4 = 1/4
    B_matrix = ([[B_1_matrix[0,0],0,B_1_matrix[0,1],0,B_1_matrix[0,2],0,B_1_matrix[0,3],0],
                  [0,B_1_matrix[1,0],0,B_1_matrix[1,1],0,B_1_matrix[1,2],0,B_1_matrix[1,3]],
                  [N_1/radius,0,N_2/radius,0,N_3/radius,0,N_4/radius,0],
                  [B_1_matrix[1,0],B_1_matrix[0,0],B_1_matrix[1,1],B_1_matrix[0,1],B_1_matrix[1,2],B_1_matrix[0,2],B_1_matrix[1,3],B_1_matrix[0,3]]])
    c_matrix = material_rotuine(poission_ratio, youngs_modulus,initial_displacement,global_plastic_strain,yield_stress,mu,lamda)
    K = weight*2*3.14*radius*area*(np.transpose(B_matrix))@(c_matrix)@(B_matrix)*np.linalg.det(jacobi_1)
    # internal_force_matrix = weight*2*3.14*radius*area*np.transpose(B_matrix)*stress_matrix*np.linalg.det(jacobi_1)
    # internal_force_matrix_all_ele[i] = internal_force_matrix
    k_all_ele[i] = K
    
    return B_matrix


initial_displacement = np.random.rand((isoparametric_edge*2),1)
global_plastic_strain = np.zeros((isoparametric_edge,1))

# k_all = element_rotuine(E1,E2,all_ele_coord,N,area,weight)

# c = material_rotuine(poission_ratio, youngs_modulus,initial_displacement,global_plastic_strain,yield_stress,mu,lamda)

# global_stiffness_matrix = 0
# for i in range(nelm):
#     assembly_1 = np.dot((np.transpose(all_a[i])),c)
#     assembly = np.dot(assembly_1,all_a[i])
#     global_stiffness_matrix = global_stiffness_matrix + assembly



# print(np.transpose(np.transpose(k_all)))


c_1 = (youngs_modulus)/((1+poission_ratio)*(1-(2*poission_ratio)))
C = c_1*np.matrix([[1-poission_ratio,poission_ratio,poission_ratio,0],
                       [poission_ratio,1-poission_ratio,poission_ratio,0],
                       [poission_ratio,poission_ratio,1-poission_ratio,0],
                       [0,0,0,((1-2*poission_ratio)/2)]])


k_all = np.zeros((nelm,isoparametric_edge*2,isoparametric_edge*2))
k = np.zeros((nelm,isoparametric_edge*2,isoparametric_edge*2))
for i in range(nelm):
    derivative_N = 1/4*np.matrix([[-(1-E2),(1-E2),(1+E2),-(1+E2)],
                                  [-(1-E1),-(1+E1),(1+E1),(1-E1)]])
    x_y_ele = all_ele_coord[i]
    jacobi_1 = derivative_N*x_y_ele
    jacobi_inverse = np.linalg.inv(jacobi_1)
    B_1_matrix = jacobi_inverse*derivative_N
    for i in range(nelm):
        if i == 1 & 3:
            radius = 1.25
        else:
            radius = 3.75
    
    N_1 = N_2 = N_3 = N_4 = 1/4
    B_matrix = np.matrix([[B_1_matrix[0,0],0,B_1_matrix[0,1],0,B_1_matrix[0,2],0,B_1_matrix[0,3],0],
                  [0,B_1_matrix[1,0],0,B_1_matrix[1,1],0,B_1_matrix[1,2],0,B_1_matrix[1,3]],
                  [N_1/radius,0,N_2/radius,0,N_3/radius,0,N_4/radius,0],
                  [B_1_matrix[1,0],B_1_matrix[0,0],B_1_matrix[1,1],B_1_matrix[0,1],B_1_matrix[1,2],B_1_matrix[0,2],B_1_matrix[1,3],B_1_matrix[0,3]]])

    k_1 = 2*3.14*radius
    k_2 = np.dot(C,B_matrix)
    k_3 = np.dot(np.transpose(B_matrix),k_2)
    k_all[i] = k_1*k_3*weight*weight*np.linalg.det(jacobi_1)
    # k[i] = weight*2*3.14*radius*np.transpose(B_matrix)*C*B_matrix*np.linalg.det(jacobi_1)
    
    
    
# print((k_all[1].shape))
# # print(k)

# global_stiffness_matrix = 0
# for i in range(nelm):
#       assembly_1 = np.dot((np.transpose(all_a[i])),k_all[i])
#       assembly = np.dot(assembly_1,all_a[i])
#       global_stiffness_matrix = global_stiffness_matrix + assembly
    #   print((assembly_1))

# print(np.around(global_stiffness_matrix))

# print(all_a)



# x_cells = int(np.sqrt(nelm))
# y_cells = int(np.sqrt(nelm))
# elements = np.arange((x_cells+1)*(y_cells+1)).reshape(y_cells+1,x_cells+1)
# for i in range(0,(x_cells+1),1):
#     if i % 2 != 0:
#         elements_1 = elements[i]
#         elements[i] = elements_1[::-1]
# # print(elements)
# ID = np.array([])
# for i in range(y_cells):
#     for j in range(x_cells):
#         id = np.matrix([elements[i,j],elements[i,j+1],elements[i+1,j+1],elements[i+1,j]])
#         ID = np.append(ID,id)

# print(ID)



# a = np.array([[0,1,2],
#               [3,4,5],
#               [6,7,8]])
# b = a[1]
# a[1] = b[::-1]
# print(a)


b= np.ones((4,8,8))
global_assembly = 0
for i in range(nelm):
    
    assembly_1 = np.dot((np.transpose(all_a[i])),k_all[i])
    assembly = np.dot(assembly_1,all_a[i])
    global_assembly = global_assembly+assembly
    # print(assembly)




# a_1 = np.matrix([[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0]])

# a_2 = np.matrix([[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]])

# a_3 = np.matrix([[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0]])

# a_4 = np.matrix([[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0]])

# assembly_1 = np.dot((np.transpose(a_1)),k_all[0])
# assembly_one = np.dot(assembly_1,a_1)

# assembly_2 = np.dot((np.transpose(a_2)),k_all[1])
# assembly_two = np.dot(assembly_2,a_2)

# assembly_3 = np.dot((np.transpose(a_3)),k_all[2])
# assembly_three = np.dot(assembly_3,a_3)

# assembly_4 = np.dot((np.transpose(a_4)),k_all[3])
# assembly_four = np.dot(assembly_4,a_4)


# # print(np.linalg.det(assembly_one+assembly_two+assembly_three+assembly_four))


# b_1 = np.matrix([[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]])

# b_2 = np.matrix([[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0]])

# b_3 = np.matrix([[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0]])

# b_4 = np.matrix([[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]])

# b_assembly_1 = np.dot((np.transpose(b_1)),k_all[0])
# b_assembly_one = np.dot(b_assembly_1,b_1)

# b_assembly_2 = np.dot((np.transpose(b_2)),k_all[1])
# b_assembly_two = np.dot(b_assembly_2,b_2)

# b_assembly_3 = np.dot((np.transpose(b_3)),k_all[2])
# b_assembly_three = np.dot(b_assembly_3,b_3)

# b_assembly_4 = np.dot((np.transpose(b_4)),k_all[3])
# b_assembly_four = np.dot(b_assembly_4,b_4)

# print(np.linalg.det(b_assembly_one+b_assembly_two+b_assembly_three+b_assembly_four))


# print(np.linalg.det(assembly_one),np.linalg.det(b_assembly_one))
# a_1 = np.matrix([[1,0,0,0,0,0,0,0,0,0,0,0],
#                  [0,1,0,0,0,0,0,0,0,0,0,0],
#                  [0,0,1,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,1,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,0,1,0],
#                  [0,0,0,0,0,0,0,0,0,0,0,1],
#                  [0,0,0,0,0,0,0,0,1,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,1,0,0]])

# a_2 = np.matrix([[0,0,1,0,0,0,0,0,0,0,0,0],
#                  [0,0,0,1,0,0,0,0,0,0,0,0],
#                  [0,0,0,0,1,0,0,0,0,0,0,0],
#                  [0,0,0,0,0,1,0,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,0,1,0,0,0],
#                  [0,0,0,0,0,0,0,0,0,1,0,0],
#                  [0,0,0,0,0,0,1,0,0,0,0,0],
#                  [0,0,0,0,0,0,0,1,0,0,0,0]])



# assembly_1 = np.dot((np.transpose(a_1)),k_all[0])
# assembly_one = np.dot(assembly_1,a_1)

# assembly_2 = np.dot((np.transpose(a_2)),k_all[1])
# assembly_two = np.dot(assembly_2,a_2)


# print((assembly_one+assembly_two))




'''# https://www.ccg.msm.cam.ac.uk/images/FEMOR_Lecture_1.pdf #'''




# a = np.dot(np.transpose(a_1),a_1)+np.dot(np.transpose(a_2),a_2)+np.dot(np.transpose(a_3),a_3)+np.dot(np.transpose(a_4),a_4)
# print(np.linalg.det(a))


# print((k_all[0]))



# print(all_a)

'''
===================================================================================================================================
                                                     four elements
===================================================================================================================================
'''

                ##1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9
a_1 = np.matrix([[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], #1
                 [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], #1
                 [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], #2
                 [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0], #2
                 [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0], #4
                 [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0], #4
                 [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0], #5
                 [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]])#5

                ##1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9
a_2 = np.matrix([[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], #2
                 [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0], #2
                 [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0], #3
                 [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0], #3
                 [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0], #5
                 [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0], #5
                 [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0], #6
                 [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0]])#6

                ##1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9
a_3 = np.matrix([[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0], #4
                 [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0], #4
                 [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0], #5
                 [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0], #5
                 [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0], #7
                 [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0], #7
                 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0], #8
                 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0]])#8
                

                ##1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9
a_4 = np.matrix([[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0], #5
                 [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0], #5
                 [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0], #6
                 [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0], #6
                 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0], #8
                 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0], #8
                 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0], #9
                 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]])#9



assembly_1 = np.dot((np.transpose(a_1)),k_all[0])
assembly_one = np.dot(assembly_1,a_1)

assembly_2 = np.dot((np.transpose(a_2)),k_all[1])
assembly_two = np.dot(assembly_2,a_2)

assembly_3 = np.dot((np.transpose(a_3)),k_all[2])
assembly_three = np.dot(assembly_3,a_3)

assembly_4 = np.dot((np.transpose(a_4)),k_all[3])
assembly_four = np.dot(assembly_4,a_4)

total = assembly_one + assembly_two + assembly_three + assembly_four
# print((total[8,8]))
# print((total[8,9]))
# print((total[9,8]))
# print((total[9,9]))
# print(np.linalg.det(total))

# b = np.ones((8,8))
# assembly_1 = np.dot(b,a_1)
# assembly_one = np.dot(np.transpose(a_1),assembly_1)

# assembly_2 = np.dot(b,a_2)
# assembly_two = np.dot(np.transpose(a_2),assembly_2)

# assembly_3 = np.dot(b,a_3)
# assembly_three = np.dot(np.transpose(a_3),assembly_3)

# assembly_4 = np.dot(b,a_4)
# assembly_four = np.dot(np.transpose(a_4),assembly_4)
# total_assembly = (assembly_one)
# print((total_assembly[8,8]))
# print((total_assembly[8,9]))
# print((total_assembly[9,8]))
# print((total_assembly[9,9]))
# print()
# print(assembly_one)
# print((assembly_one+assembly_two+assembly_three+assembly_four))
# print(np.linalg.det(assembly_one+assembly_two+assembly_three+assembly_four))
'''
(1,2,4,5),(2,3,5,6),(4,5,7,8),(5,6,8,9)
and 
(1,2,5,4),(2,3,6,5),(4,5,8,7),(5,6,9,8)
have same assembly of (18,18) with np.ones
'''



'''
===============================================================================================================================
                                     two elements
===============================================================================================================================
'''
                ##1,1,2,2,3,3,4,4,5,5,6,6
a_1 = np.matrix([[1,0,0,0,0,0,0,0,0,0,0,0], #1
                 [0,1,0,0,0,0,0,0,0,0,0,0], #1
                 [0,0,1,0,0,0,0,0,0,0,0,0], #2
                 [0,0,0,1,0,0,0,0,0,0,0,0], #2
                 [0,0,0,0,0,0,1,0,0,0,0,0], #4
                 [0,0,0,0,0,0,0,1,0,0,0,0], #4
                 [0,0,0,0,0,0,0,0,1,0,0,0], #5
                 [0,0,0,0,0,0,0,0,0,1,0,0]])#5

                ##1,1,2,2,3,3,4,4,5,5,6,6
a_2 = np.matrix([[0,0,1,0,0,0,0,0,0,0,0,0], #2
                 [0,0,0,1,0,0,0,0,0,0,0,0], #2                                                                                                  
                 [0,0,0,0,1,0,0,0,0,0,0,0], #3
                 [0,0,0,0,0,1,0,0,0,0,0,0], #3
                 [0,0,0,0,0,0,0,0,1,0,0,0], #5
                 [0,0,0,0,0,0,0,0,0,1,0,0], #5
                 [0,0,0,0,0,0,0,0,0,0,1,0], #6                                                                                                                                                                             
                 [0,0,0,0,0,0,0,0,0,0,0,1]])#6

                                                           

assembly_1 = np.dot((np.transpose(a_1)),k_all[0])                                                                                                                         
assembly_one = np.dot(assembly_1,a_1)

assembly_2 = np.dot((np.transpose(a_2)),k_all[1])
assembly_two = np.dot(assembly_2,a_2)


# print(np.linalg.det(assembly_one+assembly_two))


'''
if there any change of id(node numbering) also there is no change in the matrix still det is zero
for two elements its working properly, global stiffness is singular before reducing.

'''

'''
==============================================================================================================================================
'''

# print(np.linalg.det(k_all))
# print(np.allclose(k_all[0],k_all[3]))


a = all_ele_coord[0]

print(a[0,0])
print(a[1,0])

