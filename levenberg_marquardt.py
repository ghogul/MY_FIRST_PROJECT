import numpy as np
import matplotlib.pyplot as plt


filename = "STRESS STRAIN.txt"

exp_strain = np.loadtxt(filename, usecols=(1))
exp_stress = np.loadtxt(filename, usecols=(0))

#exp_strain = exp_strain[12:-7]
#exp_stress = exp_stress[12:-7]

E = (exp_stress[3]-exp_stress[2])/(exp_strain[3]-exp_strain[2])
print(E)
#sim_stress = sim_stress[5:]
#sim_strain = sim_strain[5:]

# print(exp_strain.size)
# print(exp_stress.size)
# print(sim_stress.size)
# print(sim_strain.size)
#print(sim_stress-exp_stress)
# plt.plot(exp_strain,exp_stress)
# plt.plot(sim_strain,sim_stress)
# plt.show()
# error = np.zeros((len(exp_stress)))
# for i in range(len(exp_stress)):
#     error[i] = np.sqrt(((exp_stress[i] - sim_stress[i])/exp_stress[i])**2)
#     print(error[i])
# print((np.sum(error[1:])))

sigma_0 = 400
sigma_1 = 200
sigma_2 = 200
n       = 20

x_0 = np.array([sigma_0,sigma_1,sigma_2,n])
epsilon = exp_strain[12:]
stress = exp_stress[12:]
# print(x_0)
'''============================================================================================'''
def voce_fun(x,epsilon):
    sigma_0 = x[0]
    sigma_1 = x[1]
    sigma_2 = x[2]
    n       = x[3]
    sigma_yield = sigma_0 + (sigma_1*epsilon) + (sigma_2*(1-np.exp(-n*epsilon)))
    return sigma_yield
# voce = voce_fun(x_0,epsilon)
# print(voce)
# print(epsilon)
# plt.plot(epsilon,voce)
# plt.plot(epsilon,stress)
# plt.show()
'''============================================================================================'''
# def obj_fun(x,epsilon,stress):
#     m = len(epsilon)
#     mse = 1/m*((voce_fun(x,epsilon)-stress)**2).sum()
#     return mse
'''============================================================================================'''
def jacobian(x,epsilon):
    sigma_0 = x[0]
    sigma_1 = x[1]
    sigma_2 = x[2]
    n       = x[3]
    no_of_data = len(epsilon)
    no_of_para = len(x)
    a = np.exp(-n*epsilon)
    b = 1 - a
    partial_sigma_0 = np.ones(no_of_data)
    partial_sigma_1 = epsilon
    partial_sigma_2 = b
    partial_n       = sigma_2*epsilon*a
    jacobian = np.array([partial_sigma_0,partial_sigma_1,partial_sigma_2,partial_n])
    return jacobian


# jaco = jacobian(x_0,epsilon)
# print(jaco)
'''============================================================================================'''

def quadratic_model(delta_x,f_0,gradient_0,hessian_0,xs=None):
    xs = np.sqrt(xs/xs.min())
    delta_x = delta_x/xs
    f_1 = f_0 + np.dot(gradient_0,delta_x)+0.5*np.dot(delta_x,np.dot(hessian_0,delta_x))
    return f_1
'''============================================================================================'''

lamda = 0.001
f_tol = 1e-6
eta = 0.25
max_iteration = 100
lamda_iteration = 10

# epsilon = exp_strain[12:]
# stress = exp_stress[12:]
#print(epsilon)
#print(stress)

def levenberg_marquardt(x,epsilon,stress,lamda, f_tol, eta, max_iteration, lamda_iteration):
    no_of_data = len(epsilon)
    no_of_para = len(x)

    xs = x

    standard_deviation = np.std(stress)
    x_0 = x
    residual_0 = (voce_fun(x_0,epsilon)-stress)/standard_deviation
    f_0 = 0.5*(np.linalg.norm(residual_0)**2)
    jacobian_0 = jacobian(x_0,epsilon)/standard_deviation 
    gradient_0 = np.inner(jacobian_0,residual_0)
    norm_gradient_0 = np.linalg.norm(gradient_0)
    print(standard_deviation)
    print(x)
    print(residual_0)
    print(f_0)
    print(jacobian_0)
    print(gradient_0)
    print(norm_gradient_0)
    i = 1
    j = 1
    d_f = 1
    while (abs(d_f) > f_tol) and (i <= max_iteration) and (j <= lamda_iteration):
        square_root_lamda = np.sqrt(lamda)
        square_root_lamda_identity = square_root_lamda*np.eye(no_of_para)

        jacobian_square_root_lamda_identity = np.hstack([jacobian_0,square_root_lamda_identity])

        z = np.zeros(no_of_para) 
        residual_z = np.hstack([residual_0,z])

        (delta_x,residual,rank,sv) = np.linalg.lstsq(jacobian_square_root_lamda_identity.T,-residual_z, rcond=-1)
        norm_delta_x = np.linalg.norm(delta_x)

        x_1 = x_0 + delta_x
        residual_1 = (voce_fun(x_1,epsilon)-stress)/standard_deviation
        f_1 = 0.5*(np.linalg.norm(residual_1)**2)
        jacobian_1 = jacobian(x_1,epsilon)/standard_deviation
        gradient_1 = np.inner(jacobian_1,residual_1)
        norm_gradient = np.linalg.norm(gradient_1)

        gradient_0 = np.inner(jacobian_0,residual_0)
        hessian_0 = np.inner(jacobian_0,jacobian_0)
        m_0 = f_0
        m_1 = quadratic_model(delta_x,f_0,gradient_0,hessian_0,xs=xs)

        a = (f_0 - f_1)/(m_0 - m_1)
        d_f = f_1 - f_0


        print('='*80)
        print('iteration',i,'lamda iteration',j)
        print('lamda',lamda)
        print('x_0',x_0)
        print('f_0',f_0)
        print('m_0',m_0)
        print('.'*50)
        print('x_1',x_1)
        print('f_1',f_1)
        print('m_1',m_1)
        print('.'*50)
        print('a',a)
        print('.'*50)
        print('delta_x',delta_x)
        print('norm_delta_x',norm_delta_x)
        print('df',d_f)
        print('norm_gradient',norm_gradient)
        print('='*80)

        if a < 0.25:
            lamda = lamda*10
            print('increase lamda to',lamda)

        else:
            if a > 0.75:
                lamda = lamda/10
                print('decrease lamda to',lamda)
            else:
                print('keep lamda',lamda)

        if a > eta:
            print('step accepted')
            x_0 = x_1
            residual_0 = residual_1
            f_0 = f_1
            jacobian_0 = jacobian_1 
            j = 1
            i += 1

        else:
            print('step rejected')
            d_f = 1
            j += 1

        if i > max_iteration:
            print('maximum number of i iteration is reached')
        if j > lamda_iteration:
            print('maximum number of j iteration is reached')

    print('x_opt',x_1)
    return x_1

'''============================================================================================'''

#levenberg_marquardt = levenberg_marquardt(x_0,epsilon,stress,lamda, f_tol, eta, max_iteration, lamda_iteration)


