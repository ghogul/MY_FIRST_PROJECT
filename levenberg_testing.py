import numpy as np
import matplotlib.pyplot as plt

np.random.seed()


'''
=================================================================================================================================================================================
                                             Function of true polynomial and creating noisy data
=================================================================================================================================================================================
'''

def poly_fun():
    '''
    =====================================================================
    defined a polynomial function
    ---------------------------------------------------------------------
    size : number of data points in X-axis
    X : array of data points in X-axis
    y_true : defining a function with known co-efficients
    y_noise : for algorithm we need a noisy data here we are creating a noisy data
    ---------------------------------------------------------------------
    return
    y_noise  
    x 
    y_true 
    =====================================================================
    '''
    size = 80
    x = np.linspace(1,2,size)
    y_true = (0.2*x**3)+(0.5*x**2)+(0.8*x)
    noise = np.random.normal(0,0.08,size)
    y_noise = y_true + noise
    
    return y_noise,x,y_true
'''
=================================================================================================================================================================================
'''

'''
=================================================================================================================================================================================
                                                Function for jacobian
=================================================================================================================================================================================
'''
def jacobian(y,x):
    '''
    =====================================================================
    jacobian of a function
    ---------------------------------------------------------------------
    partial_1 : partial derivative of first term
    partial_2 : partial derivative of second term
    partial_3 : partial derivative of third term
    ---------------------------------------------------------------------
    return
    jacobian : array of partial deriavtives of the function
    =====================================================================
    '''

    partial_1 = x**3
    partial_2 = x**2
    partial_3 = x
    jacobian = np.array([partial_1,partial_2,partial_3])
    return jacobian
'''
=================================================================================================================================================================================
'''




'''
=================================================================================================================================================================================
                                                Function for model_function(m)
=================================================================================================================================================================================
'''
def quadratic_model(delta_y,f_0,gradient_0,hessian_0,ys=None):
    '''
    ======================================================================
    quadratic model
    ----------------------------------------------------------------------
    delta_y : step
    f_0 : function value f(delta_y=0)
    gradient_0 : gradient g(delta_y=0)
    hessian_0 : Hessian h(delta_y=0)
    ----------------------------------------------------------------------
    return
    m_1 : function f(delta_y)
    ======================================================================
    '''
    m_1 = f_0 + np.dot(gradient_0,delta_y)+0.5*np.dot(delta_y,np.dot(hessian_0,delta_y))
    return m_1

'''
=================================================================================================================================================================================
'''



'''
=================================================================================================================================================================================
                                                polynomial Function for algorithm
=================================================================================================================================================================================
'''
def poly_fun_unknown(y,x):
    '''
    =====================================================================
    This function defined with unkown co-efficient 
    ---------------------------------------------------------------------
    a : unkown co-efficient of the first term
    b : unkown co-efficient of the second term
    c : unkown co-efficient of the third term
    ---------------------------------------------------------------------
    return
    y : algorithm computed co=efficient of a function
    =====================================================================
    '''
    a = y[0]
    b = y[1]
    c = y[2]
    y = (a*x**3)+(b*x**2)+(c*x)
    return y
'''
=================================================================================================================================================================================
'''



'''
=================================================================================================================================================================================
                                                Levenberg Marquardt algorithm
=================================================================================================================================================================================
'''


def levenberg_marquardt(y,poly_fun_unknown,jacobian,x,y_noise,lamda,f_tol,eta,max_iteration,lamda_iteration):
    '''
    ======================================================================
    Levenberg-Marquardt method for solving 
    nonlinear least square problems
    ----------------------------------------------------------------------
    y : start parameter np.array([y1,...,yn])
    poly_fun_unknown : function for for computing residuals r=poly_fun_unknown(y,x)-y
    jacobian : function for jacobian which returns np.array([[dri/dyj,...,drm/dyn)
    x : data points np.array([x1,...,xm])
    y_noise : data values np.array([y_noise1,...,y_noisem])
    lamda : initial value of lamda
    ftol : termination criterion
    eta : accuracy to accept step
    max_iteration : maximum_number of iterations
    lambda_itereration : maximum_number of lambda adjustions
    ----------------------------------------------------------------------
    return:
    y_1 : algorithm computed co-efficients
    ======================================================================
    '''
    
    
    
    no_of_data = len(x)
    no_of_para = len(y)
    

    ys = y

    standard_deviation = np.std(y_noise)
    y_0 = y
    
    residual_0 = (poly_fun_unknown(y_0,x) - y_noise)/standard_deviation
    
    f_0 = 0.5*(np.linalg.norm(residual_0)**2)
    jacobian_0 = jacobian(y_0,x)/standard_deviation
    gradient_0 = np.inner(jacobian_0,residual_0)
    norm_gradient_0 = np.linalg.norm(gradient_0)

    i = 1
    j = 1
    delta_f = 1
    loop = 1
    while ((abs(delta_f) > f_tol) and (i <= max_iteration) and (j <= lamda_iteration)):
        f.write("="*80)
        f.write("\n")
        f.write(".."*10)
        f.write("\n")
        f.write("loop = ")
        f.write("\t")
        f.write(str(loop))
        f.write("\n")
        f.write(".."*10)
        f.write("\n")

        square_root_lamda = np.sqrt(lamda)
        square_root_lamda_identity = square_root_lamda*np.eye(no_of_para)

        jacobian_square_root_lamda_identity = np.hstack([jacobian_0,square_root_lamda_identity])

        z = np.zeros(no_of_para) 
        residual_z = np.hstack([residual_0,z])

        (delta_y,residual,rank,sv) = np.linalg.lstsq(jacobian_square_root_lamda_identity.T,-(residual_z), rcond=-1)
        norm_delta_y = np.linalg.norm(delta_y)

        y_1 = y_0 + delta_y
        

        residual_1 = (poly_fun_unknown(y_1,x) - y_noise)/standard_deviation
        f_1 = 0.5*(np.linalg.norm(residual_1)**2)
        jacobian_1 = jacobian(y_1,x)/standard_deviation
        gradient_1 = np.inner(jacobian_1,residual_1)
        norm_gradient_1 = np.linalg.norm(gradient_1)

        gradient_0 = np.inner(jacobian_0,residual_0)
        #print('gradient_0',gradient_0)
        hessian_0 = np.inner(jacobian_0,jacobian_0)
        #print('hessian_0',hessian_0)
        m_0 = f_0
        m_1 = quadratic_model(delta_y,f_0,gradient_0,hessian_0,ys=ys)

        a = np.divide((f_0 - f_1),(m_0 - m_1),out=np.zeros_like((m_0 - m_1)),where=(m_0 - m_1)!=0)
        

        delta_f = f_1 - f_0

        if a <= 0.25:
            lamda = lamda*10
            f.write('lamda is increased to = ')
            f.write("")
            f.write(str(lamda))
            f.write("\n")

        else:
            if a > 0.75:
                lamda = lamda/10
                f.write('lamda is decreased to = ')
                f.write("")
                f.write(str(lamda))
                f.write("\n")
            # else:
            #     print('keeping same lamda',lamda)

        if a > eta:
            f.write('step accepted')
            f.write("\n")
            y_0 = y_1
            residual_0 = residual_1
            f_0 = f_1
            jacobian_0 = jacobian_1 
            j = 1
            i += 1

        else:
            f.write('step rejected')
            f.write("\n")
            delta_f = 1
            j += 1

        if i > max_iteration:
            f.write('maximum number of iteration is reached')
        if j > lamda_iteration:
            f.write('maximum number of lamda iteration is reached')

        

        
        loop = loop +1
        
    f.write('='*80)
    f.write("\n")
    f.write('Algorithim computed coeeficient')
    f.write("\t")
    f.write(str(y_1))
    
    
     
    return y_1

    


'''
=================================================================================================================================================================================
'''



'''
=================================================================================================================================================================================
                                                calling algorithm and writing an output file
=================================================================================================================================================================================
'''


'''y : initial guess for the algorithm'''
y = np.array([10,20,30])

'''Calling a function which generates data for algorithm and true values'''
y_noise,x,y_true = poly_fun()

'''writing an output in a file and calling algorithm'''

with open('Algorithm_testing.txt', 'w') as f:
    f.write("="*80)
    f.write("\n")
    f.write("\t\t\tprogram test for algorithm\n")
    f.write("="*80)
    f.write("\n")
    f.write(".."*20)
    f.write("\n")
    f.write("True polynomial function = 0.2X^3+0.5X^2+0.8X")
    f.write("\n")
    f.write(".."*20)
    f.write("\n")
    f.write("initial guess = ")
    f.write(str(y))
    f.write("\n")
    lm = levenberg_marquardt(y,poly_fun_unknown,jacobian,x,y_noise,lamda = 0.0001,f_tol = 1e-10,eta = 0.001,max_iteration = 100,lamda_iteration = 20)




'''
=================================================================================================================================================================================
                                                Post processing
=================================================================================================================================================================================
'''   

'''y_sim : final co-efficients computed by algorithm and for testing data and training data i am suffling the data'''


y_sim = poly_fun_unknown(lm,x)
total_shuffle = np.array([y_noise,x]).T
np.random.shuffle(total_shuffle)

training = total_shuffle[:int(len(x)*0.6)]
testing  = total_shuffle[int(len(x)*0.6):]


'''plots for the true function and algorithm computed co-efficient function'''
fig,ax = plt.subplots(2,constrained_layout=True)
ax[0].plot(x,y_true,color='blue',label='True data')
ax[0].scatter(x,y_noise,color='red',marker='.')
ax[1].plot(x,y_true,color='blue',label='True data')
ax[1].plot(x,y_sim,color='black',linestyle ='--',label='simulated data')
ax[1].scatter(x,y_noise,color='red',marker='.')
ax[0].grid()
ax[0].legend()
ax[1].grid()
ax[1].legend()
ax[0].set_title('True function curve ')
ax[0].set_xlabel('X')
ax[0].set_ylabel('Y')
ax[1].set_title('True and Algorithm function curves ')
ax[1].set_xlabel('X')
ax[1].set_ylabel('Y')
plt.savefig("Testing algorithm true and simulated value")
fig,ax = plt.subplots(2,constrained_layout=True)
ax[0].scatter(training[:,1],training[:,0],color='red',marker='.',label='Training data')
ax[0].plot(x,y_sim,color='blue',label='Predicted data')
ax[1].scatter(testing[:,1],testing[:,0],color='red',marker='.',label='Testing data')
ax[1].plot(x,y_sim,color='blue',label='Predicted data')
ax[0].grid()
ax[0].legend()
ax[1].grid()
ax[1].legend()
ax[0].set_title('Training data')
ax[0].set_xlabel('X')
ax[0].set_ylabel('Y')
ax[1].set_title('Testing data')
ax[1].set_xlabel('X')
ax[1].set_ylabel('Y')
plt.savefig("Testing algorithm training-testing datas")

'''
=================================================================================================================================================================================
'''

