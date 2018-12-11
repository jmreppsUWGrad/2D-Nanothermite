#########################################################################################################
################## TEMPORAL SCHEMES RUNGE KUTTA FOR TRANSCRITICAL COMBUSTION ############################
#########################################################################################################

import numpy as np
#import os
#import sys
#
#from pdb import set_trace as keyboard

scheme_data = {
		"RK2":{}, "RK3":{}, "RK4":{}, "RK4_CLASSICAL":{}, "RK4_LOW":{}, "RK6":{}
		}

scheme_data["RK2"]["rk_substep_fraction"] = np.array([1./2., 1./2.]).astype("float64")
scheme_data["RK2"]["rk_coeff"] = [
               			  [0.       ,     0. ],
               			  [1.       ,     0. ]
               			 ]
scheme_data["RK2"]["information"] = "This is the general RK2 method"


'''
scheme_data["RK3"]["rk_substep_fraction"] = np.array([1./6., 2./3., 1./6.]).astype("float64")
scheme_data["RK3"]["rk_coeff"] = [
               			  [0.       ,     0. ,    0. ],
               			  [1./2.    ,     0. ,    0. ],
               			  [-1       ,     2. ,    0. ]
               			 ]
scheme_data["RK3"]["information"] = "This is the general RK3 method"
'''
scheme_data["RK3"]["rk_substep_fraction"] = np.array([8./15.,2./15.,1./3.]).astype("float64")
scheme_data["RK3"]["rk_coeff"] = [
               			  [8./15. ,     0. ,     0.  ],
               			  [-17./60., 5./12. ,     0. ],
               			  [  0.     , -5./12., 3./4. ]
               			 ]
scheme_data["RK3"]["information"] = "This is Wray's low-storage method"


scheme_data["RK4"]["rk_substep_fraction"] = np.array([1./6., 1./3., 1./3., 1./6.]).astype("float64")
scheme_data["RK4"]["rk_coeff"] = [
               			  [0.   ,   0. ,    0. ,  0.],
               			  [1./2.,   0. ,    0. ,  0.],
               			  [0.   , 1./2.,    0. ,  0.],
				        [0.   ,    0.,    1. ,  0.]
               			 ]
scheme_data["RK4"]["information"] = "This is the original method"

scheme_data["RK4_CLASSICAL"]["rk_substep_fraction"] = np.array([1./8., 3./8., 3./8., 1./8.]).astype("float64")
scheme_data["RK4_CLASSICAL"]["rk_coeff"] = [
               			  [0.    ,   0. ,    0. ,  0.],
               			  [1./3. ,   0. ,    0. ,  0.],
               			  [-1./3.,   1. ,    0. ,  0.],
				        [  1.  ,  -1. ,    1. ,  0.]
               			 ]
scheme_data["RK4_CLASSICAL"]["information"] = "This is the classical 3/8ths method"

scheme_data["RK4_LOW"]["rk_substep_fraction"] = np.array([1./6., 1./3., 1./3., 1./6.]).astype("float64")
scheme_data["RK4_LOW"]["rk_coeff"] = [
					      [ 0.        ,  0.        ,  0.        ,  0.        ],
       					[ 0.69631521,  0.        ,  0.        ,  0.        ],
       					[ 0.07801568,  0.21640084,  0.        ,  0.        ],
       					[ 0.07801568,  0.0470887 ,  0.69991726,  0.        ]
				    ]
scheme_data["RK4_LOW"]["information"] = "This is the low storage 4 step method"
scheme_data["RK6"]["rk_substep_fraction"] = np.array([1./12., 0., 0., 0., 5./12., 5./12., 1./12.]).astype("float64")
scheme_data["RK6"]["rk_coeff"] = \
[
 [ 0.         ,  0.         ,  0.         ,  0.         ,  0.         ,  0.         ,  0.        ],
 [ 0.39344262 ,  0.         ,  0.         ,  0.         ,  0.         ,  0.         ,  0.        ],
 [-0.2710704  ,  0.87821326 ,  0.         ,  0.         ,  0.         ,  0.         ,  0.        ],
 [ 0.13937841 ,  0.40771164 ,  0.16989108 ,  0.         ,  0.         ,  0.         ,  0.        ],
 [ 0.15627693 ,  0.12653933 ,  0.06358578 , -0.07000884 ,  0.         ,  0.         ,  0.        ],
 [-0.04816371 , -0.33771295 , -0.17556093 ,  0.33157308 ,  0.9534713  ,  0.         ,  0.        ],
 [ 0.45943389 ,  1.0558681  ,  0.55987576 , -1.30782124 , -1.14932252 ,  1.38196601 ,  0.        ]
]

scheme_data["RK6"]["information"] = "This is the 7 step 6th order RUNGE KUTTA method 7th order accuracy "
#'''
#scheme_data["RK6_TEMP"] = {}
#scheme_data["RK6_TEMP"]["rk_substep_fraction"] = np.array([1./12., 0., 0., 0., 5./12., 5./12., 1./12.]).astype("float64")
#scheme_data["RK6_TEMP"]["rk_coeff"] = \
#[
# [ 0.         ,  0.         ,  0.         ,  0.         ,  0.         ,  0.         ,  0.        ],
# [ 0.4        ,  0.         ,  0.         ,  0.         ,  0.         ,  0.         ,  0.        ],
# [-0.54352131 ,  1.14352131 ,  0.         ,  0.         ,  0.         ,  0.         ,  0.        ],
# [ 0.34457997 ,  0.19518001 ,  0.26024002 ,  0.         ,  0.         ,  0.         ,  0.        ],
# [ 0.15838049 ,  0.12761311 ,  0.02584164 , -0.03544203 ,  0.         ,  0.         ,  0.        ],
# [-0.04240538 , -0.25308909 , -0.02093251 ,  0.16833295 ,  0.87170083 ,  0.         ,  0.        ],
# [ 0.42012444 ,  0.62737993 , -0.02454563 , -0.66445459 , -0.74047016 ,  1.38196601 ,  0.        ]
#]
#'''

#scheme_data["RK6_TEMP"]["information"] = "This is the 7 step 6th order RUNGE KUTTA method 6th order accuracy "

scheme_data["RK8"] = {}
scheme_data["RK8"]["rk_substep_fraction"] = np.array([1./20. ,0., 0., 0., 0., 13./180., 1./5., 16./45., 1./5., 13./180., 1./20.]).astype("float64")
scheme_data["RK8"]["rk_coeff"] = \
[
       [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00],
       [  5.20833333e-03,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00],
       [ -4.62366395e-01,   5.37169580e-01,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00],
       [  2.80511943e-02,   0.00000000e+00,   8.41535829e-02,
          0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00],
       [  8.18978670e+00,   0.00000000e+00,  -2.87750723e+01,
          2.12694961e+01,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00],
       [  4.03740608e-02,   0.00000000e+00,   0.00000000e+00,
          1.32188232e-01,   1.10872008e-04,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00],
       [  2.58739003e-01,   0.00000000e+00,   0.00000000e+00,
         -1.07897045e+00,   3.49556838e-01,   1.29800144e+00,
          0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00],
       [  1.95596959e-01,   0.00000000e+00,   0.00000000e+00,
         -7.71355141e-01,   8.86776408e-02,   1.01696927e+00,
         -2.98887262e-02,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00],
       [ -5.64537398e-03,   0.00000000e+00,   0.00000000e+00,
          3.66641509e-01,  -4.76781447e-02,  -2.20463266e-01,
          1.38999909e-02,   6.59184499e-02,   0.00000000e+00,
          0.00000000e+00,   0.00000000e+00],
       [ -1.61869945e+00,   0.00000000e+00,   0.00000000e+00,
          8.48593620e+00,  -1.60007462e+00,  -1.67069639e+01,
          5.67072969e-01,   2.68925147e+00,   9.01080414e+00,
          0.00000000e+00,   0.00000000e+00],
       [  8.76517117e-01,   0.00000000e+00,   0.00000000e+00,
         -4.11389408e+00,   4.72947418e-01,   1.37853746e+01,
          2.85293485e-02,  -5.92592593e-01,  -9.70629863e+00,
          2.49416793e-01,   0.00000000e+00]
]


scheme_data["RK8"]["information"] = "This is the 11 step 8th order RUNGE KUTTA method 9th order accuracy "



class runge_kutta():
    ''' Pass scheme name as a string '''
    def __init__(self, scheme):
        self.scheme = scheme.upper()
        if scheme.upper() not in scheme_data.keys():
            print "Unrecognized Runge Kutta scheme " + scheme
            print "Select from : ", sorted(self.listSupportedSchemes())
            print "Or 'Euler' "
            self.Nk=-1
            return
        self.rk_substep_fraction =   scheme_data[self.scheme]["rk_substep_fraction"]
        self.Nk			 =   self.rk_substep_fraction.size	# Total number of steps
        self.rk_coeff 		 =   scheme_data[self.scheme]["rk_coeff"]

    def listSupportedSchemes(self):
        return scheme_data.keys()
