# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 19:02:12 2021

@author: diego.mera
"""

##########################################################################################################################################################################
import openseespy.opensees as op
import openseespy.postprocessing.ops_vis as opsv
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
#import the os module
import os
import math
from tabulate import tabulate

import scipy.linalg
import socket, struct
import sys

op.wipe()

#%%%%%%%%%%%%%%%%%%%%%%% Define Host and Port %%%%%%%%%%%%%%%%%%%%%%%
TCP_IP = '0.0.0.0'
# TCP_IP = '127.0.0.1'
TCP_PORT = 8090
BUFFER_SIZE = 4096

#########################################################################################################################################################################
#All results in Inch, Kip and Sec
# Define ELEMENTS & SECTIONS 
inch = 1.0
kip = 1.0
sec = 1.0
LunitTXT = 'inch'
FunitTXT = 'kip'
TunitTXT = 'sec'
ft = 12*inch
ksi = kip/math.pow(inch,2)
psi = ksi/1000
lbf = psi*inch*inch
pcf = lbf/math.pow(ft,3)
inch2 = inch*inch
inch4 = math.pow(inch,4)
cm = inch/2.54
PI = 2 * math.asin(1.0)
g = 32.2 * ft/math.pow(sec,2)
Ubig = 1e10
Usmall = 1/Ubig


op.model('basic', '-ndm', 2, '-ndf', 3) 
LCol = 36.0*ft       # column length
LBeam = 42.0*ft      # beam length
Weight = 200.0*kip   # superstructure weight

# define section geometry
HCol = 5.0*ft   # Column Depth
BCol = 5.0*ft   # Column Width
HBeam = 8.0*ft  # Beam Depth
BBeam = 5.0*ft  # Beam Width

# calculated parameters
PCol =Weight/2  # nodal dead-load weight per column
g = 386.4       # Gravity in/s2
Mass =  PCol/g  # mass lb
MCol = ((Weight/LBeam)*math.pow(LBeam,2))/12

# calculated geometry parameters
ACol = HCol*BCol                        # cross-sectional area column
ABeam = HBeam*BBeam                     # cross-sectional area beam
IzCol = (BCol*math.pow(HCol,3))/12      # Column moment of inertia
IzBeam =  (BBeam*math.pow(HBeam,3))/12  # Beam moment of inertia


# Nodes
op.node(1, 0.0, 0.0)
op.node(2, LBeam, 0.0);
# op.node(3, 2*LBeam, 0.0);
op.node(4, 0.0, LCol)
op.node(5, LBeam, LCol)
op.node(6, 2*LBeam, LCol)
op.node(7, 0.0, 2*LCol)
op.node(8, LBeam, 2*LCol)
op.node(9, 2*LBeam, 2*LCol)
op.node(10, 0.0, 3*LCol)
op.node(11, LBeam, 3*LCol)
op.node(12, 2*LBeam, 3*LCol)

# Create a homogeneous SP constriant.
op.fix(1, 1, 1, 0)
op.fix(2, 1, 1, 0); 
op.fix(6, 0, 1, 0); 

IDctrlNode = 2
IDctrlDOF = 1

# set the mass at a node
op.mass(4, Mass, 0.0, 0.0)
op.mass(5, Mass, 0.0, 0.0)
op.mass(6, Mass, 0.0, 0.0)
op.mass(7, Mass, 0.0, 0.0)
op.mass(8, Mass, 0.0, 0.0)
op.mass(9, Mass, 0.0, 0.0)
op.mass(10, Mass, 0.0, 0.0)
op.mass(11, Mass, 0.0, 0.0)
op.mass(12, Mass, 0.0, 0.0)

ColSecTag = 1			 # assign a tag number to the column section
BeamSecTag = 2           # assign a tag number to the beam section
	
seccion = 'Hormigon'

if seccion == 'Hormigon':
    coverCol = 6.0*inch      # Column cover to reinforcing steel NA.
    numBarsCol = 10          # number of longitudinal-reinforcement bars in column. (symmetric top & bot)
    barAreaCol = 2.25*inch2  # area of longitudinal-reinforcement bars
    
    # MATERIAL parameters
    IDconcU = 1 			 # material ID tag -- unconfined cover concrete (here used for complete section)
    IDreinf = 2 			 # material ID tag -- reinforcement
    
    # nominal concrete compressive strength
    fc = -4.0*ksi 				   # CONCRETE Compressive Strength (+Tension, -Compression)
    Ec = 57*ksi*math.sqrt(-fc/psi) # Concrete Elastic Modulus (the term in sqr root needs to be in psi
    
    # unconfined concrete
    fc1U = fc			    # UNCONFINED concrete (todeschini parabolic model), maximum stress
    eps1U = -0.003			# strain at maximum strength of unconfined concrete
    fc2U =  0.2*fc1U		# ultimate stress
    eps2U = -0.05			# strain at ultimate stress
    Lambda = 0.1		    # ratio between unloading slope at $eps2 and initial slope $Ec
    
    # tensile-strength properties
    ftU = -0.14* fc1U		# tensile strength +tension
    Ets = ftU/0.002			# tension softening stiffness
    
    Fy = 66.8*ksi			# STEEL yield stress
    Es = 29000.0*ksi		# modulus of steel
    Bs = 0.01				# strain-hardening ratio 
    R0 = 18.0				# control the transition from elastic to plastic branches
    cR1 = 0.925				# control the transition from elastic to plastic branches
    cR2 = 0.15				# control the transition from elastic to plastic branches
    
    op.uniaxialMaterial('Concrete02', IDconcU, fc1U, eps1U, fc2U, eps2U, Lambda, ftU, Ets) # build cover concrete (unconfined)
    op.uniaxialMaterial('Steel02', IDreinf, Fy, Es, Bs, R0,cR1,cR2) # build reinforcement material
    # FIBER SECTION properties -------------------------------------------------------------
    # symmetric section
    #                        y
    #                        ^
    #                        |     
    #             ---------------------     --   --
    #             |   o     o     o    |     |    -- cover
    #             |                    |     |
    #             |                    |     |
    #    z <---   |          +         |     H
    #             |                    |     |
    #             |                    |     |
    #             |   o     o     o    |     |    -- cover
    #             ---------------------     --   --
    #             |-------- B --------|
    #
    # RC section: 
    coverY = HCol/2.0	# The distance from the section z-axis to the edge of the cover concrete -- outer edge of cover concrete
    coverZ = BCol/2.0	# The distance from the section y-axis to the edge of the cover concrete -- outer edge of cover concrete
    coreY = coverY-coverCol
    coreZ = coverZ-coverCol
    nfY = 16        # number of fibers for concrete in y-direction
    nfZ = 4			# number of fibers for concrete in z-direction
    
    op.section('Fiber', ColSecTag)
    op.patch('quad', IDconcU, nfZ, nfY, -coverY,coverZ, -coverY,-coverZ, coverY,-coverZ, coverY,coverZ) # Define the concrete patch
    op.layer('straight', IDreinf, numBarsCol, barAreaCol, -coreY,coreZ,-coreY,-coreZ)
    op.layer('straight', IDreinf, numBarsCol, barAreaCol, coreY,coreZ, coreY,-coreZ)
    
    # BEAM section:
    op.section('Elastic', BeamSecTag,Ec,ABeam,IzBeam)	# elastic beam section)

elif seccion == 'Agregator':
    #Define Elements and Sections
    ColMatTagFlex  = 2
    ColMatTagAxial = 3
    
    Fy = 66.8*ksi			# STEEL yield stress
    E0 = 29000.0*ksi		# modulus of steel
    b = 0.01
    R0=18.5
    cR1=0.925
    cR2=0.15
    params=[R0,cR1,cR2]
    a2 = 1.0
    a1=a2*Fy/E0
    a4 = 1.0
    a3=a4*Fy/E0
    
    op.uniaxialMaterial('Steel02', ColMatTagFlex, Fy, E0, b, R0, cR1, cR2, a1, a2, a3, a4) #steel moment curvature isused for Mz of the section only, # bilinear behavior for flexure
    op.uniaxialMaterial('Elastic', ColMatTagAxial, E0) # this is not used as a material, this is an axial-force-strain response
    op.section('Aggregator', ColSecTag, ColMatTagAxial, 'P', ColMatTagFlex, 'Mz')  # combine axial and flexural behavior into one section (no P-M interaction here)

    # BEAM section:
    op.section('Elastic', BeamSecTag,E0,ABeam,IzBeam)	# elastic beam section)


ColTransfTag = 1
BeamTransfTag = 2
op.geomTransf('Linear', ColTransfTag)
op.geomTransf('Linear', BeamTransfTag)

numIntgrPts = 5

IntegTagCol = 1
IntegTagBeam = 2

# We are using gauss-Legendre  integration as it is the default integration scheme used in opensees tcl
op.beamIntegration('Legendre', IntegTagCol , ColSecTag, 3)  
op.beamIntegration('Legendre', IntegTagBeam , BeamSecTag, 3)

# Select the section Type
sectiontype = 'elastic' # elastic - nonlinear - disp - force

if sectiontype == 'nonlinear':
    # columns
    op.element('nonlinearBeamColumn', 1, 1, 4, numIntgrPts, ColSecTag, ColTransfTag)
    op.element('nonlinearBeamColumn', 2, 2, 5, numIntgrPts, ColSecTag, ColTransfTag)
    op.element('nonlinearBeamColumn', 3, 4, 7, numIntgrPts, ColSecTag, ColTransfTag)
    op.element('nonlinearBeamColumn', 4, 5, 8, numIntgrPts, ColSecTag, ColTransfTag)
    op.element('nonlinearBeamColumn', 5, 6, 9, numIntgrPts, ColSecTag, ColTransfTag)
    op.element('nonlinearBeamColumn', 6, 7, 10, numIntgrPts, ColSecTag, ColTransfTag)
    op.element('nonlinearBeamColumn', 7, 8, 11, numIntgrPts, ColSecTag, ColTransfTag)
    op.element('nonlinearBeamColumn', 8, 9, 12, numIntgrPts, ColSecTag, ColTransfTag)
    # beam
    op.element('nonlinearBeamColumn', 9, 4, 5, numIntgrPts, BeamSecTag, BeamTransfTag)
    op.element('nonlinearBeamColumn', 10, 5, 6, numIntgrPts, BeamSecTag, BeamTransfTag)
    op.element('nonlinearBeamColumn', 11, 7, 8, numIntgrPts, BeamSecTag, BeamTransfTag)
    op.element('nonlinearBeamColumn', 12, 8, 9, numIntgrPts, BeamSecTag, BeamTransfTag)
    op.element('nonlinearBeamColumn', 13, 10, 11, numIntgrPts, BeamSecTag, BeamTransfTag)
    op.element('nonlinearBeamColumn', 14, 11, 12, numIntgrPts, BeamSecTag, BeamTransfTag)
    
elif sectiontype == 'disp':
    # columns
    op.element('dispBeamColumn', 1, 1, 4, ColTransfTag, IntegTagCol)
    op.element('dispBeamColumn', 2, 2, 5, ColTransfTag, IntegTagCol)
    op.element('dispBeamColumn', 3, 4, 7, ColTransfTag, IntegTagCol)
    op.element('dispBeamColumn', 4, 5, 8, ColTransfTag, IntegTagCol)
    op.element('dispBeamColumn', 5, 6, 9, ColTransfTag, IntegTagCol)
    op.element('dispBeamColumn', 6, 7, 10, ColTransfTag, IntegTagCol)
    op.element('dispBeamColumn', 7, 8, 11, ColTransfTag, IntegTagCol)
    op.element('dispBeamColumn', 8, 9, 12, ColTransfTag, IntegTagCol)
    # beam
    op.element('dispBeamColumn', 9, 4, 5, BeamTransfTag, IntegTagBeam)
    op.element('dispBeamColumn', 10, 5, 6, BeamTransfTag, IntegTagBeam)
    op.element('dispBeamColumn', 11, 7, 8, BeamTransfTag, IntegTagBeam)
    op.element('dispBeamColumn', 12, 8, 9, BeamTransfTag, IntegTagBeam)
    op.element('dispBeamColumn', 13, 10, 11, BeamTransfTag, IntegTagBeam)
    op.element('dispBeamColumn', 14, 11, 12, BeamTransfTag, IntegTagBeam)
    
elif sectiontype == 'force':
    # columns
    # op.element('forceBeamColumn', eleTag, *eleNodes, transfTag, integrationTag, '-iter', maxIter=10, tol=1e-12, '-mass', mass=0.0)
    op.element('forceBeamColumn', 1, 1, 4, ColTransfTag, IntegTagCol)
    op.element('forceBeamColumn', 2, 2, 5, ColTransfTag, IntegTagCol)
    op.element('forceBeamColumn', 3, 4, 7, ColTransfTag, IntegTagCol)
    op.element('forceBeamColumn', 4, 5, 8, ColTransfTag, IntegTagCol)
    op.element('forceBeamColumn', 5, 6, 9, ColTransfTag, IntegTagCol)
    op.element('forceBeamColumn', 6, 7, 10, ColTransfTag, IntegTagCol)
    op.element('forceBeamColumn', 7, 8, 11, ColTransfTag, IntegTagCol)
    op.element('forceBeamColumn', 8, 9, 12, ColTransfTag, IntegTagCol)
    # beam
    op.element('forceBeamColumn', 9, 4, 5, BeamTransfTag, IntegTagBeam)
    op.element('forceBeamColumn', 10, 5, 6, BeamTransfTag, IntegTagBeam)
    op.element('forceBeamColumn', 11, 7, 8, BeamTransfTag, IntegTagBeam)
    op.element('forceBeamColumn', 12, 8, 9, BeamTransfTag, IntegTagBeam)
    op.element('forceBeamColumn', 13, 10, 11, BeamTransfTag, IntegTagBeam)
    op.element('forceBeamColumn', 14, 11, 12, BeamTransfTag, IntegTagBeam)
elif sectiontype == 'elastic':
    # columns
    op.element('elasticBeamColumn', 1, 1, 4, ACol, Es, IzCol, ColTransfTag)
    op.element('elasticBeamColumn', 2, 2, 5, ACol, Es, IzCol, ColTransfTag)
    op.element('elasticBeamColumn', 3, 4, 7, ACol, Es, IzCol, ColTransfTag)
    op.element('elasticBeamColumn', 4, 5, 8, ACol, Es, IzCol, ColTransfTag)
    op.element('elasticBeamColumn', 5, 6, 9, ACol, Es, IzCol, ColTransfTag)
    op.element('elasticBeamColumn', 6, 7, 10, ACol, Es, IzCol, ColTransfTag)
    op.element('elasticBeamColumn', 7, 8, 11, ACol, Es, IzCol, ColTransfTag)
    op.element('elasticBeamColumn', 8, 9, 12, ACol, Es, IzCol, ColTransfTag)
    # beam
    op.element('elasticBeamColumn', 9, 4, 5, ABeam, Es, IzBeam, BeamTransfTag)
    op.element('elasticBeamColumn', 10, 5, 6, ABeam, Es, IzBeam, BeamTransfTag)
    op.element('elasticBeamColumn', 11, 7, 8, ABeam, Es, IzBeam, BeamTransfTag)
    op.element('elasticBeamColumn', 12, 8, 9, ABeam, Es, IzBeam, BeamTransfTag)
    op.element('elasticBeamColumn', 13, 10, 11, ABeam, Es, IzBeam, BeamTransfTag)
    op.element('elasticBeamColumn', 14, 11, 12, ABeam, Es, IzBeam, BeamTransfTag)
else:
    print('Undefined section')


#defining gravity loads
WzBeam = Weight/LBeam
op.timeSeries('Linear', 1)
op.pattern('Plain', 1, 1)
op.eleLoad('-ele', [9, 10, 11, 12, 13, 14], '-type', '-beamUniform', -WzBeam, 0.0, 0.0)
Tol = 1e-8                       # convergence tolerance for test
op.integrator('LoadControl', 1)  # determine the next time step for an analysis
op.numberer('Plain')             # renumber dof's to minimize band-width (optimization), if you want to
op.system('BandGeneral')         # how to store and solve the system of equations in the analysis
op.constraints('Plain')          # how it handles boundary conditions
op.test('NormDispIncr', Tol, 6)  # determine if convergence has been achieved at the end of an iteration step
op.algorithm('Newton')           # use Newton's solution algorithm: updates tangent stiffness at every iteration
op.analysis('Static')            # define type of analysis static or transient
op.analyze(1)                    # apply gravity

op.loadConst('-time', 0.0)       # maintain constant gravity loads and reset time to zero


#%%   SEISMIC RECORD APPLIED IN A FOR LOOP
op.wipeAnalysis()

#applying Dynamic Ground motion analysis
Tol = 1e-8
maxNumIter = 10
GMdirection = 1
GMfact = 1.5
GMfatt = g*GMfact

Lambda = op.eigen('-fullGenLapack', 9) # eigenvalue mode 1
Omega = math.pow(Lambda[0], 0.5)
omega = np.sqrt(Lambda)
period = 2.0*np.pi/omega
freq = 1/period
betaKcomm = 2 * (0.02/Omega)

Damping = 'Modal' # Modal - Rayleigh

if Damping == 'Rayleigh':
    xDamp = 0.05				                        # 5% damping ratio
    alphaM = 0.0				                        # M-prop. damping; D = alphaM*M	
    betaKcurr = 0.0		                                # K-proportional damping;      +beatKcurr*KCurrent
    betaKinit = 0.0                                     # initial-stiffness proportional damping      +beatKinit*Kini
    op.rayleigh(alphaM,betaKcurr, betaKinit, betaKcomm) # RAYLEIGH damping
elif Damping == 'Modal':    
    zeta = 0.05
    op.modalDamping(zeta) # Modal Damping
else:
    print('Unspecified damping')

######################### Load sismic record #################################

# register = np.loadtxt('elcentro_0.01.txt',skiprows=0,usecols=[0])
# dt = 0.01
# nPts = 3119

register = np.loadtxt('elcentro_0.02.txt',skiprows=0,usecols=[0])
dt = 0.02
nPts = 1560

# register = np.loadtxt('elcentro_0.04.txt',skiprows=0,usecols=[0])
# dt = 0.04
# nPts = 780

# register = np.loadtxt('elcentro_0.05.txt',skiprows=0,usecols=[0])
# dt = 0.05
# nPts = 624

# register = np.loadtxt('elcentro_0.1.txt',skiprows=0,usecols=[0])
# dt = 0.1
# nPts = 312




i=1
motion = register[0:-1]
values = list(-1 * motion)  # should be negative
op.timeSeries('Path', i+20000, '-dt', dt, '-values', *values, '-factor', GMfatt)
IDloadTag = i+4000		# load tag
op.pattern('UniformExcitation', IDloadTag, GMdirection, '-accel', i+20000)

index = list(np.linspace(0, (len(motion)-1)*dt, len(motion)))

time2 = [0.0]

u44 = [0.0]
u77= [0.0]

okk = 0


s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind((TCP_IP, TCP_PORT))
print('Waiting for Simulink to start')
s.listen(1)
conn, addr = s.accept()
print('Connection address: ', addr)






op.constraints('Transformation')
# op.numberer('RCM')
op.numberer('Plain')
op.system('BandGeneral')
op.test('EnergyIncr', Tol, maxNumIter)
# op.test('FixedNumIter',maxNumIter)
op.algorithm('Linear')
op.integrator('AlphaOS',0.9)
# gamma = 0.5
# beta = 0.25
# op.integrator('Newmark', gamma, beta)
op.analysis('Transient')

for i in range(0,nPts-1):  
    print('Vamos en el paso : ', i)
    # receive data from Simulink
    # data = conn.recv(BUFFER_SIZE)
    recv_msg, send_addr = conn.recvfrom(BUFFER_SIZE)
    # print('Los datos recibidos son: ', recv_msg)
    # print('buffer de los datos recibidos es: ', len(recv_msg))
    # datos = struct.unpack('>d', data)
    if len(recv_msg)==8:
        Px = struct.unpack('>d', recv_msg)[0]
        print('Los datos recibidos son: ', Px)
    else:
        Px = 0
        print('Los datos recibidos son: ', 0)
    
    # Px = 100+i*1
    
    # if i>200:
    #     Px = 0
        
    op.timeSeries('Constant', i+2)
    op.pattern('Plain', i+2, i+2)
    op.load(6, Px, 0., 0.)
    
    okk = op.analyze(1, dt)
    op.remove('loadPattern', i+2)
    op.remove('timeSeries', i+2)
    
    # send data to Simulink
    msg1 = struct.pack('>d', float(op.nodeDisp(6,1)))
    conn.send(msg1)
    # print('buffer de los datos enviados es: ', len(msg1))
    print('sent data:', op.nodeDisp(6,1)) 
    
    if okk == 0 :
        tCurrent = dt*i              
        time2.append(tCurrent)
        u44.append(op.nodeDisp(4,1))
        u77.append(op.nodeDisp(7,1))

# 1. plot model with tag lebels
opsv.plot_model()

# 2. plot deformed model
sfac = 15

plt.figure()
# plot_defo with optional arguments
# sfac = opsv.plot_defo()
opsv.plot_defo(sfac, fmt_interp='b.-')
opsv.plot_defo(sfac, 3, endDispFlag=0, fmt_interp='r.--')
opsv.plot_defo(sfac, 2, fmt_interp='k.-')
opsv.plot_defo(sfac)


plt.figure()
plt.plot(time2, u77)
plt.ylabel('Horizontal Displacement of node 7 (in)')
plt.xlabel('Time (s)')
# plt.savefig('Horizontal Disp at Node 3 vs time.jpeg', dpi = 500)
plt.grid()


plt.figure()
plt.plot(time2, u44)
plt.ylabel('Horizontal Displacement of node 4 (in)')
plt.xlabel('Time (s)')
# plt.savefig('Horizontal Disp at Node 4 vs time.jpeg', dpi = 500)
plt.grid()

# plt.show()


# Displat Period and Frequency of the model
d = [ ["Mode 1", period[0], freq[0]],
      ["Mode 2", period[1], freq[1]],
      ["Mode 3", period[2], freq[2]],
      ["Mode 4", period[3], freq[3]],
      ["Mode 5", period[4], freq[4]],
      ["Mode 6", period[5], freq[5]],
      ["Mode 7", period[6], freq[6]],
      ["Mode 8", period[7], freq[7]],
      ["Mode 9", period[8], freq[8]]]

print(tabulate(d, headers=["Mode", "Period [s]", "Frequency [Hz]"]))

op.wipe()
