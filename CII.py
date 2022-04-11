#!/usr/bin/env python 
from __future__ import print_function
from tempfile import tempdir
import pinocchio
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.linalg
from sys import argv
from os.path import dirname, join, abspath

np.set_printoptions(threshold=sys.maxsize)

def modelInitialize():
    global model, foot_distance, data, LFframe_id, RFframe_id, PELVjoint_id, LHjoint_id, RHjoint_id, LFjoint_id, q_init, RFjoint_id, LFcframe_id, RFcframe_id, q, qdot, qddot, LF_tran, RF_tran, PELV_tran, LF_rot, RF_rot, PELV_rot, qdot_z, qddot_z, HRR_rot_init, HLR_rot_init, HRR_tran_init, HLR_tran_init, LF_rot_init, RF_rot_init, LF_tran_init, RF_tran_init, PELV_tran_init, PELV_rot_init, CPELV_tran_init, q_command, qdot_command, qddot_command, robotIginit, q_c
    #model = pinocchio.buildModelFromUrdf("/home/jhk/legAnalysis/dyros_tocabi_with_redhands_hipheavy.urdf",pinocchio.JointModelFreeFlyer())     
    model = pinocchio.buildModelFromUrdf("/home/jhk/legAnalysis/dyros_tocabi_with_redhands.urdf",pinocchio.JointModelFreeFlyer())  
    #model = pinocchio.buildModelFromUrdf("/home/jhk/legAnalysis/atlas.urdf",pinocchio.JointModelFreeFlyer())     
    #model = pinocchio.buildModelFromUrdf("/home/jhk/legAnalysis/dyros_tocabi_with_redhands_footheavy.urdf",pinocchio.JointModelFreeFlyer())      

    pi = 3.14159265359
    data = model.createData()
    q = pinocchio.randomConfiguration(model)
    q_command = pinocchio.randomConfiguration(model)
    q_init = [0.0, 0.0, 0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0, 0, 0, 0, 0, 0.0, 0.0, 0, -0, -1.5, 0, 0, 0.0, 0, 0.0]
    for i in range(0, len(q)):
        q[i] = q_init[i]

    q_c = np.zeros((40,900))
    for i in range(0,30):
        for j in range(0,30):
            q_c[:,30*i+j] = q_init
            q_c[8,30*i+j] = pi*(j)/58 - pi/4
            q_c[9,30*i+j] = -pi*(i)/(29*3)
            q_c[10,30*i+j] = pi*(i)/(29*1.5) 
            q_c[14,30*i+j] = -pi*(j)/58 + pi/4
            q_c[15,30*i+j] = -pi*(i)/(29*3)
            q_c[16,30*i+j] = pi*(i)/(29*1.5)

    modeldof = model.nq - 7

    foot_distance = np.zeros(3)

    qdot = pinocchio.utils.zero(model.nv)
    qdot_command = pinocchio.utils.zero(model.nv)
    qddot = pinocchio.utils.zero(model.nv)
    qddot_command = pinocchio.utils.zero(model.nv)
    qdot_z = pinocchio.utils.zero(model.nv)
    qddot_z = pinocchio.utils.zero(model.nv)
    pinocchio.forwardKinematics(model, data, q, qdot, qddot)
    pinocchio.updateFramePlacements(model,data)
    pinocchio.updateGlobalPlacements(model,data)
    pinocchio.computeJointJacobians(model, data, q)
    pinocchio.computeMinverse(model, data, q)
    pinocchio.ccrba(model, data, q, qdot)

    robotIginit = data.Ig.inertia


def modelUpdate(q_desired, qdot_desired, qddot_desired):
    global contactnum, M, G, COR, Minv, b, robotJac, robotdJac, LF_j, RF_j, LF_cj, RF_cj, LF_cdj, RF_cdj, robotLambdac, robotJcinvT, robotNc, robotPc, robotmuc, robotW, robothc, LF_tran_cur, RF_tran_cur, PELV_tran_cur, COM_tran_cur, RFc_tran_cur, LFc_tran_cur, robotWinv, robotCAM, robotIg, robotCII
    pinocchio.forwardKinematics(model, data, q_desired, qdot_desired, qddot_desired)
    pinocchio.updateFramePlacements(model,data)
    pinocchio.updateGlobalPlacements(model,data)
    pinocchio.crba(model, data, q_desired)
    
    pinocchio.computeJointJacobians(model, data, q_desired)
    pinocchio.computeMinverse(model, data, q_desired)

    pinocchio.crba(model, data, q_desired)
    pinocchio.computeCoriolisMatrix(model, data, q_desired, qdot_desired)
    pinocchio.rnea(model, data, q_desired, qdot_z, qddot_z)
    pinocchio.computeMinverse(model,data,q_desired)
    pinocchio.centerOfMass(model,data,False)
    pinocchio.computeCentroidalMomentum(model,data, q_desired, qdot_desired)
    robotCAM = data.hg.angular
    Ag = pinocchio.computeCentroidalMap(model, data, q_desired)
    pinocchio.ccrba(model, data, q_desired, qdot_desired)
    robotIg = data.Ig.inertia
    robotCII = np.linalg.det(np.subtract(np.matmul(robotIg, scipy.linalg.pinv2(robotIginit)), np.identity(3)))

def talker():
    f = open("newfile.txt", 'w')
    f1 = open("newfile1.txt", 'w')
    f2 = open("newfile2.txt", 'w')
    
    modelInitialize()
    for i in range(0,30):
        for j in range(0,30):
            q_command = q_c[:,30*i+j]
            modelUpdate(q_command,qdot_command,qddot_command)
            f.write('%f %f' % (q_c[8,30*i+j],q_c[9,30*i+j]))
            f.write("\n")
            f1.write('%f %f' % (30*i+j, robotCII))
            f1.write("\n")

    f.close()
    f1.close()
    f2.close()

if __name__=='__main__':
    talker()