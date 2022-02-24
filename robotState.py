#!/usr/bin/env python 
import pinocchio
import numpy as np
from sys import argv
from os.path import dirname, join, abspath


def talker():
    model = pinocchio.buildModelFromUrdf("/home/jhk/legAnalysis/dyros_tocabi_sim.urdf",pinocchio.JointModelFreeFlyer())      
    data = model.createData()
    q = pinocchio.randomConfiguration(model)
    q_init = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.55, 1.26, -0.71, 0.0, 0.0, 0.0, -0.55, 1.26, -0.71, 0.0, 0.0, 0.0, 0.0, 0.2, 0.6, 1.5, -1.47, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -0.2, -0.6, -1.5, 1.47, 1.0, 0.0, 1.0, 0.0]

    for i in range(0, len(q)):
        q[i] = q_init[i]
    qdot = pinocchio.utils.zero(model.nv)
    qddot = pinocchio.utils.zero(model.nv)
    print(q)
    print(qdot)
    print(qddot)
    pinocchio.forwardKinematics(model, data, q, qdot, qddot)
    pinocchio.computeJointJacobians(model, data, q)

    LF_tran = data.oMi[7].translation
    RF_tran = data.oMi[13].translation
    LF_rot = data.oMi[7].rotation
    RF_rot = data.oMi[13].rotation

    PELV_tran = data.oMi[1].translation
    PELV_rot = data.oMi[1].rotation


   # pinocchio.crba(model data q)
   # pinocchio.computeCoriolisMatrix(model data q qdot)
   # pinocchio.rnea(model data q qdot qddot)
    

    COR_slice = data.C[:6:39]
    M_slice = data.M[:6:39]
    g_slice = data.tau[6:39]

  


if __name__=='__main__':
    talker()