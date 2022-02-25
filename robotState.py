#!/usr/bin/env python 
from __future__ import print_function
import pinocchio
import sys
import numpy as np
from sys import argv
from os.path import dirname, join, abspath


def talker():
    np.set_printoptions(threshold=sys.maxsize)
    model = pinocchio.buildModelFromUrdf("/home/jhk/legAnalysis/dyros_tocabi_sim.urdf",pinocchio.JointModelFreeFlyer())      
    data = model.createData()
    q = pinocchio.randomConfiguration(model)
    q_init = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.55, 1.26, -0.71, 0.0, 0.0, 0.0, -0.55, 1.26, -0.71, 0.0, 0.0, 0.0, 0.0, 0.2, 0.6, 1.5, -1.47, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -0.2, -0.6, -1.5, 1.47, 1.0, 0.0, 1.0, 0.0]

    for i in range(0, len(q)):
        q[i] = q_init[i]
    qdot = pinocchio.utils.zero(model.nv)
    qddot = pinocchio.utils.zero(model.nv)
    qdot_z = pinocchio.utils.zero(model.nv)
    qddot_z = pinocchio.utils.zero(model.nv)
    pinocchio.forwardKinematics(model, data, q, qdot, qddot)
    pinocchio.updateFramePlacements(model,data)
    pinocchio.computeJointJacobians(model, data, q)

    LF_tran = data.oMi[7].translation
    RF_tran = data.oMi[13].translation
    LF_rot = data.oMi[7].rotation
    RF_rot = data.oMi[13].rotation

    PELV_tran = data.oMi[1].translation
    PELV_rot = data.oMi[1].rotation


    print(RF_tran)
    print(PELV_tran)

    pinocchio.crba(model, data, q)
    pinocchio.computeCoriolisMatrix(model, data, q, qdot)
    pinocchio.rnea(model, data, q, qdot_z, qddot_z)

    frame_id = model.getFrameId("L_Foot_Link")
    LF_j = pinocchio.computeFrameJacobian(model,data,q,frame_id,pinocchio.LOCAL_WORLD_ALIGNED)
    
    frame_id = model.getFrameId("R_Foot_Link")
    RF_j = pinocchio.computeFrameJacobian(model,data,q,frame_id,pinocchio.LOCAL_WORLD_ALIGNED)

    contactState = 1

    #contactPoint[0] = [0.03, 0.0, -0.0135]
    #contactPoint[1] = [0.03, 0.0, -0.0135]
    #contactPoint[2] = [0.03, 0.0, -0.035]
    #contactPoint[3] = [0.03, 0.0, -0.035]
    

    #COR_slice = data.C[:6:39]
    #M_slice = data.M[:6:39]
    #g_slice = data.tau[6:39]
    

    M = data.M
    COR = data.C
    G = data.tau

    #for i in range(0,6):
     #   for j in range(0,39):
      #      print(data.J[i][j], end=' ')
       # print("\n")    

  


if __name__=='__main__':
    talker()