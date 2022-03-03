#!/usr/bin/env python 
from __future__ import print_function
import pinocchio
import sys
import numpy as np
from sys import argv
from os.path import dirname, join, abspath


def talker():
    np.set_printoptions(threshold=sys.maxsize)
    model = pinocchio.buildModelFromUrdf("/home/jhk/legAnalysis/dyros_tocabi_with_redhands.urdf",pinocchio.JointModelFreeFlyer())      
   
    LFframe_id = model.getFrameId("L_Foot_Link")
    RFframe_id = model.getFrameId("R_Foot_Link")

    contactPointLF = pinocchio.SE3.Identity()
    contactPointRF = pinocchio.SE3.Identity()
    
    contactPointLF.translation.T.flat = [0.03, 0, -0.1585]
    contactPointRF.translation.T.flat = [0.03, 0, -0.1585]

    joint_id = model.getJointId("L_AnkleRoll_Joint")
    model.addBodyFrame("LF_contact", joint_id, contactPointLF, LFframe_id)

    joint_id = model.getJointId("R_AnkleRoll_Joint")
    model.addBodyFrame("RF_contact", joint_id, contactPointRF, RFframe_id)
 
    LFcframe_id = model.getFrameId("LF_contact")
    RFcframe_id = model.getFrameId("RF_contact")

    data = model.createData()
    q = pinocchio.randomConfiguration(model)
    q_init = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.55, 1.26, -0.71, 0.0, 0.0, 0.0, -0.55, 1.26, -0.71, 0.0, 0.0, 0.0, 0.0, 0.2, 0.6, 1.5, -1.47, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -0.2, -0.6, -1.5, 1.47, 1.0, 0.0, 1.0, 0.0]

    for i in range(0, len(q)):
        q[i] = q_init[i]

    modeldof = model.nq - 7

    qdot = pinocchio.utils.zero(model.nv)
    qddot = pinocchio.utils.zero(model.nv)
    qdot_z = pinocchio.utils.zero(model.nv)
    qddot_z = pinocchio.utils.zero(model.nv)
    pinocchio.forwardKinematics(model, data, q, qdot, qddot)
    pinocchio.updateFramePlacements(model,data)
    pinocchio.computeJointJacobians(model, data, q)
    pinocchio.computeMinverse(model, data, q)

    LF_tran = data.oMi[7].translation
    RF_tran = data.oMi[13].translation
    LF_rot = data.oMi[7].rotation
    RF_rot = data.oMi[13].rotation

    PELV_tran = data.oMi[1].translation
    PELV_rot = data.oMi[1].rotation

    pinocchio.crba(model, data, q)
    pinocchio.computeCoriolisMatrix(model, data, q, qdot)
    pinocchio.rnea(model, data, q, qdot_z, qddot_z)


    contactState = 1
    contactnum = 0

    if contactState == 1:
        print("double support")
        contactnum = 2
        robotJac = np.zeros((2 * 6, model.nv))
        robotIc = np.zeros((2 * 6, 2 * 6))
    elif contactState == 2:
        print("RF support")
        contactnum = 1
        robotJac = np.zeros((1 * 6, model.nv))
        robotIc = np.zeros((1 * 6, 1 * 6))
    else:
        print("LF support")    
        contactnum = 1
        robotJac = np.zeros((1 * 6, model.nv))
        robotIc = np.zeros((1 * 6, 1 * 6))

    #COR_slice = data.C[:6:39]
    #M_slice = data.M[:6:39]
    #g_slice = data.tau[6:39]


    M = data.M
    COR = data.C
    G = data.tau
    Minv = data.Minv

    LF_j = pinocchio.computeFrameJacobian(model,data,q,LFframe_id,pinocchio.LOCAL_WORLD_ALIGNED)    
    RF_j = pinocchio.computeFrameJacobian(model,data,q,RFframe_id,pinocchio.LOCAL_WORLD_ALIGNED)
    LF_cj = pinocchio.computeFrameJacobian(model,data,q,LFcframe_id,pinocchio.LOCAL_WORLD_ALIGNED)    
    RF_cj = pinocchio.computeFrameJacobian(model,data,q,RFcframe_id,pinocchio.LOCAL_WORLD_ALIGNED)

    for i in range(0, contactnum):
        if i == 0:
            robotJac[0:6,0:model.nv] = LF_cj
        elif i == 1:
            robotJac[6:12,0:model.nv] = RF_cj
        print(i)

    robotLambdac = np.linalg.inv(np.matmul(np.matmul(robotJac,Minv),np.transpose(robotJac)))
    robotJcinvT = np.matmul(np.matmul(robotLambdac, robotJac),Minv)
    robotNc = np.subtract(np.identity(model.nv),np.matmul(np.transpose(robotJac),robotJcinvT))
    robotPc = np.matmul(robotJcinvT,G)
    robotW = np.matmul(Minv[6:6+modeldof,0:model.nq],robotNc[0:model.nq,6:6+modeldof])
    robotWinv = np.linalg.inv(robotW)

    robotContactForce = pinocchio.utils.zero(12)

    if contactState == 1:
        print("ss")
        #robotContactForce = robotJcinvT[0:12,6:6+modeldof]*robotTorque - robotPc
    elif contactState == 2:
        print("RF support")
        #robotContactForce[6:12] = robotJcinvT[0:6,6:6+modeldof]*robotTorque - robotPc
    else:
        print("LF support")    
        #robotContactForce[0:6] = robotJcinvT[0:6,6:6+modeldof]*robotTorque - robotPc
    
    print(robotNc.shape)
    
    #  joint_id = model.getFrameId("L_Foot_Joint")
  #  model.addJoint(joint_id, pinocchio.JointModel, pinocchio.SE3.Identity(),"L_Foot_Joint1")
  
  #  print("\nJoint placements:")
   # for name, oMi in zip(model.names, data.oMi):
    #    print(("{:<24} : {: .2f} {: .2f} {: .2f}"
     #       .format( name, *oMi.translation.T.flat )))

   # print("\nframe placements:")
   # for name, oMf in zip(model.frames, data.oMf):
    #    print(("{:<24} : {: .2f} {: .2f} {: .2f}"
     #       .format( name, *oMf.translation.T.flat )))

    #for i in range(0,6):
     #   for j in range(0,39):
      #      print(data.J[i][j], end=' ')
       # print("\n")    

  
if __name__=='__main__':
    talker()