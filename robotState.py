#!/usr/bin/env python 
from __future__ import print_function
import pinocchio
import sys
import numpy as np
from sys import argv
from os.path import dirname, join, abspath

def walkingSetup():
    global x_direction, step_length, hz, total_tick, foot_step_number, final_step_length, t_total_t, t_temp_t, t_double, t_start, t_total, t_temp, t_last, t_double_1, t_double_2
    global zc, wn, current_step_num, ref_zmp, walking_tick, total_tick
    hz = 1000
    x_direction = 1.0
    step_length = 0.1
    foot_step_number = int(round(x_direction / step_length))
    final_step_length = x_direction % step_length

    if x_direction % step_length != 0:
        foot_step_number = foot_step_number + 1   
    
    t_total_t = 1.2
    t_temp_t = 1.0
    t_double = 0.1
    total_tick = 80000

    t_total = t_total_t * hz
    t_temp = t_temp_t * hz
    t_start = t_temp + 1
    t_last = t_temp + t_total
    t_double_1 = 0.1 * hz
    t_double_2 = 0.1 * hz

    current_step_num = 0
    zc = 0.727822
    wn = np.sqrt(9.81/zc)

    ref_zmp = np.zeros((total_tick,2));
    walking_tick = np.zeros((1,total_tick));

def footStep(): 
    global foot_step
    foot_step = np.zeros((foot_step_number + 1, 2))

    foot_step[0,0] = LF_tran[0]
    foot_step[0,1] = LF_tran[1]

    for i in range(1,foot_step_number + 1):
        foot_step[i,0] = step_length + foot_step[i-1,0]
        foot_step[i,1] = foot_step[i-1,1] * (-1)

    if final_step_length == 0:
        foot_step[foot_step_number,0] = step_length + foot_step[foot_step_number - 1,0]
        foot_step[foot_step_number,1] = foot_step[foot_step_number - 1,1] * (-1)
    else:
        foot_step[foot_step_number,0] = final_step_length + foot_step[foot_step_number - 1,0]
        foot_step[foot_step_number,1] = foot_step[foot_step_number - 1,1] * (-1)

#def zmpGenerator():


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
    q_init = [0.0, 0.0, 0.814175, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.56, 1.3, -0.73, 0.0, 0.0, 0.0, -0.56, 1.3, -0.73, 0.0, 0.0, 0.0, 0.0, 0.2, 0.6, 1.5, -1.47, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -0.2, -0.6, -1.5, 1.47, 1.0, 0.0, 1.0, 0.0]

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
    
    global LF_tran
    global RF_tran

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
        robotdJac = np.zeros((2 * 6, model.nv))
        robotIc = np.zeros((2 * 6, 2 * 6))
    elif contactState == 2:
        print("RF support")
        contactnum = 1
        robotJac = np.zeros((1 * 6, model.nv))
        robotdJac = np.zeros((1 * 6, model.nv))
        robotIc = np.zeros((1 * 6, 1 * 6))
    else:
        print("LF support")    
        contactnum = 1
        robotJac = np.zeros((1 * 6, model.nv))
        robotdJac = np.zeros((1 * 6, model.nv))
        robotIc = np.zeros((1 * 6, 1 * 6))

    M = data.M
    COR = data.C
    G = data.tau
    Minv = data.Minv
    b = np.matmul(COR,qdot)

    LF_j = pinocchio.computeFrameJacobian(model,data,q,LFframe_id,pinocchio.LOCAL_WORLD_ALIGNED)    
    RF_j = pinocchio.computeFrameJacobian(model,data,q,RFframe_id,pinocchio.LOCAL_WORLD_ALIGNED)
    LF_cj = pinocchio.computeFrameJacobian(model,data,q,LFcframe_id,pinocchio.LOCAL_WORLD_ALIGNED)    
    RF_cj = pinocchio.computeFrameJacobian(model,data,q,RFcframe_id,pinocchio.LOCAL_WORLD_ALIGNED)
    RF_cdj = pinocchio.frameJacobianTimeVariation(model,data,q,qdot,RFcframe_id,pinocchio.LOCAL_WORLD_ALIGNED)
    LF_cdj = pinocchio.frameJacobianTimeVariation(model,data,q,qdot,LFcframe_id,pinocchio.LOCAL_WORLD_ALIGNED)

    for i in range(0, contactnum):
        if i == 0:
            robotJac[0:6,0:model.nv] = LF_cj
            robotdJac[0:6,0:model.nv] = LF_cdj
        elif i == 1:
            robotJac[6:12,0:model.nv] = RF_cj
            robotdJac[6:12,0:model.nv] = RF_cdj
        print(i)

    robotLambdac = np.linalg.inv(np.matmul(np.matmul(robotJac,Minv),np.transpose(robotJac)))
    robotJcinvT = np.matmul(np.matmul(robotLambdac, robotJac),Minv)
    robotNc = np.subtract(np.identity(model.nv),np.matmul(np.transpose(robotJac),robotJcinvT))
    robotPc = np.matmul(robotJcinvT,G)
   # robotW = np.matmul(Minv[6:6+modeldof,0:model.nq],robotNc[0:model.nq,6:6+modeldof])
   # robotWinv = np.linalg.inv(robotW)

    robotmuc = np.matmul(robotLambdac,np.subtract(np.matmul(np.matmul(robotJac,Minv),b),np.matmul(robotdJac,qdot)))
    robothc = np.matmul(np.transpose(robotJac),np.add(robotmuc, robotPc))

    robotContactForce = pinocchio.utils.zero(12)

    robotTorque = np.matmul(np.linalg.pinv(robotNc),np.subtract(np.add(np.add(np.matmul(M,qddot),b),G),robothc))

    if contactState == 1:
        robotContactForce = np.subtract(np.subtract(np.matmul(robotJcinvT,robotTorque), robotPc),robotmuc)
    elif contactState == 2:
        robotContactForce[6:12] = np.subtract(np.subtract(np.matmul(robotJcinvT,robotTorque), robotPc),robotmuc)
    else:   
        robotContactForce[0:6] = np.subtract(np.subtract(np.matmul(robotJcinvT,robotTorque), robotPc),robotmuc)

    walkingSetup()
    footStep()

    print(robotTorque)
    print(robotContactForce)
  
if __name__=='__main__':
    talker()