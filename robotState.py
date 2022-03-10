#!/usr/bin/env python 
from __future__ import print_function
import pinocchio
import sys
import numpy as np
from sys import argv
from os.path import dirname, join, abspath

def rotateWithY(pitch_angle):
    rotate_with_y = np.zeros((3,3))

    rotate_with_y[0, 0] = np.cos(pitch_angle)
    rotate_with_y[1, 0] = 0.0
    rotate_with_y[2, 0] = -1 * np.sin(pitch_angle)

    rotate_with_y[0, 1] = 0.0
    rotate_with_y[1, 1] = 1.0
    rotate_with_y[2, 1] = 0.0

    rotate_with_y[0, 2] = np.sin(pitch_angle)
    rotate_with_y[1, 2] = 0.0
    rotate_with_y[2, 2] = np.cos(pitch_angle)

    return rotate_with_y 

def rotateWithX(roll_angle):
    rotate_with_x = np.zeros((3,3))

    rotate_with_x[0, 0] = 1.0
    rotate_with_x[1, 0] = 0.0
    rotate_with_x[2, 0] = 0.0

    rotate_with_x[0, 1] = 0.0
    rotate_with_x[1, 1] = np.cos(roll_angle)
    rotate_with_x[2, 1] = np.sin(roll_angle)

    rotate_with_x[0, 2] = 0.0
    rotate_with_x[1, 2] = -1 * np.sin(roll_angle)
    rotate_with_x[2, 2] = np.cos(roll_angle)

    return rotate_with_x     

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

    ref_zmp = np.zeros((total_tick,2))
    walking_tick = np.zeros((1,total_tick))

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

def zmpGenerator():
    print("ZMP")

def comGenerator():
    print("COM")

def inverseKinematics(LF_rot_c, RF_rot_c, PELV_rot_c, LF_tran_c, RF_tran_c, PELV_tran_c, HRR_tran_init_c, HLR_tran_init_c, HRR_rot_init_c, HLR_rot_init_c, PELV_tran_init_c, PELV_rot_init_c, CPELV_tran_init_c):
    global leg_q
    M_PI = 3.14159265358979323846
    leg_q = np.zeros(12)

    l_upper = 0.35
    l_lower = 0.254

    offset_hip_pitch = 0.0
    offset_knee_pitch = 0.0
    offset_ankle_pitch = 0.0

    lpt = np.subtract(PELV_tran_c, LF_tran_c)
    rpt = np.subtract(PELV_tran_c, RF_tran_c)
    lp = np.matmul(np.transpose(LF_rot_c), np.transpose(lpt))
    rp = np.matmul(np.transpose(RF_rot_c), np.transpose(rpt))
    
    PELF_rot = np.matmul(np.transpose(PELV_rot_c), np.transpose(LF_rot_c))
    PERF_rot = np.matmul(np.transpose(PELV_rot_c), np.transpose(LF_rot_c))

    ld = np.zeros(3)  
    rd = np.zeros(3)

    ld[0] = HLR_tran_init_c[0] - PELV_tran_init_c[0]
    ld[1] = HLR_tran_init_c[1] - PELV_tran_init_c[1]
    ld[2] = -(CPELV_tran_init_c[2] - HLR_tran_init_c[2]) + (CPELV_tran_init_c[2] - PELV_tran_init_c[2])

    rd[0] = HRR_tran_init_c[0] - PELV_tran_init_c[0]
    rd[1] = HRR_tran_init_c[1] - PELV_tran_init_c[1]
    rd[2] = -(CPELV_tran_init_c[2] - HRR_tran_init_c[2]) + (CPELV_tran_init_c[2] - PELV_tran_init_c[2])

    ld = np.matmul(np.transpose(LF_rot_c), ld)
    rd = np.matmul(np.transpose(RF_rot_c), rd)

    lr = lp + ld
    rr = rp + rd

    lc = np.linalg.norm(lr)*np.linalg.norm(lr)

    leg_q[3] = -1 * np.arccos((l_upper * l_upper + l_lower * l_lower - lc * lc) / (2 * l_upper * l_lower)) + M_PI
    l_ankle_pitch = np.arcsin((l_upper * np.sin(M_PI - leg_q[3])) / lc)
    
    leg_q[4] = -1 * np.arctan2(lr[0], np.sqrt(lr[1] * lr[1] + lr[2] * lr[2])) - l_ankle_pitch
    leg_q[5] = np.arctan2(lr[1], lr[2])

    r_tl2 = np.zeros((3,3))
    r_l2l3 = np.zeros((3,3))
    r_l3l4 = np.zeros((3,3))
    r_l4l5 = np.zeros((3,3))

    r_l2l3 = rotateWithY(leg_q[3])
    r_l3l4 = rotateWithY(leg_q[4])
    r_l4l5 = rotateWithX(leg_q[5])

    r_tl2 = np.matmul(np.matmul(np.matmul(PELV_rot, np.transpose(r_l4l5)),np.transpose(r_l3l4)),np.transpose(r_l2l3))
    leg_q[1] = np.arcsin(r_tl2[2, 1])

    c_lq5 = np.divide(-r_tl2[0, 1], np.cos(leg_q[1]))

    if c_lq5 > 1.0:
        c_lq5 = 1.0
    elif c_lq5 < -1.0:
        c_lg5 = -1.0

    leg_q[0] = -1 * np.arcsin(c_lq5)
    leg_q[2] = -1 * np.arcsin(r_tl2[2, 0] / np.cos(leg_q[1])) + offset_hip_pitch
    leg_q[3] = leg_q[3] - offset_knee_pitch
    leg_q[4] = leg_q[4] - offset_ankle_pitch

    rc = np.linalg.norm(rr)*np.linalg.norm(rr)
    leg_q[9] = -1 * np.arccos((l_upper * l_upper + l_lower * l_lower - rc * rc) / (2 * l_upper * l_lower)) + M_PI

    r_ankle_pitch = np.arcsin((l_upper * np.sin(M_PI - leg_q[9])) / rc)
    leg_q[10] = -1 * np.arctan2(rr[0], np.sqrt(rr[1] * rr[1] + rr[2] * rr[2])) - r_ankle_pitch
    leg_q[11] = np.arctan2(rr[1], rr[2])
    r_tr2 = np.zeros((3,3))
    r_r2r3 = np.zeros((3,3))
    r_r3r4 = np.zeros((3,3))
    r_r4r5 = np.zeros((3,3))

    r_r2r3 = rotateWithY(leg_q[0])
    r_r3r4 = rotateWithY(leg_q[10])
    r_r4r5 = rotateWithX(leg_q[11])

    r_tl2 = np.matmul(np.matmul(np.matmul(PELV_rot, np.transpose(r_r4r5)),np.transpose(r_r3r4)),np.transpose(r_r2r3))
    leg_q[7] = np.arcsin(r_tr2[2,1])
    c_rq5 = -r_tr2[0, 1] / np.cos(leg_q[7])

    if c_rq5 > 1.0:
        c_rq5 = 1.0
    elif c_rq5 < -1.0:
        c_rg5 = -1.0
    
    leg_q[6] = -1* np.arcsin(c_rq5)
    leg_q[8] = np.arcsin(r_tr2[2, 0] / np.cos(leg_q[7])) - offset_hip_pitch
    leg_q[9] = -1 * leg_q[9] + offset_knee_pitch
    leg_q[10] = -1 * leg_q[10] + offset_ankle_pitch

    leg_q[0] = leg_q[0] * (-1)
    leg_q[6] = leg_q[6] * (-1)
    leg_q[8] = leg_q[8] * (-1)
    leg_q[9] = leg_q[9] * (-1)
    leg_q[10] = leg_q[10] * (-1)

    print(leg_q)
    print("IK")

def modelInitialize():
    global model, data, LFframe_id, RFframe_id, LFcframe_id, RFcframe_id, q, qdot, qddot, LF_tran, RF_tran, PELV_tran, LF_rot, RF_rot, PELV_rot, qdot_z, qddot_z, HRR_rot_init, HLR_rot_init, HRR_tran_init, HLR_tran_init, PELV_tran_init, PELV_rot_init, CPELV_tran_init
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
    pinocchio.updateGlobalPlacements(model,data)
    pinocchio.computeJointJacobians(model, data, q)
    pinocchio.computeMinverse(model, data, q)

    LF_tran = data.oMi[7].translation
    RF_tran = data.oMi[13].translation
    LF_rot = data.oMi[7].rotation
    RF_rot = data.oMi[13].rotation

    PELV_tran = data.oMi[1].translation
    PELV_rot = data.oMi[1].rotation

    LF_tran_init = data.oMi[7].translation
    RF_tran_init = data.oMi[13].translation
    HLR_tran_init = data.oMi[4].translation
    HRR_tran_init = data.oMi[10].translation
    LF_rot_init = data.oMi[7].rotation
    RF_rot_init = data.oMi[13].rotation
    HLR_rot_init = data.oMi[4].rotation
    HRR_rot_init = data.oMi[10].rotation

    PELV_tran_init = data.oMi[1].translation
    CPELV_tran_init = data.oMi[1].translation
    PELV_rot_init = data.oMi[1].rotation

    print("C")
    print(CPELV_tran_init)
    print(PELV_tran_init)
    print(RF_tran_init)
    print(LF_tran_init)
    print(LFframe_id)

  
  #  CPELV_tran_init[0] = -0.00739
  #  CPELV_tran_init[1] = 0.0
  #  CPELV_tran_init[2] = 0.849

  #  PELV_tran_init[0] = 0.0695
   # PELV_tran_init[1] = 0.0
   # PELV_tran_init[2] = 0.849

   # PELV_tran[0] = 0.0695
   # PELV_tran[1] = 0.0
   # PELV_tran[2] = 0.849

   # LF_tran_init[0] = -0.017
   # LF_tran_init[1] = 0.102
   # LF_tran_init[2] = -0.696

    #RF_tran_init[0] = -0.017
    #RF_tran_init[1] = -0.102
    #RF_tran_init[2] = -0.696

    #LF_tran[0] = -0.017
    #LF_tran[1] = 0.102
    #LF_tran[2] = -0.696

    #RF_tran[0] = -0.017
    #RF_tran[1] = -0.102
    #RF_tran[2] = -0.696

    #HLR_tran_init[0] = 0.105
    #HLR_tran_init[1] = 0.1025
    #HLR_tran_init[2] = 0.707

    #HRR_tran_init[0] = 0.105
    #HRR_tran_init[1] = -0.1025
    #HRR_tran_init[2] = 0.707


def modelUpdate():
    global contactState, contactnum, M, G, COR, Minv, b, robotJac, robotdJac, robotIc, LF_j, RF_j, LF_cj, RF_cj, LF_cdj, RF_cdj, robotLambdac, robotJcinvT, robotNc, robotPc, robotmuc, robothc
    pinocchio.forwardKinematics(model, data, q, qdot, qddot)
    pinocchio.updateFramePlacements(model,data)
    pinocchio.updateGlobalPlacements(model,data)
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

    robotLambdac = np.linalg.inv(np.matmul(np.matmul(robotJac,Minv),np.transpose(robotJac)))
    robotJcinvT = np.matmul(np.matmul(robotLambdac, robotJac),Minv)
    robotNc = np.subtract(np.identity(model.nv),np.matmul(np.transpose(robotJac),robotJcinvT))
    robotPc = np.matmul(robotJcinvT,G)

    robotmuc = np.matmul(robotLambdac,np.subtract(np.matmul(np.matmul(robotJac,Minv),b),np.matmul(robotdJac,qdot)))
    robothc = np.matmul(np.transpose(robotJac),np.add(robotmuc, robotPc))

def estimateContactForce():
    robotContactForce = pinocchio.utils.zero(12)

    robotTorque = np.matmul(np.linalg.pinv(robotNc),np.subtract(np.add(np.add(np.matmul(M,qddot),b),G),robothc))

    if contactState == 1:
        robotContactForce = np.subtract(np.subtract(np.matmul(robotJcinvT,robotTorque), robotPc),robotmuc)
    elif contactState == 2:
        robotContactForce[6:12] = np.subtract(np.subtract(np.matmul(robotJcinvT,robotTorque), robotPc),robotmuc)
    else:   
        robotContactForce[0:6] = np.subtract(np.subtract(np.matmul(robotJcinvT,robotTorque), robotPc),robotmuc)

    print(robotTorque)
    print(robotContactForce)

def talker():
    modelInitialize()
   # walkingSetup()
    #footStep()
    #zmpGenerator()
    #comGenerator()

    inverseKinematics(LF_rot, RF_rot, PELV_rot, LF_tran, RF_tran, PELV_tran, HRR_tran_init, HLR_tran_init, HRR_rot_init, HLR_rot_init, PELV_tran_init, PELV_rot_init, CPELV_tran_init)
  #  modelUpdate()
   # estimateContactForce()
  
if __name__=='__main__':
    talker()