#!/usr/bin/env python 
from __future__ import print_function
from tempfile import tempdir
import pinocchio
from pinocchio.utils import npToTuple
from pinocchio.rpy import matrixToRpy, rpyToMatrix
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.linalg
from sys import argv
from os.path import dirname, join, abspath

np.set_printoptions(threshold=sys.maxsize)

def quinticSpline(time, time_0, time_f, x_0, x_dot_0, x_ddot_0, x_f, x_dot_f, x_ddot_f):
    time_s = time_f - time_0
    a1 = x_0
    a2 = x_dot_0
    a3 = x_ddot_0 / 2.0

    Temp = np.zeros((3,3))
    R_temp = np.zeros(3)

    Temp[0,0] = math.pow(time_s, 3)
    Temp[0,1] = math.pow(time_s, 4)
    Temp[0,2] = math.pow(time_s, 5)
    Temp[1,0] = 3.0 * math.pow(time_s, 2)
    Temp[1,1] = 4.0 * math.pow(time_s, 3)
    Temp[1,2] = 5.0 * math.pow(time_s, 4)
    Temp[2,0] = 6.0 * time_s
    Temp[2,1] = 12.0 * math.pow(time_s, 2)
    Temp[2,2] = 20.0 * math.pow(time_s, 3)

    R_temp[0] = x_f - x_0 - x_dot_0 * time_s - x_ddot_0 * math.pow(time_s, 2) / 2.0
    R_temp[1] = x_dot_f - x_dot_0 - x_ddot_0 * time_s
    R_temp[2] = x_ddot_f - x_ddot_0

    RES = np.matmul(np.linalg.inv(Temp), R_temp)

    a4 = RES[0]
    a5 = RES[1]
    a6 = RES[2]

    time_fs = time - time_0

    position = a1 + a2 * math.pow(time_fs, 1) + a3 * math.pow(time_fs, 2) + a4 * math.pow(time_fs, 3) + a5 * math.pow(time_fs, 4) + a6 * math.pow(time_fs, 5)
    
    result = position

    if time < time_0:
      result = x_0
    elif time > time_f:
      result = x_f

    return result
def quinticSplineDot(time, time_0, time_f, x_0, x_dot_0, x_ddot_0, x_f, x_dot_f, x_ddot_f):
    time_s = time_f - time_0
    a1 = x_0
    a2 = x_dot_0
    a3 = x_ddot_0 / 2.0

    Temp = np.zeros((3,3))
    R_temp = np.zeros(3)

    Temp[0,0] = math.pow(time_s, 3)
    Temp[0,1] = math.pow(time_s, 4)
    Temp[0,2] = math.pow(time_s, 5)
    Temp[1,0] = 3.0 * math.pow(time_s, 2)
    Temp[1,1] = 4.0 * math.pow(time_s, 3)
    Temp[1,2] = 5.0 * math.pow(time_s, 4)
    Temp[2,0] = 6.0 * time_s
    Temp[2,1] = 12.0 * math.pow(time_s, 2)
    Temp[2,2] = 20.0 * math.pow(time_s, 3)

    R_temp[0] = x_f - x_0 - x_dot_0 * time_s - x_ddot_0 * math.pow(time_s, 2) / 2.0
    R_temp[1] = x_dot_f - x_dot_0 - x_ddot_0 * time_s
    R_temp[2] = x_ddot_f - x_ddot_0

    RES = np.matmul(np.linalg.inv(Temp), R_temp)

    a4 = RES[0]
    a5 = RES[1]
    a6 = RES[2]

    time_fs = time - time_0

    position = a1 + a2 * math.pow(time_fs, 1) + a3 * math.pow(time_fs, 2) + a4 * math.pow(time_fs, 3) + a5 * math.pow(time_fs, 4) + a6 * math.pow(time_fs, 5)
    velocity = a2 + 2.0 * a3 * math.pow(time_fs, 1) + 3.0 * a4 * math.pow(time_fs, 2) + 4.0 * a5 * math.pow(time_fs, 3) + 5.0 * a6 * math.pow(time_fs, 4);
    
    result = velocity

    if time < time_0:
      result = x_dot_0
    elif time > time_f:
      result = x_dot_f

    return result

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

def rotateWithZ(yaw_angle):
    rotate_with_z = np.zeros((3,3))

    rotate_with_z[0, 0] = np.cos(yaw_angle)
    rotate_with_z[1, 0] = np.sin(yaw_angle)
    rotate_with_z[2, 0] = 0.0

    rotate_with_z[0, 1] = -1 * np.sin(yaw_angle)
    rotate_with_z[1, 1] = np.cos(yaw_angle)
    rotate_with_z[2, 1] = 0.0

    rotate_with_z[0, 2] = 0.0
    rotate_with_z[1, 2] = 0.0
    rotate_with_z[2, 2] = 1.0

    return rotate_with_z

def skew(vector):
    return np.array([[0, -vector[2], vector[1]], 
                     [vector[2], 0, -vector[0]], 
                     [-vector[1], vector[0], 0]])   

def winvCalc(W):
    global V2, thres
    thres = 0
    u, s, v = np.linalg.svd(W, full_matrices=False)
    diag_s = np.diag(s)
 
    for i in range(0, (np.shape(s))[0]):
        if(diag_s[i, i] < 0.00001):
            thres = i
            break

    #print(thres)
    V2 = v[thres:(np.shape(s))[0], 0:(np.shape(s))[0]]
 
def walkingSetup():
    global x_direction, y_direction, yaw_direction, step_length, hz, total_tick, t_total_t, t_start_real, t_temp_t, t_double, t_rest_1, t_rest_2, t_start, t_total, t_temp, t_last, t_double_1, t_double_2
    global zc, wn, current_step_num, ref_zmp, ref_com, time, total_tick, phase_variable, lfoot, rfoot, foot_height, foot_step_dir, A_lipm, B_lipm, C_lipm, K_lipm, f_lipm, Ks_lipm, Kx_lipm, G_lipm, N_preview
    hz = 500
    x_direction = 0.5
    y_direction = 0.00
    yaw_direction = 0.00
    step_length = 0.1
    
    velocity_control = 1.0
    t_total_t = 1.0 * velocity_control
    t_temp_t = 2.0
    t_double = 0.2 * velocity_control
    
    t_total = t_total_t * hz
    t_temp = t_temp_t * hz
    t_start = t_temp + 1
    t_last = t_temp + t_total
    t_double_1 = 0.1 * hz * velocity_control
    t_double_2 = 0.1 * hz * velocity_control
    t_rest_1 = 0.1 * hz * velocity_control
    t_rest_2 = 0.1 * hz * velocity_control
    t_start_real = t_start + t_rest_1

    foot_step_dir = 1

    foot_height = 0.05

    current_step_num = 0
    zc = 0.727822
    wn = np.sqrt(9.81/zc)

    dt = 1/float(hz)

    A_lipm = np.mat(([1, dt, 0],
                [0, 1, dt],
                [0, 0, 1]))
    B_lipm = np.mat((0, 0, dt)).T
    C_lipm = np.mat((1, 0, -zc/9.81))

    Q_lipm = 1
    R_lipm = 1e-6

    t_preview = 1.6 * velocity_control
    N_preview = int(t_preview/dt)
    K_lipm, f_lipm = calculatePreviewControlParams(A_lipm, B_lipm, C_lipm, Q_lipm, R_lipm, N_preview)
    Ks_lipm, Kx_lipm, G_lipm = calculatePreviewControlParams2(A_lipm, B_lipm, C_lipm, Q_lipm, R_lipm, N_preview)

def footStep(): 
    global foot_step_number, foot_step
    
    initial_rot = math.atan2(y_direction, x_direction)
    if initial_rot > 0.0:
        initial_drot = np.deg2rad(10.0)
    else:
        initial_drot = np.deg2rad(-10.0)

    init_totalstep_num = int(initial_rot / initial_drot)
    init_residual_angle = initial_rot - init_totalstep_num * initial_drot

    final_rot = yaw_direction - initial_rot
    if final_rot > 0.0:
        final_drot = np.deg2rad(10.0)
    else:
        final_drot = np.deg2rad(-10.0)

    final_foot_step_number = int(final_rot / final_drot)

    final_residual_angle = final_rot - final_foot_step_number * final_drot
    l = np.sqrt(x_direction * x_direction + y_direction * y_direction)
    dlength = step_length
    middle_foot_step_number = int(l / dlength)
    middle_residual_length = l - middle_foot_step_number * dlength
    del_size = 1
    numberOfFootstep = init_totalstep_num * del_size + middle_foot_step_number * del_size + final_foot_step_number * del_size
    
    if init_totalstep_num != 0 or np.abs(init_residual_angle) >= 0.0001:
        if init_totalstep_num % 2 == 0:
            numberOfFootstep = numberOfFootstep + 2
        else:
            if (np.abs(init_residual_angle) >= 0.0001):
                numberOfFootstep = numberOfFootstep + 3
            else:
                numberOfFootstep = numberOfFootstep + 2

    if (middle_foot_step_number != 0 or np.abs(middle_residual_length) >= 0.0001):
        if (middle_foot_step_number % 2 == 0):
            numberOfFootstep = numberOfFootstep + 2
        else:
            if (np.abs(middle_residual_length) >= 0.0001):
                numberOfFootstep = numberOfFootstep + 3
            else:
                numberOfFootstep = numberOfFootstep + 2

    if (final_foot_step_number != 0 or abs(final_residual_angle) >= 0.0001):
        if (abs(final_residual_angle) >= 0.0001):
            numberOfFootstep = numberOfFootstep + 2
        else:
            numberOfFootstep = numberOfFootstep + 2

    numberOfFootstep = numberOfFootstep + 1
    foot_step = np.zeros((numberOfFootstep, 7))
    index = 0

    if (foot_step_dir == 1):
        is_right = 1
    else:
        is_right = -1

    temp = -is_right
    temp2 = -is_right
    temp3 = is_right

    if (init_totalstep_num != 0 or np.abs(init_residual_angle) >= 0.0001):
        for i in range(0, init_totalstep_num):
            temp = -1 * temp
            foot_step[index, 0] = temp * foot_distance[1] / 2.0 * np.sin((i + 1) * initial_drot)
            foot_step[index, 1] = -temp * foot_distance[1] / 2.0 * np.cos((i + 1) * initial_drot)
            foot_step[index, 5] = (i + 1) * initial_drot
            foot_step[index, 6] = 0.5 + 0.5 * temp
            index = index + 1

        if (temp == is_right):
            if (abs(init_residual_angle) >= 0.0001):
                temp = -1 * temp

                foot_step[index, 0] = temp * foot_distance[1] / 2.0 * np.sin((init_totalstep_num)*initial_drot + init_residual_angle)
                foot_step[index, 1] = -temp * foot_distance[1] / 2.0 * np.cos((init_totalstep_num)*initial_drot + init_residual_angle)
                foot_step[index, 5] = (init_totalstep_num)*initial_drot + init_residual_angle
                foot_step[index, 6] = 0.5 + 0.5 * temp
                index = index + 1

                temp = -1 * temp

                foot_step[index, 0] = temp * foot_distance[1] / 2.0 * np.sin((init_totalstep_num)*initial_drot + init_residual_angle)
                foot_step[index, 1] = -temp * foot_distance[1] / 2.0 * np.cos((init_totalstep_num)*initial_drot + init_residual_angle)
                foot_step[index, 5] = (init_totalstep_num)*initial_drot + init_residual_angle
                foot_step[index, 6] = 0.5 + 0.5 * temp
                index = index + 1

                temp = -1 * temp

                foot_step[index, 0] = temp * foot_distance[1] / 2.0 * np.sin((init_totalstep_num)*initial_drot + init_residual_angle)
                foot_step[index, 1] = -temp * foot_distance[1] / 2.0 * np.cos((init_totalstep_num)*initial_drot + init_residual_angle)
                foot_step[index, 5] = (init_totalstep_num)*initial_drot + init_residual_angle
                foot_step[index, 6] = 0.5 + 0.5 * temp
                index = index + 1
            else:
                temp = -1 * temp

                foot_step[index, 0] = temp * foot_distance[1] / 2.0 * np.sin((init_totalstep_num)*initial_drot + init_residual_angle)
                foot_step[index, 1] = -temp * foot_distance[1] / 2.0 * np.cos((init_totalstep_num)*initial_drot + init_residual_angle)
                foot_step[index, 5] = (init_totalstep_num)*initial_drot + init_residual_angle
                foot_step[index, 6] = 0.5 + 0.5 * temp
                index = index + 1
        elif (temp == -is_right):
            temp = -1 * temp

            foot_step[index, 0] = temp * foot_distance[1] / 2.0 * np.sin((init_totalstep_num)*initial_drot + init_residual_angle)
            foot_step[index, 1] = -temp * foot_distance[1] / 2.0 * np.cos((init_totalstep_num)*initial_drot + init_residual_angle)
            foot_step[index, 5] = (init_totalstep_num)*initial_drot + init_residual_angle
            foot_step[index, 6] = 0.5 + 0.5 * temp
            index = index + 1

            temp = -1 * temp

            foot_step[index, 0] = temp * foot_distance[1] / 2.0 * np.sin((init_totalstep_num)*initial_drot + init_residual_angle)
            foot_step[index, 1] = -temp * foot_distance[1] / 2.0 * np.cos((init_totalstep_num)*initial_drot + init_residual_angle)
            foot_step[index, 5] = (init_totalstep_num)*initial_drot + init_residual_angle
            foot_step[index, 6] = 0.5 + 0.5 * temp
            index = index + 1

    if (middle_foot_step_number != 0 or abs(middle_residual_length) >= 0.0001):
        temp2 = -1 * temp2

        foot_step[index, 0] = 0.0
        foot_step[index, 1] = -temp2 * (foot_distance[1] / 2.0)
        foot_step[index, 5] = 0.0
        foot_step[index, 6] = 0.5 + 0.5 * temp2

        index = index + 1

        for i in range(0, middle_foot_step_number):
            temp2 = -1 * temp2

            foot_step[index, 0] = np.cos(initial_rot) * (dlength * (i + 1)) + temp2 * np.sin(initial_rot) * (foot_distance[1] / 2.0)
            foot_step[index, 1] = np.sin(initial_rot) * (dlength * (i + 1)) - temp2 * np.cos(initial_rot) * (foot_distance[1] / 2.0)
            foot_step[index, 5] = initial_rot
            foot_step[index, 6] = 0.5 + 0.5 * temp2
            index = index + 1
            
        if (temp2 == -is_right):
            if (np.abs(middle_residual_length) >= 0.0001):
                temp2 = -1 * temp2

                foot_step[index, 0] = np.cos(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length) + temp2 * np.sin(initial_rot) * (foot_distance[1] / 2.0)
                foot_step[index, 1] = np.sin(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length) - temp2 * np.cos(initial_rot) * (foot_distance[1] / 2.0)
                foot_step[index, 5] = initial_rot
                foot_step[index, 6] = 0.5 + 0.5 * temp2
                index = index + 1

                temp2 = -1 * temp2

                foot_step[index, 0] = np.cos(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length) + temp2 * np.sin(initial_rot) * (foot_distance[1] / 2.0)
                foot_step[index, 1] = np.sin(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length) - temp2 * np.cos(initial_rot) * (foot_distance[1] / 2.0)
                foot_step[index, 5] = initial_rot
                foot_step[index, 6] = 0.5 + 0.5 * temp2
                index = index + 1

                temp2 = -1 * temp2

                foot_step[index, 0] = np.cos(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length) + temp2 * np.sin(initial_rot) * (foot_distance[1] / 2.0)
                foot_step[index, 1] = np.sin(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length) - temp2 * np.cos(initial_rot) * (foot_distance[1] / 2.0)
                foot_step[index, 5] = initial_rot
                foot_step[index, 6] = 0.5 + 0.5 * temp2
                index = index + 1
            else:
                temp2 = -1 * temp2
                foot_step[index, 0] = np.cos(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length) + temp2 * np.sin(initial_rot) * (foot_distance[1] / 2.0)
                foot_step[index, 1] = np.sin(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length) - temp2 * np.cos(initial_rot) * (foot_distance[1] / 2.0)
                foot_step[index, 5] = initial_rot
                foot_step[index, 6] = 0.5 + 0.5 * temp2
                index = index + 1
                temp2 = -1 * temp2

                foot_step[index, 0] = np.cos(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length) + temp2 * np.sin(initial_rot) * (foot_distance[1] / 2.0)
                foot_step[index, 1] = np.sin(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length) - temp2 * np.cos(initial_rot) * (foot_distance[1] / 2.0)
                foot_step[index, 5] = initial_rot
                foot_step[index, 6] = 0.5 + 0.5 * temp2
                index = index + 1
        elif (temp2 == is_right):
            temp2 = -1 * temp2
            foot_step[index, 0] = np.cos(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length) + temp2 * np.sin(initial_rot) * (foot_distance[1] / 2.0)
            foot_step[index, 1] = np.sin(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length) - temp2 * np.cos(initial_rot) * (foot_distance[1] / 2.0)
            foot_step[index, 5] = initial_rot
            foot_step[index, 6] = 0.5 + 0.5 * temp2
            index = index + 1

            temp2 = -1 * temp2

            foot_step[index, 0] = np.cos(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length) + temp2 * np.sin(initial_rot) * (foot_distance[1] / 2.0)
            foot_step[index, 1] = np.sin(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length) - temp2 * np.cos(initial_rot) * (foot_distance[1] / 2.0)
            foot_step[index, 5] = initial_rot
            foot_step[index, 6] = 0.5 + 0.5 * temp2
            index = index + 1

    final_position_x = np.cos(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length)
    final_position_y = np.sin(initial_rot) * (dlength * (middle_foot_step_number) + middle_residual_length)

    if (final_foot_step_number != 0 or abs(final_residual_angle) >= 0.0001):
        for i in range(0, final_foot_step_number):
            temp3 = -1 * temp3

            foot_step[index, 0] = final_position_x + temp3 * foot_distance[1] / 2.0 * np.sin((i + 1) * final_drot + initial_rot)
            foot_step[index, 1] = final_position_y - temp3 * foot_distance[1] / 2.0 * np.cos((i + 1) * final_drot + initial_rot)
            foot_step[index, 5] = (i + 1) * final_drot + initial_rot
            foot_step[index, 6] = 0.5 + 0.5 * temp3
            index = index + 1

        if (abs(final_residual_angle) >= 0.0001):
            temp3 = -1 * temp3

            foot_step[index, 0] = final_position_x + temp3 * foot_distance[1] / 2.0 * np.sin(yaw_direction)
            foot_step[index, 1] = final_position_y - temp3 * foot_distance[1] / 2.0 * np.cos(yaw_direction)
            foot_step[index, 5] = yaw_direction
            foot_step[index, 6] = 0.5 + 0.5 * temp3
            index = index + 1

            temp3 = -1 * temp3

            foot_step[index, 0] = final_position_x + temp3 * foot_distance[1] / 2.0 * np.sin(yaw_direction)
            foot_step[index, 1] = final_position_y - temp3 * foot_distance[1] / 2.0 * np.cos(yaw_direction)
            foot_step[index, 5] = yaw_direction
            foot_step[index, 6] = 0.5 + 0.5 * temp3
            index = index + 1
        else:
            temp3 = -1 * temp3

            foot_step[index, 0] = final_position_x + temp3 * foot_distance[1] / 2.0 * np.sin(yaw_direction)
            foot_step[index, 1] = final_position_y - temp3 * foot_distance[1] / 2.0 * np.cos(yaw_direction)
            foot_step[index, 5] = yaw_direction
            foot_step[index, 6] = 0.5 + 0.5 * temp3
            index = index + 1

    for i in range(0, numberOfFootstep):
        if (foot_step[i, 6] == 1):
            foot_step[i, 0] = foot_step[i, 0] + RF_tran[0]
            foot_step[i, 1] = RF_tran[1]
        else:
            foot_step[i, 0] = foot_step[i, 0] + LF_tran[0]
            foot_step[i, 1] = LF_tran[1]
    
    foot_step_number = numberOfFootstep 

def cpGenerator():
    global zmp_refx, zmp_refy, com_refx, com_refy, time, capturePoint_refx, capturePoint_refy, com_refdx, com_refdy, com_refddx, com_refddy, total_tick, current_step_num, COM_x_1, COM_y_1, ux_1, uy_1
    capturePoint_ox = np.zeros(foot_step_number + 3)
    capturePoint_oy = np.zeros(foot_step_number + 3)
    zmp_dx = np.zeros(foot_step_number + 2)
    zmp_dy = np.zeros(foot_step_number + 2)
    b_offset = np.zeros(foot_step_number + 2)
    capturePoint_offsetx = np.zeros(foot_step_number + 3)
    capturePoint_offsety = np.zeros(foot_step_number + 3)

    total_tick = t_total * (foot_step_number + 1) + t_temp - 1

    COM_x_1 = np.asmatrix(np.zeros((3, int(total_tick + 1))))
    COM_y_1 = np.asmatrix(np.zeros((3, int(total_tick + 1))))
    COM_x_1[0,0] = PELV_tran_init[0]
    COM_y_1[0,0] = PELV_tran_init[1]
    ux_1 = np.asmatrix(np.zeros((int(total_tick + 1),1)))
    uy_1 = np.asmatrix(np.zeros((int(total_tick + 1),1)))

    time = np.zeros(int(total_tick))

    for i in range(0, foot_step_number + 2):
        b_offset[i] = math.exp(wn * t_total_t)

    for i in range(0, foot_step_number + 2):
        capturePoint_offsety[i] = 0.00
        capturePoint_offsetx[i] = 0.00
    
    capturePoint_ox[0] = PELV_tran_init[0]
    capturePoint_oy[0] = PELV_tran_init[1]
    capturePoint_ox[foot_step_number + 1] = (foot_step[foot_step_number - 1, 0] + foot_step[foot_step_number - 2, 0]) / 2 + capturePoint_offsetx[foot_step_number + 1]
    capturePoint_oy[foot_step_number + 1] = 0.0
    capturePoint_ox[foot_step_number + 2] = (foot_step[foot_step_number - 1, 0] + foot_step[foot_step_number - 2, 0]) / 2 + capturePoint_offsetx[foot_step_number + 2]
    capturePoint_oy[foot_step_number + 2] = 0.0

    for i in range(0,foot_step_number):
        if i == 0:
            capturePoint_ox[1] = PELV_tran_init[0] + capturePoint_offsetx[1]
            capturePoint_oy[1] = LF_tran[1] - capturePoint_offsety[1]
         
        else:
            if i % 2 == 0:
                capturePoint_ox[i + 1] = foot_step[i - 1, 0] + capturePoint_offsetx[i + 1]
                capturePoint_oy[i + 1] = foot_step[i - 1, 1] - capturePoint_offsety[i + 1]
            else:
                capturePoint_ox[i + 1] = foot_step[i - 1, 0] + capturePoint_offsetx[i + 1]
                capturePoint_oy[i + 1] = foot_step[i - 1, 1] + capturePoint_offsety[i + 1]   

    for i in range(0, foot_step_number + 2):
        zmp_dx[i] = capturePoint_ox[i + 1] / (1 - b_offset[i]) - (b_offset[i] * capturePoint_ox[i]) / (1 - b_offset[i])
        zmp_dy[i] = capturePoint_oy[i + 1] / (1 - b_offset[i]) - (b_offset[i] * capturePoint_oy[i]) / (1 - b_offset[i])

    capturePoint_refx = np.zeros(int(t_total * (foot_step_number + 1) + t_temp - 1))
    capturePoint_refy = np.zeros(int(t_total * (foot_step_number + 1) + t_temp - 1))
    com_refx = np.zeros(int(t_total * (foot_step_number + 1) + t_temp - 1))
    com_refy = np.zeros(int(t_total * (foot_step_number + 1) + t_temp - 1))
    com_refdx = np.zeros(int(t_total * (foot_step_number + 1) + t_temp - 1))
    com_refdy = np.zeros(int(t_total * (foot_step_number + 1) + t_temp - 1))
    com_refddx = np.zeros(int(t_total * (foot_step_number + 1) + t_temp - 1))
    com_refddy = np.zeros(int(t_total * (foot_step_number + 1) + t_temp - 1))
    zmp_refx = np.zeros(int(t_total * (foot_step_number + 1) + t_temp - 1)) 
    zmp_refy = np.zeros(int(t_total * (foot_step_number + 1) + t_temp - 1))
   
    for i in range(0,(int(t_total * (foot_step_number + 1) + t_temp - 1))):
        time[i] = i
        
        if i < t_temp - 1:
            current_step = int(i / (t_temp + t_total))
            if t_temp - t_total <= i:
                tick = (i - (t_temp - t_total - 1)) / float(hz)
            else:
                tick = i/float(hz)

            capturePointChange = int(i / (t_temp - 1))
        else:
            current_step = int((i - t_temp - t_total) / (t_total)) + 1
            capturePointChange = int((i - t_temp + 1) / (t_total)) + 1
            tick = i/float(hz) - float(t_total) * (int(capturePointChange) - 1) / float(hz) - (float(t_temp) - 1) / float(hz)

        if (int(capturePointChange) == foot_step_number + 1 and tick > (t_total) / hz) == False:
            if int(capturePointChange) == foot_step_number + 2:
                capturePoint_refx[i] = math.exp(wn * tick) * capturePoint_ox[int(capturePointChange) - 1] + (1 - math.exp(wn * tick)) * zmp_dx[int(capturePointChange) - 1]
                capturePoint_refy[i] = math.exp(wn * tick) * capturePoint_oy[int(capturePointChange) - 1] + (1 - math.exp(wn * tick)) * zmp_dy[int(capturePointChange) - 1]
            else:
                capturePoint_refx[i] = math.exp(wn * tick) * capturePoint_ox[int(capturePointChange)] + (1 - math.exp(wn * tick)) * zmp_dx[int(capturePointChange)]
                capturePoint_refy[i] = math.exp(wn * tick) * capturePoint_oy[int(capturePointChange)] + (1 - math.exp(wn * tick)) * zmp_dy[int(capturePointChange)]
        else:
            capturePoint_refx[i] = math.exp(wn * t_total / hz) * capturePoint_ox(int(capturePointChange)) + (1 - math.exp(wn * t_total / hz)) * zmp_dx(int(capturePointChange))
            capturePoint_refy[i] = math.exp(wn * t_total / hz) * capturePoint_oy(int(capturePointChange)) + (1 - math.exp(wn * t_total / hz)) * zmp_dy(int(capturePointChange))

        if int(capturePointChange) == 0 and i < t_temp - t_total:
            capturePoint_refx[i] = capturePoint_ox[0]
            capturePoint_refy[i] = capturePoint_oy[0]
        elif int(capturePointChange) == 0 and t_temp - t_total <= i:
            capturePoint_refx[i] = math.exp(wn * tick) * capturePoint_ox[int(capturePointChange)] + (1 - math.exp(wn * tick)) * zmp_dx[int(capturePointChange)]
            capturePoint_refy[i] = math.exp(wn * tick) * capturePoint_oy[int(capturePointChange)] + (1 - math.exp(wn * tick)) * zmp_dy[int(capturePointChange)]
        
        if i == 0:
            zmp_refx[0] = PELV_tran_init[0]
            zmp_refy[0] = PELV_tran_init[1]
        else:
            zmp_refx[i] = (capturePoint_refx[i - 1]) - (capturePoint_refx[i] - capturePoint_refx[i - 1]) * hz / (wn)
            zmp_refy[i] = (capturePoint_refy[i - 1]) - (capturePoint_refy[i] - capturePoint_refy[i - 1]) * hz / (wn)
    

        '''    
        if (i>t_temp):
            tick = (i-t_temp)/hz - t_total_t * current_step_num
            if(current_step_num == 0):
                A = foot_step[current_step_num + 1, 1] - 0.02
                Ky = (A * t_double * wn * np.tanh(wn*(t_total_t/2 - t_double)))/(1+t_double*wn*np.tanh(wn *(t_total_t/2 - t_double)))
                slopey = Ky/t_double
                B = (foot_step[current_step_num, 0] - RF_tran[0])/2
                Kx = (B * t_double * wn )/(t_double*wn+np.tanh(wn *(t_total_t/2 - t_double)))
                slopex = Kx/t_double
                if(tick <= t_double):
                    zmp_refy[i] = slopey * tick
                    zmp_refx[i] = slopex * tick + PELV_tran_init[0] 
                elif(t_double < tick and tick < t_total_t - t_double):
                    zmp_refy[i] = A
                    zmp_refx[i] = B + PELV_tran_init[0]
                else:
                    zmp_refy[i] = slopey * (t_total_t - tick)
                    zmp_refx[i] = (2*B - Kx) + slopex * (tick - (t_total_t - t_double)) + PELV_tran_init[0]
            elif(current_step_num != foot_step_number):
                if(foot_step[current_step_num - 1, 1] < 0):
                    A = foot_step[current_step_num - 1, 1] + 0.02
                else:
                    A = foot_step[current_step_num - 1, 1] - 0.02
                Ky = (A * t_double * wn * np.tanh(wn*(t_total_t/2 - t_double)))/(1+t_double*wn*np.tanh(wn *(t_total_t/2 - t_double)))
                slopey = Ky/t_double
                if(current_step_num == 1):
                    B = (foot_step[current_step_num, 0] - RF_tran[0])/2   
                else:
                    B = (foot_step[current_step_num, 0] - foot_step[current_step_num - 1, 0])/2
                Kx = (B * t_double * wn )/(t_double*wn+np.tanh(wn *(t_total_t/2 - t_double)))
                slopex = Kx/t_double
                if(tick <= t_double):
                    zmp_refy[i] = slopey * tick
                    zmp_refx[i] = slopex * tick + PELV_tran_init[0] + 2 * (current_step_num - 1) * B
                elif(t_double < tick and tick < t_total_t - t_double):
                    zmp_refy[i] = A
                    zmp_refx[i] = B + PELV_tran_init[0] + 2 * (current_step_num - 1) * B
                else:
                    zmp_refy[i] = slopey * (t_total_t - tick)
                    zmp_refx[i] = (2*B - Kx) + slopex * (tick - (t_total_t - t_double)) + PELV_tran_init[0] + 2 * (current_step_num - 1)* B
                if(i > t_total * (foot_step_number - 2) + t_temp - 1):
                   zmp_refx[i] = zmp_refx[int(t_total * (foot_step_number - 2)  + t_temp - 1)]
            else:
                zmp_refx[i] = zmp_refx[int(t_total * (foot_step_number - 2)  + t_temp - 1)]
            if(tick  > t_total_t - 1.5/hz) and i > t_temp:
                current_step_num = current_step_num + 1
        else:
            zmp_refx[i] = PELV_tran_init[0]
            zmp_refy[i] = PELV_tran_init[1]     
    current_step_num = 0
    '''
    
    
def calculatePreviewControlParams(A, B, C, Q, R, N):
    P = scipy.linalg.solve_discrete_are(A, B, C.T*Q*C, R)
    K = (R + B.T*P*B).I*(B.T*P*A)

    f = np.zeros((1, N))
    for i in range(N):
        f[0,i] = (R+B.T*P*B).I*B.T*(((A-B*K).T)**i)*C.T*Q

    return K, f

def calculatePreviewControlParams2(A, B, C, Q, R, N):
    C_dot_A = C*A
    C_dot_B = C*B

    A_tilde = np.matrix([[1, C_dot_A[0,0], C_dot_A[0,1], C_dot_A[0,2]],
                            [0, A[0,0], A[0,1], A[0,2]],
                            [0, A[1,0], A[1,1], A[1,2]],
                            [0, A[2,0], A[2,1], A[2,2]]])
    B_tilde = np.matrix([[C_dot_B[0,0]],
                            [B[0,0]],
                            [B[1,0]],
                            [B[2,0]]])
    C_tilde = np.matrix([[1, 0, 0, 0]])

    P_tilde = scipy.linalg.solve_discrete_are(A_tilde, B_tilde, C_tilde.T*Q*C_tilde, R)
    K_tilde = (R + B_tilde.T*P_tilde*B_tilde).I*(B_tilde.T*P_tilde*A_tilde)

    Ks = K_tilde[0, 0]
    Kx = K_tilde[0, 1:]

    Ac_tilde = A_tilde - B_tilde*K_tilde

    G = np.zeros((1, N))

    G[0] = -Ks
    I_tilde = np.matrix([[1],[0],[0],[0]])
    X_tilde = -Ac_tilde.T*P_tilde*I_tilde

    for i in range(N):
        G[0,i] = (R + B_tilde.T*P_tilde*B_tilde).I*(B_tilde.T)*X_tilde
        X_tilde = Ac_tilde.T*X_tilde

    return Ks, Kx, G

    
def comGenerator():
    global current_step_num
    for i in range(0, int(t_total * (foot_step_number + 1) + t_temp - 1)):
        if (i > t_temp):
            tick = (i-t_temp)/hz - t_total_t * current_step_num
            if(current_step_num == 0):
                A = foot_step[current_step_num + 1, 1] - 0.02
                Ky = (A * t_double * wn * np.tanh(wn*(t_total_t/2 - t_double)))/(1+t_double*wn*np.tanh(wn *(t_total_t/2 - t_double)))
                slopey = Ky/t_double
                B = (foot_step[current_step_num, 0] - RF_tran[0])/2
                Kx = (B * t_double * wn )/(t_double*wn+np.tanh(wn *(t_total_t/2 - t_double)))                
                slopex = Kx/t_double
                Cy1 = Ky - A
                Cy2 = Ky/(t_double * wn)
                if(tick <= t_double):
                    com_refy[i] = slopey * tick
                    com_refx[i] = slopex * tick + PELV_tran_init[0] 
                elif(t_double < tick and tick < t_total_t - t_double):
                    com_refy[i] = Cy1 * np.cosh(wn * (tick - t_double)) + Cy2 * np.sinh(wn * (tick - t_double)) + A
                    com_refx[i] = B + PELV_tran_init[0]
                else:
                    com_refy[i] = slopey * (t_total_t - tick)
                    com_refx[i] = (2*B - Kx) + slopex * (tick - (t_total_t - t_double)) + PELV_tran_init[0]
            elif(current_step_num != foot_step_number):
                if(foot_step[current_step_num - 1, 1] < 0):
                    A = foot_step[current_step_num - 1, 1] + 0.02
                else:
                    A = foot_step[current_step_num - 1, 1] - 0.02
                Ky = (A * t_double * wn * np.tanh(wn*(t_total_t/2 - t_double)))/(1+t_double*wn*np.tanh(wn *(t_total_t/2 - t_double)))
                slopey = Ky/t_double
                if(current_step_num == 1):
                    B = (foot_step[current_step_num, 0] - RF_tran[0])/2   
                else:
                    B = (foot_step[current_step_num, 0] - foot_step[current_step_num - 1, 0])/2
                Kx = (B * t_double * wn )/(t_double*wn+np.tanh(wn *(t_total_t/2 - t_double)))
                slopex = Kx/t_double
                Cy1 = Ky - A
                Cy2 = Ky/(t_double * wn)
                Cx1 = Kx - B
                Cx2 = Kx/(t_double * wn)
                if(tick <= t_double):
                    com_refy[i] = slopey * tick
                    com_refx[i] = slopex * tick + PELV_tran_init[0] + 2 * (current_step_num - 1)* B
                elif(t_double < tick and tick < t_total_t - t_double):
                    com_refy[i] = Cy1 * np.cosh(wn * (tick - t_double)) + Cy2 * np.sinh(wn * (tick - t_double)) + A
                    com_refx[i] = Cx1 * np.cosh(wn * (tick - t_double)) + Cx2 * np.sinh(wn * (tick - t_double)) + B + 2 * (current_step_num - 1)* B + PELV_tran_init[0]
                else:
                    com_refy[i] = slopey * (t_total_t - tick)
                    com_refx[i] = (2*B - Kx) + slopex * (tick - (t_total_t - t_double)) + PELV_tran_init[0] + 2 * (current_step_num - 1)* B
                if(i > t_total * (foot_step_number - 2) + t_temp - 1):
                    com_refx[i] = com_refx[int(t_total * (foot_step_number - 2)  + t_temp - 1)]
            else:
                com_refx[i] = com_refx[int(t_total * (foot_step_number - 2)  + t_temp - 1)]
            
            if(tick  > t_total_t - 1.5/hz ) and i > t_temp:
                current_step_num = current_step_num + 1
            
            '''
            com_refx[i] = wn / hz * capturePoint_refx[i] + (1 - wn / hz) * com_refx[i - 1]
            com_refy[i] = wn / hz * capturePoint_refy[i] + (1 - wn / hz) * com_refy[i - 1]

            com_refdx[i] = (com_refx[i] - com_refx[i - 1]) * hz
            com_refdy[i] = (com_refy[i] - com_refy[i - 1]) * hz
            com_refddx[i] = (com_refdx[i] - com_refdx[i - 1]) * hz
            com_refddy[i] = (com_refdy[i] - com_refdy[i - 1]) * hz
            '''
            
            if(i > int(total_tick) - N_preview):
                for j in range(0,N_preview):
                    if(i + j >= int(total_tick)):
                        ZMP_x_preview[j] = zmp_refx[int(total_tick) - 1]
                        ZMP_y_preview[j] = zmp_refy[int(total_tick) - 1]
                    else:
                        ZMP_x_preview[j] = zmp_refx[i + j]
                        ZMP_y_preview[j] = zmp_refy[i + j]
            else:
                ZMP_x_preview = np.asmatrix(zmp_refx[i:i+N_preview]).T
                ZMP_y_preview = np.asmatrix(zmp_refy[i:i+N_preview]).T

            ZMP_x = C_lipm*COM_x_1[:,i]
            ZMP_y = C_lipm*COM_y_1[:,i]

            ux_1[i] = -K_lipm*COM_x_1[:, i] + f_lipm*ZMP_x_preview
            uy_1[i] = -K_lipm*COM_y_1[:, i] + f_lipm*ZMP_y_preview

            COM_x_1[:,i+1] = A_lipm*COM_x_1[:, i] + B_lipm*ux_1[i]
            COM_y_1[:,i+1] = A_lipm*COM_y_1[:, i] + B_lipm*uy_1[i]
            com_refx[i] = COM_x_1[0,i]
            com_refy[i] = COM_y_1[0,i]
            com_refdx[i] = COM_x_1[1,i]
            com_refdy[i] = COM_y_1[1,i]
            com_refddx[i] = COM_x_1[2,i]
            com_refddy[i] = COM_y_1[2,i]

        else:
            '''
            com_refx[i] = PELV_tran_init[0]
            com_refy[i] = PELV_tran_init[1]
            COM_x_1[0,i + 1] = PELV_tran_init[0]
            COM_y_1[0,i + 1] = PELV_tran_init[1]
            '''
            com_refx[i] = COM_tran_init[0]
            com_refy[i] = COM_tran_init[1]
            COM_x_1[0,i + 1] = COM_tran_init[0]
            COM_y_1[0,i + 1] = COM_tran_init[1]
            
            com_refdx[i] = 0.0
            com_refdy[i] = 0.0
            com_refddx[i] = 0.0
            com_refddy[i] = 0.0
            #com_refy[i] = wn / hz * capturePoint_refy[i] + (1 - wn / hz) * com_refy[i - 1]
            #com_refdy[i] = (com_refy[i] - com_refy[i - 1]) * hz
            #com_refddy[i] = (com_refdy[i] - com_refdy[i - 1]) * hz
    
    current_step_num = 0                

def swingFootGenerator():
    global lfoot, rfoot, lfootd, rfootd, phase_variable
    phase_variable = np.zeros(int(total_tick))
    lfoot = np.zeros((int(total_tick),3))
    rfoot = np.zeros((int(total_tick),3))

    lfootd = np.zeros((int(total_tick),3))
    rfootd = np.zeros((int(total_tick),3))
    for i in range(0, int(t_total * (foot_step_number + 1) + t_temp - 1)):
        phase_variable[i] = 1
        if (i < t_start_real + t_double_1):
            lfoot[i,1] = LF_tran[1]
            rfoot[i,1] = RF_tran[1]
            lfoot[i,0] = LF_tran[0]
            rfoot[i,0] = RF_tran[0]
            lfoot[i,2] = LF_tran[2]
            rfoot[i,2] = RF_tran[2]
        elif (i < t_start_real + t_double_1 + t_total):
            if (foot_step[1, 6] == 1):
                lfoot[i,1] = foot_step[0, 1]
                rfoot[i,1] = RF_tran[1]
                lfoot[i,0] = foot_step[0, 0]
                rfoot[i,0] = RF_tran[0]
                lfoot[i,2] = LF_tran[2]
                rfoot[i,2] = RF_tran[2]
            else:
                lfoot[i,1] = LF_tran[1]
                rfoot[i,1] = foot_step[0, 1]
                lfoot[i,0] = LF_tran[0]
                rfoot[i,0] = foot_step[0, 0]
                lfoot[i,2] = LF_tran[2]
                rfoot[i,2] = RF_tran[2]
        else:
            j = int((i - t_temp) / t_total)
            if (j == 1):
                if (i <= t_start_real + t_double_1 + t_total * 2):
                    if (foot_step[j, 6] == 1):
                        lfoot[i,0] = foot_step[j - 1, 0]
                        lfoot[i,1] = foot_step[j - 1, 1]
                        lfoot[i,2] = LF_tran[2]
                        rfoot[i,1] = RF_tran[1]
                        rfoot[i,0] = quinticSpline(i, t_start_real + t_total + t_double_1, t_start + t_total * 2 - t_rest_2 - t_double_2, RF_tran[0], 0.0, 0.0, foot_step[j, 0], 0.0, 0.0)
                        rfootd[i,0] = quinticSplineDot(i, t_start_real + t_total + t_double_1, t_start + t_total * 2 - t_rest_2 - t_double_2, RF_tran[0], 0.0, 0.0, foot_step[j, 0], 0.0, 0.0) * hz
                        
                        if (i < t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0):
                            rfoot[i,2] = quinticSpline(i, t_start_real + t_total + t_double_1, t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2) / 2, RF_tran[2], 0.0, 0.0, RF_tran[2] + foot_height, 0.0, 0.0)               
                            rfootd[i,2] = quinticSplineDot(i, t_start_real + t_total + t_double_1, t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2) / 2, RF_tran[2], 0.0, 0.0, RF_tran[2] + foot_height, 0.0, 0.0) * hz            
                        else:
                            rfoot[i,2] = quinticSpline(i, t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0, t_start + t_total + t_total - t_rest_2 - t_double_2, RF_tran[2] + foot_height, 0.0, 0.0, RF_tran[2], 0.0, 0.0)
                            rfootd[i,2] = quinticSplineDot(i, t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0, t_start + t_total + t_total - t_rest_2 - t_double_2, RF_tran[2] + foot_height, 0.0, 0.0, RF_tran[2], 0.0, 0.0) * hz
                        if(i >= t_start_real + t_total + t_double_1) and ( i <= t_start + t_total + t_total - t_rest_2 - t_double_2):
                            phase_variable[i] = 3
                    else:
                        rfoot[i,0] = foot_step[j - 1, 0]
                        rfoot[i,1] = foot_step[j - 1, 1]
                        rfoot[i,2] = RF_tran[2]

                        lfoot[i,1] = LF_tran[1]
                        lfoot[i,0] = quinticSpline(i, t_start_real + t_total + t_double_1, t_start + t_total * 2 - t_rest_2 - t_double_2, LF_tran[0], 0.0, 0.0, foot_step[j, 0], 0.0, 0.0)
                        lfootd[i,0] = quinticSplineDot(i, t_start_real + t_total + t_double_1, t_start + t_total * 2 - t_rest_2 - t_double_2, LF_tran[0], 0.0, 0.0, foot_step[j, 0], 0.0, 0.0) * hz

                        if (i < t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0):
                            lfoot[i,2] = quinticSpline(i, t_start_real + t_total + t_double_1, t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2, LF_tran[2], 0.0, 0.0, LF_tran[2] + foot_height, 0.0, 0.0)
                            lfootd[i,2] = quinticSplineDot(i, t_start_real + t_total + t_double_1, t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2, LF_tran[2], 0.0, 0.0, LF_tran[2] + foot_height, 0.0, 0.0) * hz
                        else:
                            lfoot[i,2] = quinticSpline(i, t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0, t_start + t_total + t_total - t_rest_2 - t_double_2, LF_tran[2] + foot_height, 0.0, 0.0, LF_tran[2], 0.0, 0.0)
                            lfootd[i,2] = quinticSplineDot(i, t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0, t_start + t_total + t_total - t_rest_2 - t_double_2, LF_tran[2] + foot_height, 0.0, 0.0, LF_tran[2], 0.0, 0.0) * hz
            
                        if(i >= t_start_real + t_total + t_double_1) and ( i <= t_start + t_total + t_total - t_rest_2 - t_double_2):
                            phase_variable[i] = 2

            elif (j > 1 and j < foot_step_number):
                if (i <= t_start + t_double_1 + t_total * (j) and i >= t_start + t_total * (j) -1):    
                    if (foot_step[j, 6] == 1):
                        rfoot[i,0] = foot_step[j - 2, 0]
                        lfoot[i,0] = foot_step[j - 1, 0]
                        rfoot[i,1] = foot_step[j - 2, 1]
                        lfoot[i,1] = foot_step[j - 1, 1]
                        lfoot[i,2] = LF_tran[2]
                        rfoot[i,2] = RF_tran[2]
                    else:
                        lfoot[i,0] = foot_step[j - 2, 0]
                        rfoot[i,0] = foot_step[j - 1, 0]
                        lfoot[i,1] = foot_step[j - 2, 1]
                        rfoot[i,1] = foot_step[j - 1, 1]
                        rfoot[i,2] = RF_tran[2]
                        lfoot[i,2] = LF_tran[2]
                else:
                    if (foot_step[j, 6] == 1):
                        rfoot[i,1] = foot_step[j - 2, 1]
                        lfoot[i,1] = foot_step[j - 1, 1]
                        lfoot[i,0] = foot_step[j - 1, 0]
                        lfoot[i,2] = LF_tran[2]
                        rfoot[i,0] = quinticSpline(i, t_start_real + t_total * j + t_double_1, t_start + t_total * (j + 1) - t_rest_2 - t_double_2, foot_step[j - 2, 0], 0.0, 0.0, foot_step[j, 0], 0.0, 0.0)
                        rfootd[i,0] = quinticSplineDot(i, t_start_real + t_total * j + t_double_1, t_start + t_total * (j + 1) - t_rest_2 - t_double_2, foot_step[j - 2, 0], 0.0, 0.0, foot_step[j, 0], 0.0, 0.0) * hz
                        
                        if (i < t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0):
                            rfoot[i,2] = quinticSpline(i, t_start_real + t_total * j + t_double_1, t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2, LF_tran[2], 0.0, 0.0, RF_tran[2] + foot_height, 0.0, 0.0)
                            rfootd[i,2] = quinticSplineDot(i, t_start_real + t_total * j + t_double_1, t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2, LF_tran[2], 0.0, 0.0, RF_tran[2] + foot_height, 0.0, 0.0) * hz
                        else:
                            rfoot[i,2] = quinticSpline(i, t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0, t_start + t_total * j + t_total - t_rest_2 - t_double_2, RF_tran[2] + foot_height, 0.0, 0.0, RF_tran[2], 0.0, 0.0)
                            rfootd[i,2] = quinticSplineDot(i, t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0, t_start + t_total * j + t_total - t_rest_2 - t_double_2, RF_tran[2] + foot_height, 0.0, 0.0, RF_tran[2], 0.0, 0.0) * hz
                        if(i >= t_start_real + t_total * j + t_double_1) and ( i <=  t_start + t_total * j + t_total - t_rest_2 - t_double_2):
                            phase_variable[i] = 3
                    else:
                        lfoot[i,1] = foot_step[j - 2, 1]
                        rfoot[i,1] = foot_step[j - 1, 1]
                        rfoot[i,0] = foot_step[j - 1, 0]
                        rfoot[i,2] = RF_tran[2]
                        lfoot[i,0] = quinticSpline(i, t_start_real + t_total * j + t_double_1, t_start + t_total * (j + 1) - t_rest_2 - t_double_2, foot_step[j - 2, 0], 0.0, 0.0, foot_step[j, 0], 0.0, 0.0)
                        lfootd[i,0] = quinticSplineDot(i, t_start_real + t_total * j + t_double_1, t_start + t_total * (j + 1) - t_rest_2 - t_double_2, foot_step[j - 2, 0], 0.0, 0.0, foot_step[j, 0], 0.0, 0.0) * hz
            
                        if (i < t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0):
                            lfoot[i,2] = quinticSpline(i, t_start_real + t_total * j + t_double_1, t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2, LF_tran[2], 0.0, 0.0, LF_tran[2] + foot_height, 0.0, 0.0)
                            lfootd[i,2] = quinticSplineDot(i, t_start_real + t_total * j + t_double_1, t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2, LF_tran[2], 0.0, 0.0, LF_tran[2] + foot_height, 0.0, 0.0) * hz
                        else:
                            lfoot[i,2] = quinticSpline(i, t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0, t_start + t_total * j + t_total - t_rest_2 - t_double_2, LF_tran[2] + foot_height, 0.0, 0.0, LF_tran[2], 0.0, 0.0)
                            lfootd[i,2] = quinticSplineDot(i, t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0, t_start + t_total * j + t_total - t_rest_2 - t_double_2, LF_tran[2] + foot_height, 0.0, 0.0, LF_tran[2], 0.0, 0.0) * hz
                        if(i >= t_start_real + t_total * j + t_double_1) and ( i <=  t_start + t_total * j + t_total - t_rest_2 - t_double_2):
                            phase_variable[i] = 2
            elif (j == foot_step_number):
                if (i >= t_start + t_total * (j)-1):
                    if (foot_step[foot_step_number - 1, 6] == 1):
                        rfoot[i,0] = foot_step[foot_step_number - 1, 0]
                        lfoot[i,0] = foot_step[foot_step_number - 2, 0]
                        rfoot[i,1] = foot_step[foot_step_number - 1, 1]
                        lfoot[i,1] = foot_step[foot_step_number - 2, 1]
                    else:
                        lfoot[i,0] = foot_step[foot_step_number - 1, 0]
                        rfoot[i,0] = foot_step[foot_step_number - 2, 0]
                        lfoot[i,1] = foot_step[foot_step_number - 1, 1]
                        rfoot[i,1] = foot_step[foot_step_number - 2, 1]
                    lfoot[i,2] = LF_tran[2]
                    rfoot[i,2] = RF_tran[2]

def contactRedistribution(eta_cust, footwidth, footlength, staticFrictionCoeff, ratio_x, ratio_y, P1, P2, F12):
    global ResultantForce, ForceRedistribution, eta
    W1 = np.zeros((6,12))

    W1[0:6, 0:6] = np.identity(6)
    W1[0:6, 6:12] = np.identity(6)

    W1[3:6, 0:3] = skew(P1)
    W1[3:6, 6:9] = skew(P2)

    ResultantForce = np.matmul(W1, F12)
    
    eta_lb = 1.0 - eta_cust
    eta_ub = eta_cust

    A = (P1[2] - P2[2]) * ResultantForce[1] - (P1[1] - P2[1]) * ResultantForce[2]    
    B = ResultantForce[3] + P2[2] * ResultantForce[1] - P2[1] * ResultantForce[2]
    C = ratio_y * footwidth / 2.0 * abs(ResultantForce[2])
    a = A * A
    b = 2.0 * A * B
    c = B * B - C * C
    sol_eta1 = (-b + np.sqrt(b * b - 4.0 * a * c)) / 2.0 / a
    sol_eta2 = (-b - np.sqrt(b * b - 4.0 * a * c)) / 2.0 / a

    if (sol_eta1 > sol_eta2):
        if (sol_eta1 < eta_ub):
            eta_ub = sol_eta1
        
        if (sol_eta2 > eta_lb):
            eta_lb = sol_eta2
    else: 
        if (sol_eta2 < eta_ub):
            eta_ub = sol_eta2
        if (sol_eta1 > eta_lb):
            eta_lb = sol_eta1

    A = -(P1[2] - P2[2]) * ResultantForce[0] + (P1[0] - P2[0]) * ResultantForce[2]
    B = ResultantForce[4] - P2[2] * ResultantForce[0] + P2[0] * ResultantForce[2]
    C = ratio_x * footlength / 2.0 * abs(ResultantForce[2])
    a = A * A
    b = 2.0 * A * B
    c = B * B - C * C
    sol_eta1 = (-b + np.sqrt(b * b - 4.0 * a * c)) / 2.0 / a
    sol_eta2 = (-b - np.sqrt(b * b - 4.0 * a * c)) / 2.0 / a

    if (sol_eta1 > sol_eta2):
        if (sol_eta1 < eta_ub):
            eta_ub = sol_eta1

        if (sol_eta2 > eta_lb):
            eta_lb = sol_eta2
        
    else: 
        if (sol_eta2 < eta_ub):
            eta_ub = sol_eta2

        if (sol_eta1 > eta_lb):
            eta_lb = sol_eta1

    A = -(P1[0] - P2[0]) * ResultantForce[1] + (P1[1] - P2[1]) * ResultantForce[0]
    B = ResultantForce[5] + P2[1] * ResultantForce[0] - P2[0] * ResultantForce[1]
    C = staticFrictionCoeff * np.abs(ResultantForce[2])
    a = A * A
    b = 2.0 * A * B
    c = B * B - C * C
    sol_eta1 = (-b + np.sqrt(b * b - 4.0 * a * c)) / 2.0 / a
    sol_eta2 = (-b - np.sqrt(b * b - 4.0 * a * c)) / 2.0 / a
    if (sol_eta1 > sol_eta2): 
        if (sol_eta1 < eta_ub):
            eta_ub = sol_eta1
        if (sol_eta2 > eta_lb):
            eta_lb = sol_eta2
    else:
        if (sol_eta2 < eta_ub):
            eta_ub = sol_eta2
        if (sol_eta1 > eta_lb):
            eta_lb = sol_eta1

    eta_s = (-ResultantForce[3] - P2[2] * ResultantForce[1] + P2[1] * ResultantForce[2]) / ((P1[2] - P2[2]) * ResultantForce[1] - (P1[1] - P2[1]) * ResultantForce[2])

    eta = eta_s 
    if (eta_s > eta_ub):
        eta = eta_ub
    elif (eta_s < eta_lb):
        eta = eta_lb

    if ((eta > eta_cust) or (eta < 1.0 - eta_cust)):
        eta = 0.5

    ForceRedistribution = np.zeros(12)

    ForceRedistribution[0] = eta * ResultantForce[0]
    ForceRedistribution[1] = eta * ResultantForce[1]
    ForceRedistribution[2] = eta * ResultantForce[2]
    ForceRedistribution[3] = ((P1[2] - P2[2]) * ResultantForce[1] - (P1[1] - P2[1]) * ResultantForce[2]) * eta * eta + (ResultantForce[3] + P2[2] * ResultantForce[1] - P2[1] * ResultantForce[2]) * eta
    ForceRedistribution[4] = (-(P1[2] - P2[2]) * ResultantForce[0] + (P1[0] - P2[0]) * ResultantForce[2]) * eta * eta + (ResultantForce[4] - P2[2] * ResultantForce[0] + P2[0] * ResultantForce[2]) * eta
    ForceRedistribution[5] = (-(P1[0] - P2[0]) * ResultantForce[1] + (P1[1] - P2[1]) * ResultantForce[0]) * eta * eta + (ResultantForce[5] + P2[1] * ResultantForce[0] - P2[0] * ResultantForce[1]) * eta
    ForceRedistribution[6] = (1.0 - eta) * ResultantForce[0]
    ForceRedistribution[7] = (1.0 - eta) * ResultantForce[1]
    ForceRedistribution[8] = (1.0 - eta) * ResultantForce[2]
    ForceRedistribution[9] = (1.0 - eta) * (((P1[2] - P2[2]) * ResultantForce[1] - (P1[1] - P2[1]) * ResultantForce[2]) * eta + (ResultantForce[3] + P2[2] * ResultantForce[1] - P2[1] * ResultantForce[2]))
    ForceRedistribution[10] = (1.0 - eta) * ((-(P1[2] - P2[2]) * ResultantForce[0] + (P1[0] - P2[0]) * ResultantForce[2]) * eta + (ResultantForce[4] - P2[2] * ResultantForce[0] + P2[0] * ResultantForce[2]))
    ForceRedistribution[11] = (1.0 - eta) * ((-(P1[0] - P2[0]) * ResultantForce[1] + (P1[1] - P2[1]) * ResultantForce[0]) * eta + (ResultantForce[5] + P2[1] * ResultantForce[0] - P2[0] * ResultantForce[1]))    

def contactRedistributionWalking(command_torque, eta, ratio, supportFoot):
    '''
    global torque_contact, V2
    contact_dof_ = int(np.shape(robotJac)[0])
    
    if (contact_dof_ == 12):
        ContactForce_ = np.matmul(robotJcinvT[0:12,6:6+model.nq], command_torque) - robotPc

        P1_ = LFc_tran_cur - COM_tran_cur
        P2_ = RFc_tran_cur - COM_tran_cur

        Rotyaw = rotateWithZ(0)

        force_rot_yaw = np.zeros((12,12))

        for i in range(0,4):
            force_rot_yaw[i * 3:i * 3 +3 , i * 3:i * 3 +3] = Rotyaw
            
        ResultantForce_ = np.zeros(6)
        ResultRedistribution_ = np.zeros(12)
        F12 = np.matmul(force_rot_yaw, ContactForce_)
        eta_cust = 0.99
        foot_length = 0.26
        foot_width = 0.1         
        contactRedistribution(eta_cust, foot_width, foot_length, 1.0, 0.9, 0.9, np.matmul(Rotyaw, P1_), np.matmul(Rotyaw, P2_), F12)
        ResultantForce_ = ResultantForce
        ResultRedistribution_ = ForceRedistribution
        
        fc_redist_ = np.matmul(np.transpose(force_rot_yaw), ResultRedistribution_)

        desired_force = np.zeros(12)
        if (supportFoot == 0):
            right_master = 1.0
        else:
            right_master = 0.0

        if (right_master):
            desired_force[0:6] = -ContactForce_[0:6] + ratio * fc_redist_[0:6]  
            torque_contact = np.matmul(np.matmul(np.transpose(V2), np.linalg.inv(np.matmul(robotJcinvT[0:6, 6:model.nq], np.transpose(V2)))), desired_force[0:6])  
        else:
            desired_force[6:12] = -ContactForce_[6:12] + ratio * fc_redist_[6:12]
            torque_contact = np.matmul(np.matmul(np.transpose(V2), np.linalg.inv(np.matmul(robotJcinvT[6:12, 6:model.nq], np.transpose(V2)))), desired_force[6:12])

    else:
        torque_contact = np.zeros(33)

    return torque_contact    
    '''    
    global torque_contact, V2
    contact_dof_ = int(np.shape(robotJac)[0])

    if (contact_dof_ == 12):
        ContactForce_ = np.matmul(robotJcinvT, command_torque) - robotPc

        P1_ = LFc_tran_cur - COM_tran_cur
        P2_ = RFc_tran_cur - COM_tran_cur

        Rotyaw = rotateWithZ(0)

        force_rot_yaw = np.zeros((12,12))

        for i in range(0,4):
            force_rot_yaw[i * 3:i * 3 +3 , i * 3:i * 3 +3] = Rotyaw
            
        ResultantForce_ = np.zeros(6)
        ResultRedistribution_ = np.zeros(12)
        F12 = np.matmul(force_rot_yaw, ContactForce_)
        eta_cust = 0.99
        foot_length = 0.26
        foot_width = 0.1         
        contactRedistribution(eta_cust, foot_width, foot_length, 1.0, 0.9, 0.9, np.matmul(Rotyaw, P1_), np.matmul(Rotyaw, P2_), F12)
        ResultantForce_ = ResultantForce
        ResultRedistribution_ = ForceRedistribution
        
        fc_redist_ = np.matmul(np.transpose(force_rot_yaw), ResultRedistribution_)

        desired_force = np.zeros(12)
        if (supportFoot == 0):
            right_master = 1.0
        else:
            right_master = 0.0

        if (right_master):
            desired_force[0:6] = -ContactForce_[0:6] + ratio * fc_redist_[0:6]  
            torque_contact = np.matmul(np.matmul(np.transpose(V2), np.linalg.inv(np.matmul(robotJcinvT[0:6, 6:model.nq], np.transpose(V2)))), desired_force[0:6])  
        else:
            desired_force[6:12] = -ContactForce_[6:12] + ratio * fc_redist_[6:12]
            torque_contact = np.matmul(np.matmul(np.transpose(V2), np.linalg.inv(np.matmul(robotJcinvT[6:12, 6:model.nq], np.transpose(V2)))), desired_force[6:12])

    else:
        torque_contact = np.zeros(33)

    return torque_contact
                           
def inverseKinematics(time, LF_rot_c, RF_rot_c, PELV_rot_c, LF_tran_c, RF_tran_c, PELV_tran_c, HRR_tran_init_c, HLR_tran_init_c, HRR_rot_init_c, HLR_rot_init_c, PELV_tran_init_c, PELV_rot_init_c, CPELV_tran_init_c):
    global leg_q, leg_qdot, leg_qddot, leg_qs, leg_qdots, leg_qddots
    M_PI = 3.14159265358979323846
    if time == 0:
        leg_q = np.zeros(12)
        leg_qdot = np.zeros(12)
        leg_qddot = np.zeros(12)
        leg_qs = np.zeros((int(total_tick), 12))
        leg_qdots = np.zeros((int(total_tick), 12))
        leg_qddots = np.zeros((int(total_tick), 12))

    l_upper = 0.35
    l_lower = 0.35

    offset_hip_pitch = 0.0
    offset_knee_pitch = 0.0
    offset_ankle_pitch = 0.0

    lpt = np.subtract(PELV_tran_c, LF_tran_c)
    rpt = np.subtract(PELV_tran_c, RF_tran_c)
    lp = np.matmul(np.transpose(LF_rot_c), lpt)
    rp = np.matmul(np.transpose(RF_rot_c), rpt)
    
    PELF_rot = np.matmul(np.transpose(PELV_rot_c), LF_rot_c)
    PERF_rot = np.matmul(np.transpose(PELV_rot_c), RF_rot_c)

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

    lc = np.linalg.norm(lr)

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

    r_tl2 = np.matmul(np.matmul(np.matmul(PELF_rot, np.transpose(r_l4l5)),np.transpose(r_l3l4)),np.transpose(r_l2l3))
    leg_q[1] = np.arcsin(r_tl2[2, 1])

    c_lq5 = np.divide(-r_tl2[0, 1], np.cos(leg_q[1]))

    if c_lq5 > 1.0:
        c_lq5 = 1.0
    elif c_lq5 < -1.0:
        c_lq5 = -1.0
    
    leg_q[0] = -1 * np.arcsin(c_lq5)
    leg_q[2] = -1 * np.arcsin(r_tl2[2, 0] / np.cos(leg_q[1])) + offset_hip_pitch
    leg_q[3] = leg_q[3] - offset_knee_pitch
    leg_q[4] = leg_q[4] - offset_ankle_pitch

    rc = np.linalg.norm(rr)
    leg_q[9] = -1 * np.arccos((l_upper * l_upper + l_lower * l_lower - rc * rc) / (2 * l_upper * l_lower)) + M_PI

    r_ankle_pitch = np.arcsin((l_upper * np.sin(M_PI - leg_q[9])) / rc)
    leg_q[10] = -1 * np.arctan2(rr[0], np.sqrt(rr[1] * rr[1] + rr[2] * rr[2])) - r_ankle_pitch
    leg_q[11] = np.arctan2(rr[1], rr[2])
    r_tr2 = np.zeros((3,3))
    r_r2r3 = np.zeros((3,3))
    r_r3r4 = np.zeros((3,3))
    r_r4r5 = np.zeros((3,3))

    r_r2r3 = rotateWithY(leg_q[9])
    r_r3r4 = rotateWithY(leg_q[10])
    r_r4r5 = rotateWithX(leg_q[11])

    r_tr2 = np.matmul(np.matmul(np.matmul(PERF_rot, np.transpose(r_r4r5)),np.transpose(r_r3r4)),np.transpose(r_r2r3))
    leg_q[7] = np.arcsin(r_tr2[2,1])
    c_rq5 = -r_tr2[0, 1] / np.cos(leg_q[7])

    if c_rq5 > 1.0:
        c_rq5 = 1.0
    elif c_rq5 < -1.0:
        c_rq5 = -1.0 
    
    leg_q[6] = -1* np.arcsin(c_rq5)
    leg_q[8] = np.arcsin(r_tr2[2, 0] / np.cos(leg_q[7])) - offset_hip_pitch
    leg_q[9] = -1 * leg_q[9] + offset_knee_pitch
    leg_q[10] = -1 * leg_q[10] + offset_ankle_pitch

    leg_q[0] = leg_q[0] * (-1)
    leg_q[6] = leg_q[6] * (-1)
    leg_q[8] = leg_q[8] * (-1)
    leg_q[9] = leg_q[9] * (-1)
    leg_q[10] = leg_q[10] * (-1)

    leg_qs[time,:] = leg_q
    
    if(time == 0):
        leg_qdots[time,:] = np.zeros(12)
        leg_qddots[time,:] = np.zeros(12)
    else:
        leg_qdots[time,:] = np.subtract(leg_qs[time,:], leg_qs[time-1,:]) * hz
        leg_qddots[time,:] = np.subtract(leg_qdots[time,:], leg_qdots[time-1,:]) * hz
        
def comJacobianinverseKinematics(time, LF_rot_c, RF_rot_c, PELV_rot_c, LF_tran_c, RF_tran_c, PELV_tran_c, HRR_tran_init_c, HLR_tran_init_c, HRR_rot_init_c, HLR_rot_init_c, PELV_tran_init_c, PELV_rot_init_c, CPELV_tran_init_c):
    if(foot_step[current_step_num,6] == 1):
        r_c1 = np.subtract(COM_tran_cur, LF_tran_cur)
        skew_r_c1 = skew(r_c1)
        adjoint_ = np.zeros((6,6))
        adjoint_[0:3, 0:3] = PELV_rot
        adjoint_[3:6, 3:6] = PELV_rot
        j1 = np.matmul(adjoint_,LF_j[:,6:12])
        j2 = np.matmul(adjoint_,RF_j[:,12:18])

        j1_v = j1[0:3,0:6]
        j1_w = j1[3:6,0:6]

        skew_r2_r1 = skew(np.matmul(PELV_rot, np.subtract(LF_tran_cur,RF_tran_cur)))
        adjoint_21 = np.identity(6)
        adjoint_21[0:3, 3:6] = skew_r2_r1

        err_foot = np.zeros(3)
      #  err_foot = np.subtract(rfoot[time,:], RF_tran_cur)
        err_com = np.zeros(3)
      #  err_com[0] = com_refx[time] - COM_tran_cur[0]
      #  err_com[1] = com_refy[time] - COM_tran_cur[1]
      #  err_com[2] = - COM_tran_cur[2]
        
        I =  np.identity(3)
      #  err_com_w = matrixToRpy(np.matmul(I, np.transpose(PELV_rot_cur)))
      #  err_foot_w = matrixToRpy(np.matmul(I, np.transpose(RF_rot_cur)))

        err_foot_w = np.zeros(3)
        err_com_w = np.zeros(3)

        x2_d_dot_ = np.zeros(6)
        x2_d_dot_[0:3] = np.add(rfootd[time,0:3], 0.8 * err_foot)
        x2_d_dot_[3:6] = 0.8 * err_foot_w

        jlleg_com = Jcom[0:3, 6:12]
        jrleg_com = Jcom[0:3, 12:18]
        j_com_psem_ = np.add(np.add(np.add(-j1_v, np.matmul(skew_r_c1,j1_w)), jlleg_com), np.matmul(np.matmul(np.matmul(jrleg_com, np.linalg.inv(j2)),adjoint_21),j1))
        
        desired_u_dot = np.zeros(3)
        desired_u_dot[0] = com_refdx[time] + 0.8 * err_com[0]
        desired_u_dot[1] = com_refdy[time] + 0.8 * err_com[1]
        desired_u_dot[2] = 0.0 + 0.8 * err_com[2]
        desired_w_dot = 0.8 * err_com_w

        c_dot_psem_ = np.subtract(desired_u_dot,np.matmul(np.matmul(jrleg_com,np.linalg.inv(j2)),x2_d_dot_))

        j_total = np.zeros((6,6))
        j_total[0:3,0:6] = j_com_psem_
        j_total[3:6,0:6] = -j1_w
        
        c_total = np.zeros(6)
        c_total[0:3] = c_dot_psem_
        c_total[3:6] = desired_w_dot
        desired_leg_q_dot = np.zeros(12)
        desired_leg_q_dot[0:6] = np.matmul(scipy.linalg.pinv2(j_total),c_total)
        desired_leg_q_dot[6:12] = np.matmul(np.linalg.inv(j2),np.add(x2_d_dot_, np.matmul(np.matmul(adjoint_21,j1),desired_leg_q_dot[0:6])))

        if(time == 0):
            leg_qs[time, :] = q_init[7:19]
            leg_qs[time + 1, :] = q_init[7:19]
        else:
            leg_qs[time + 1, :] = leg_qs[time, : ] + desired_leg_q_dot/float(hz)
        
    else:
        r_c1 = np.subtract(COM_tran_cur, RF_tran_cur)
        skew_r_c1 = skew(r_c1)
        adjoint_ = np.zeros((6,6))
        adjoint_[0:3, 0:3] = PELV_rot
        adjoint_[3:6, 3:6] = PELV_rot
        j1 = np.matmul(adjoint_,RF_j[:,12:18])
        j2 = np.matmul(adjoint_,LF_j[:,6:12])

        j1_v = j1[0:3,0:6]
        j1_w = j1[3:6,0:6]

        skew_r2_r1 = skew(np.matmul(PELV_rot, np.subtract(RF_tran_cur,LF_tran_cur)))
        adjoint_21 = np.identity(6)
        adjoint_21[0:3, 3:6] = skew_r2_r1

        #err_foot = np.subtract(lfoot[time,:], LF_tran_cur)
        err_foot = np.zeros(3)
        err_com = np.zeros(3)
        
        err_foot_w = np.zeros(3)
        err_com_w = np.zeros(3)
        #err_com[0] = com_refx[time] - COM_tran_cur[0] 
        #err_com[1] = com_refy[time] - COM_tran_cur[1]
        #err_com[2] = - COM_tran_cur[2]

        I =  np.identity(3)
        #err_com_w = matrixToRpy(np.matmul(I, np.transpose(PELV_rot_cur)))
        #err_foot_w = matrixToRpy(np.matmul(I, np.transpose(LF_rot_cur)))

        x2_d_dot_ = np.zeros(6)
        x2_d_dot_[0:3] = np.add(lfootd[time,0:3], 0.8 * err_foot) 
        x2_d_dot_[3:6] = 0.8 * err_foot_w

        jlleg_com = Jcom[0:3, 6:12]
        jrleg_com = Jcom[0:3, 12:18]
        j_com_psem_ = np.add(np.add(np.add(-j1_v, np.matmul(skew_r_c1,j1_w)), jrleg_com), np.matmul(np.matmul(np.matmul(jlleg_com, np.linalg.inv(j2)),adjoint_21),j1))
        
        desired_u_dot = np.zeros(3)
        desired_u_dot[0] = com_refdx[time] + 0.8 * err_com[0]
        desired_u_dot[1] = com_refdy[time] + 0.8 * err_com[1]
        desired_u_dot[2] = 0.0 + 0.8 * err_com[2]
        
        desired_w_dot = 0.8 * err_com_w
        
        c_dot_psem_ = np.subtract(desired_u_dot,np.matmul(np.matmul(jlleg_com,np.linalg.inv(j2)),x2_d_dot_))
       
        j_total = np.zeros((6,6))
        j_total[0:3,0:6] = j_com_psem_
        j_total[3:6,0:6] = -j1_w
        

        c_total = np.zeros(6)
        c_total[0:3] = c_dot_psem_
        c_total[3:6] = desired_w_dot
        desired_leg_q_dot = np.zeros(12)
        desired_leg_q_dot[6:12] = np.matmul(scipy.linalg.pinv2(j_total),c_total)
        desired_leg_q_dot[0:6] = np.matmul(np.linalg.inv(j2),np.add(x2_d_dot_, np.matmul(np.matmul(adjoint_21,j1),desired_leg_q_dot[6:12])))

        if(time == 0):
            leg_qs[time, :] = q_init[7:19]
            leg_qs[time + 1, :] = q_init[7:19]
        else:
            leg_qs[time + 1, :] = leg_qs[time, : ] + desired_leg_q_dot/float(hz)
       
def modelInitialize():
    global model, foot_distance, data, LFframe_id, RFframe_id, PELVjoint_id, LHjoint_id, RHjoint_id, LFjoint_id, q_init, RFjoint_id, LFcframe_id, RFcframe_id, q, qdot, qddot, LF_tran, RF_tran, PELV_tran, LF_rot, RF_rot, PELV_rot, qdot_z, qddot_z, HRR_rot_init, HLR_rot_init, HRR_tran_init, HLR_tran_init, LF_rot_init, RF_rot_init, LF_tran_init, RF_tran_init, PELV_tran_init, PELV_rot_init, CPELV_tran_init, q_command, qdot_command, qddot_command, robotAginit, COM_tran_init, virtual_init, TranFVi, TranFRi, TranFLi, TranVRi, TranVLi
    #model = pinocchio.buildModelFromUrdf("/home/jhk/legAnalysis/dyros_tocabi_with_redhands_hipheavy.urdf",pinocchio.JointModelFreeFlyer())     
    model = pinocchio.buildModelFromUrdf("/home/jhk/legAnalysis/dyros_tocabi_with_redhands.urdf",pinocchio.JointModelFreeFlyer())     
    #model = pinocchio.buildModelFromUrdf("/home/jhk/legAnalysis/dyros_tocabi_with_redhands_footheavy.urdf",pinocchio.JointModelFreeFlyer())      
   
    LFframe_id = model.getFrameId("L_Foot_Link")
    RFframe_id = model.getFrameId("R_Foot_Link")
    PELVframe_id = model.getFrameId("Pelvis_Link")
    
    PELVjoint_id = model.getJointId("root_joint")
    LHjoint_id = model.getJointId("L_HipYaw_Joint")
    RHjoint_id = model.getJointId("R_HipYaw_Joint")
    RFjoint_id = model.getJointId("R_AnkleRoll_Joint")
    LFjoint_id = model.getJointId("L_AnkleRoll_Joint")

    contactPointLF = pinocchio.SE3.Identity()
    contactPointRF = pinocchio.SE3.Identity()
    
    contactPointLF.translation.T.flat = [0.03, 0, -0.1585]
    contactPointRF.translation.T.flat = [0.03, 0, -0.1585]

    model.addBodyFrame("LF_contact", LFjoint_id, contactPointLF, LFframe_id)
    model.addBodyFrame("RF_contact", RFjoint_id, contactPointRF, RFframe_id)

    LFcframe_id = model.getFrameId("LF_contact")
    RFcframe_id = model.getFrameId("RF_contact")

    data = model.createData()
    q = pinocchio.randomConfiguration(model)
    q_command = pinocchio.randomConfiguration(model)
    q_init = [0.0, 0.0, 0.7457, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.606, 1.572, -0.92, 0.0, 0.0, 0.0, -0.606, 1.572, -0.92, 0.0, 0.0, 0.0, 0.0, 0.2, 0.6, 1.5, -1.47, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -0.2, -0.6, -1.5, 1.47, 1.0, 0.0, 1.0, 0.0]
    #q_init = [0.0, 0.0, 0.814175, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.55, 1.26, -0.71, 0.0, 0.0, 0.0, -0.55, 1.26, -0.71, 0.0, 0.0, 0.0, 0.0, 0.2, 0.6, 1.5, -1.47, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -0.2, -0.6, -1.5, 1.47, 1.0, 0.0, 1.0, 0.0]
    for i in range(0, len(q)):
        q[i] = q_init[i]
        q_command[i] = q_init[i]

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
    pinocchio.centerOfMass(model,data,False)

    LF_tran = data.oMi[LFjoint_id].translation
    RF_tran = data.oMi[RFjoint_id].translation
    LF_rot = data.oMi[LFjoint_id].rotation
    RF_rot = data.oMi[RFjoint_id].rotation

    PELV_tran = np.add(data.oMi[PELVjoint_id].translation, model.inertias[PELVjoint_id].lever)
    PELV_rot = data.oMi[PELVjoint_id].rotation

    virtual_init = data.oMi[PELVjoint_id].translation

    LF_tran_init = data.oMi[LFjoint_id].translation
    RF_tran_init = data.oMi[RFjoint_id].translation
    HLR_tran_init = data.oMi[LHjoint_id].translation
    HRR_tran_init = data.oMi[RHjoint_id].translation
    LF_rot_init = data.oMi[LFjoint_id].rotation
    RF_rot_init = data.oMi[RFjoint_id].rotation
    HLR_rot_init = data.oMi[LHjoint_id].rotation
    HRR_rot_init = data.oMi[RHjoint_id].rotation

    COM_tran_init = data.com[0]
    PELV_tran_init = np.add(data.oMi[PELVjoint_id].translation, model.inertias[PELVjoint_id].lever)
    CPELV_tran_init = data.oMi[PELVjoint_id].translation 
    PELV_rot_init = data.oMi[PELVjoint_id].rotation

    virtual_init = data.oMf[PELVjoint_id].translation

    foot_distance = LF_tran_init - RF_tran_init

    TranFVi = np.zeros((4,4))
    TranFVi[0:3,0:3] = np.identity(3)
    TranFVi[0:3,3] = PELV_tran
    TranFVi[3,3] = 1.0

    TranFRi = np.zeros((4,4))
    TranFRi[0:3,0:3] = RF_rot_init
    TranFRi[0:3,3] = RF_tran_init
    TranFRi[3,3] = 1.0

    TranFLi = np.zeros((4,4))
    TranFLi[0:3,0:3] = LF_rot_init
    TranFLi[0:3,3] = LF_tran_init
    TranFLi[3,3] = 1.0

    TranVRi = np.matmul(np.linalg.inv(TranFVi),TranFRi)
    TranVLi = np.matmul(np.linalg.inv(TranFVi),TranFLi)

def modelUpdate(q_desired, qdot_desired, qddot_desired):
    global contactnum, M, G, COR, Minv, b, robotJac, robotdJac, LF_j, RF_j, LF_cj, RF_cj, LF_cdj, RF_cdj, robotLambdac, robotJcinvT, robotNc, robotPc, robotmuc, robotW, robothc, LF_tran_cur, RF_tran_cur, PELV_tran_cur, COM_tran_cur, RFc_tran_cur, LFc_tran_cur, robotWinv, robotCAM, Jcom, LF_rot_cur, RF_rot_cur, PELV_rot_cur, virtual_cur, virtual_dcur, virtual_ddcur, TranFV, TranFR, TranFL, TranVR, TranVL
    pinocchio.forwardKinematics(model, data, q_desired, qdot_desired, qddot_desired)
    pinocchio.updateFramePlacements(model,data)
    pinocchio.updateGlobalPlacements(model,data)
    pinocchio.crba(model, data, q_desired)
    
    pinocchio.computeJointJacobians(model, data, q_desired)
    pinocchio.computeMinverse(model, data, q_desired)

    LF_tran_cur = data.oMi[LFjoint_id].translation
    RF_tran_cur = data.oMi[RFjoint_id].translation
    LF_rot_cur = data.oMi[LFjoint_id].rotation
    RF_rot_cur = data.oMi[RFjoint_id].rotation

    PELV_tran_cur = np.add(data.oMi[PELVjoint_id].translation, model.inertias[PELVjoint_id].lever)
    PELV_rot_cur = data.oMi[PELVjoint_id].rotation

    pinocchio.crba(model, data, q_desired)
    pinocchio.computeCoriolisMatrix(model, data, q_desired, qdot_desired)
    pinocchio.rnea(model, data, q_desired, qdot_z, qddot_z)
    pinocchio.computeMinverse(model,data,q_desired)
    pinocchio.centerOfMass(model,data,False)
    pinocchio.computeCentroidalMomentum(model,data, q_desired, qdot_desired)
    Jcom = pinocchio.jacobianCenterOfMass(model,data,q_desired,True)
    
    robotCAM = data.hg.angular
    COM_tran_cur = data.com[0]
    RFc_tran_cur = data.oMf[RFcframe_id].translation
    LFc_tran_cur = data.oMf[LFcframe_id].translation

    if contactState == 1:
        contactnum = 2
        robotJac = np.zeros((2 * 6, model.nv))
        robotdJac = np.zeros((2 * 6, model.nv))
    elif contactState == 2:
        contactnum = 1
        robotJac = np.zeros((1 * 6, model.nv))
        robotdJac = np.zeros((1 * 6, model.nv))
    else:
        contactnum = 1
        robotJac = np.zeros((1 * 6, model.nv))
        robotdJac = np.zeros((1 * 6, model.nv))

    M = data.M
    COR = data.C
    G = data.tau
    Minv = data.Minv
    b = np.matmul(COR,qdot_desired)

    LF_j = pinocchio.computeFrameJacobian(model,data,q_desired,LFframe_id,pinocchio.LOCAL_WORLD_ALIGNED)    
    RF_j = pinocchio.computeFrameJacobian(model,data,q_desired,RFframe_id,pinocchio.LOCAL_WORLD_ALIGNED)
    LF_cj = pinocchio.computeFrameJacobian(model,data,q_desired,LFcframe_id,pinocchio.LOCAL_WORLD_ALIGNED)    
    RF_cj = pinocchio.computeFrameJacobian(model,data,q_desired,RFcframe_id,pinocchio.LOCAL_WORLD_ALIGNED)
    RF_cdj = pinocchio.frameJacobianTimeVariation(model,data,q_desired,qdot_desired,RFcframe_id,pinocchio.LOCAL_WORLD_ALIGNED)
    LF_cdj = pinocchio.frameJacobianTimeVariation(model,data,q_desired,qdot_desired,LFcframe_id,pinocchio.LOCAL_WORLD_ALIGNED)

    for i in range(0, contactnum):     
        if i == 0:
            if contactState == 2:
                robotJac[0:6,0:model.nv] = RF_cj
                robotdJac[0:6,0:model.nv] = RF_cdj
            else:
                robotJac[0:6,0:model.nv] = LF_cj
                robotdJac[0:6,0:model.nv] = LF_cdj
        elif i == 1:
            robotJac[6:12,0:model.nv] = RF_cj
            robotdJac[6:12,0:model.nv] = RF_cdj

    robotLambdac = np.linalg.inv(np.matmul(np.matmul(robotJac,Minv),np.transpose(robotJac)))
    robotJcinvT = np.matmul(np.matmul(robotLambdac, robotJac),Minv)
    robotNc = np.subtract(np.identity(model.nv),np.matmul(np.transpose(robotJac),robotJcinvT))
    
    robotW = np.matmul(Minv[6:model.nq-1,0:model.nq],robotNc[0:model.nq-1,6:6+33])
    robotWinv = scipy.linalg.pinv2(robotW)
    
    winvCalc(robotW)
    
    robotPc = np.matmul(robotJcinvT,G)
    robotmuc = np.matmul(robotLambdac,np.subtract(np.matmul(np.matmul(robotJac,Minv),b),np.matmul(robotdJac,qdot_desired)))
    robothc = np.matmul(np.transpose(robotJac),np.add(robotmuc, robotPc))

    virtual_cur = data.oMf[PELVjoint_id].translation
    virtual_dcur = data.ov[PELVjoint_id].linear
    virtual_ddcur = data.oa[PELVjoint_id].linear
    
    TranFV = np.zeros((4,4))
    TranFV[0:3,0:3] = np.identity(3)
    TranFV[0:3,3] = PELV_tran_cur
    TranFV[3,3] = 1.0

    TranFR = np.zeros((4,4))
    TranFR[0:3,0:3] = RF_rot_cur
    TranFR[0:3,3] = RF_tran_cur
    TranFR[3,3] = 1.0

    TranFL = np.zeros((4,4))
    TranFL[0:3,0:3] = LF_rot_cur
    TranFL[0:3,3] = LF_tran_cur
    TranFL[3,3] = 1.0

    TranVR = np.matmul(np.linalg.inv(TranFV),TranFR)
    TranVL = np.matmul(np.linalg.inv(TranFV),TranFL)

def jointUpdateLocal(time):
    global PELV_tran_prev, PELV_tran_dot, V_tran, Vdot_tran, Vddot_tran, V_tranG, TranVRi, TranVLi
    if(time == 0):
        PELV_tran_prev = np.zeros(3)
        PELV_tran_dot = np.zeros(3)
        V_tran = np.zeros(3)
        Vdot_tran = np.zeros(3)
        Vddot_tran = np.zeros(3)
        V_tranG = np.zeros(3)

    q_command[0] = 0.0#com_refx[time]-com_refx[0]
    qdot_command[0] = 0.0
    qddot_command[0] = 0.0
    q_command[1] = 0.0
    qdot_command[1] = 0.0
    qddot_command[1] = 0.0

    if(time != 0):
        if(current_step_num == 0):
            if(foot_step[current_step_num,6] == 1):
                V_tran[0] = (TranVLi[0:3,3] - TranVL[0:3,3])[0]
                V_tran[1] = (TranVLi[0:3,3] - TranVL[0:3,3])[1]
            else:
                V_tran[0] = (TranVRi[0:3,3] - TranVR[0:3,3])[0]
                V_tran[1] = (TranVRi[0:3,3] - TranVR[0:3,3])[1]

            if(time == t_last):
                TranFVi = np.zeros((4,4))
                TranFVi[0:3,0:3] = np.identity(3)
                TranFVi[0:3,3] = PELV_tran_cur
                TranFVi[3,3] = 1.0

                TranFRi = np.zeros((4,4))
                TranFRi[0:3,0:3] = RF_rot_cur
                TranFRi[0:3,3] = RF_tran_cur
                TranFRi[3,3] = 1.0

                TranFLi = np.zeros((4,4))
                TranFLi[0:3,0:3] = LF_rot_cur
                TranFLi[0:3,3] = LF_tran_cur
                TranFLi[3,3] = 1.0

                TranVRi = np.matmul(np.linalg.inv(TranFVi),TranFRi)
                TranVLi = np.matmul(np.linalg.inv(TranFVi),TranFLi)
                V_tranG[0] = V_tran[0]
                V_tranG[1] = V_tran[1]
        else:
            if(foot_step[current_step_num,6] == 1):
                V_tran[0] = (TranVLi[0:3,3] - TranVL[0:3,3])[0] + V_tranG[0]
                V_tran[1] = (TranVLi[0:3,3] - TranVL[0:3,3])[1] + V_tranG[1]
            else:
                V_tran[0] = (TranVRi[0:3,3] - TranVR[0:3,3])[0] + V_tranG[0]
                V_tran[1] = (TranVRi[0:3,3] - TranVR[0:3,3])[1] + V_tranG[1]

            if(time == t_last):
                TranFVi = np.zeros((4,4))
                TranFVi[0:3,0:3] = np.identity(3)
                TranFVi[0:3,3] = PELV_tran_cur
                TranFVi[3,3] = 1.0

                TranFRi = np.zeros((4,4))
                TranFRi[0:3,0:3] = RF_rot_cur
                TranFRi[0:3,3] = RF_tran_cur
                TranFRi[3,3] = 1.0

                TranFLi = np.zeros((4,4))
                TranFLi[0:3,0:3] = LF_rot_cur
                TranFLi[0:3,3] = LF_tran_cur
                TranFLi[3,3] = 1.0

                TranVRi = np.matmul(np.linalg.inv(TranFVi),TranFRi)
                TranVLi = np.matmul(np.linalg.inv(TranFVi),TranFLi)
                V_tranG[0] = V_tran[0]
                V_tranG[1] = V_tran[1]
        
        Vdot_tran[0] = (V_tran[0] - PELV_tran_prev[0]) * float(hz)#com_refdx[time]#virtual_dcur[0]#np.subtract(PELV_tran_cur[0], PELV_tran_prev[0]) *  float(hz)
        Vddot_tran[0] = (Vdot_tran[0] - PELV_tran_dot[0]) * float(hz)#com_refddx[time]#virtual_ddcur[0]#np.subtract(qdot_command[0], PELV_tran_dot[0]) * float(hz)
        Vdot_tran[1] = (V_tran[1] - PELV_tran_prev[1]) * float(hz)#com_refdy[time]#virtual_dcur[1]#np.subtract(PELV_tran_cur[1], PELV_tran_prev[1]) *  float(hz)
        Vddot_tran[1] = (Vdot_tran[1] - PELV_tran_dot[1]) * float(hz)#com_refddy[time]#virtual_ddcur[1]#np.subtract(qdot_command[1], PELV_tran_dot[1]) * float(hz)
        
    if(time == 0):
        leg_qs[time,:] = q_init[7:19]
        PELV_tran_prev[0] = V_tran[0]
        PELV_tran_dot[0] = Vdot_tran[0]
        PELV_tran_prev[1] = V_tran[1]
        PELV_tran_dot[1] = Vdot_tran[1]
    else:
        PELV_tran_prev[1] = V_tran[1]
        PELV_tran_dot[1] = Vdot_tran[1]
        PELV_tran_prev[1] = V_tran[1]
        PELV_tran_dot[1] = Vdot_tran[1]

    for i in range(7, 19):
        q_command[i] = leg_qs[time, i-7]
        
    for i in range(19,40):
        q_command[i] = q_init[i]

    for i in range(6, 18):
        qdot_command[i] = leg_qdots[time, i-6]
        qddot_command[i] = leg_qddots[time, i-6]

    for i in range(18,39):
        qdot_command[i] = 0.0
        qddot_command[i] = 0.0

def jointUpdateGlobal(time):
    global PELV_tran_prev, PELV_tran_dot
    if(time == 0):
        PELV_tran_prev = np.zeros(3)
        PELV_tran_dot = np.zeros(3)

    q_command[0] = 0.0#TranVRi[0,3] - TranVR[0,3]
    qdot_command[0] = 0.0
    qddot_command[0] = 0.0
    q_command[1] = 0.0#TranVRi[1,3] - TranVR[1,3]
    qdot_command[1] = 0.0
    qddot_command[1] = 0.0
    
    if(time != 0):
        q_command[0] = V_tran[0]
        q_command[1] = V_tran[1]

        #qdot_command[0] = (q_command[0] - PELV_tran_prev[0]) * float(hz)#com_refdx[time]#virtual_dcur[0]#np.subtract(PELV_tran_cur[0], PELV_tran_prev[0]) *  float(hz)
        #qddot_command[0] = (qdot_command[0] - PELV_tran_dot[0]) * float(hz)#com_refddx[time]#virtual_ddcur[0]#np.subtract(qdot_command[0], PELV_tran_dot[0]) * float(hz)
        #qdot_command[1] = (q_command[1] - PELV_tran_prev[1]) * float(hz)#com_refdy[time]#virtual_dcur[1]#np.subtract(PELV_tran_cur[1], PELV_tran_prev[1]) *  float(hz)
        #qddot_command[1] = (qdot_command[1] - PELV_tran_dot[1]) * float(hz)#com_refddy[time]#virtual_ddcur[1]#np.subtract(qdot_command[1], PELV_tran_dot[1]) * float(hz)
        
    if(time == 0):
        leg_qs[time,:] = q_init[7:19]
        PELV_tran_prev[0] = q_command[0]
        PELV_tran_dot[0] = qdot_command[0]
        PELV_tran_prev[1] = q_command[1]
        PELV_tran_dot[1] = qdot_command[1]
    else:
        PELV_tran_prev[1] = q_command[1]
        PELV_tran_dot[1] = qdot_command[1]
        PELV_tran_prev[1] = q_command[1]
        PELV_tran_dot[1] = qdot_command[1]

    for i in range(7, 19):
        q_command[i] = leg_qs[time, i-7]
        
    for i in range(19,40):
        q_command[i] = q_init[i]

    for i in range(6, 18):
        qdot_command[i] = leg_qdots[time, i-6]
        qddot_command[i] = leg_qddots[time, i-6]

    for i in range(18,39):
        qdot_command[i] = 0.0
        qddot_command[i] = 0.0
    

def phaseUpdate(time):
    global t_last, t_start_real, t_start, phaseChange, phaseChange1, double2Single_pre, double2Single, single2Double_pre, single2Double, current_step_num
    if(time == 0):
        current_step_num = 0

    if (time == t_last):
        if (current_step_num != foot_step_number - 1):
            t_start = t_last + 1
            t_start_real = t_start + t_rest_1
            t_last = t_start + t_total - 1
            current_step_num =current_step_num + 1
            if (current_step_num == foot_step_number):
                current_step_num = foot_step_number - 1
 
    if (time >= t_start_real + t_double_1 *0.25 and time <= t_start_real + t_double_1 and current_step_num != 0):
        phaseChange = True
        phaseChange1 = False
        double2Single_pre = t_start_real + t_double_1 * 0.25
        double2Single = t_start_real + t_double_1
    else:
        phaseChange = False

    if (time >= t_start + t_total - t_rest_2 - t_double_2 + 1 and time <= t_start + t_total - t_rest_2 - t_double_2 + t_double_1 - 1 and current_step_num != 0 and phaseChange == False):
        phaseChange1 = True
        phaseChange = False
        single2Double_pre = t_start + t_total - t_rest_2 - t_double_2 + 1
        single2Double = t_start + t_total - t_rest_2 - t_double_2 + t_double_1
    else:
        phaseChange1 = False
        if (time < t_start_real + t_double_1 * 0.25 and time > t_start_real + t_double_1):
            phaseChange = False

def estimateContactForce(time, qddot_desired):
    global robotContactForce, robotTorque, rate, TorqueContact, stateCo, f1, f2, phaseChange, phaseChange1, double2Single_pre, double2Single
    '''
    robotContactForce = pinocchio.utils.zero(12)
    robotTorque = np.matmul(scipy.linalg.pinv2(robotNc),np.subtract(np.add(np.add(np.matmul(M,qddot_desired),b),G),robothc))
    robotTorque_r = np.matmul(robotWinv,np.matmul(Minv[6:model.nq],np.matmul(robotNc, robotTorque)))
    G_r = np.matmul(robotWinv,np.matmul(Minv[6:model.nq],np.matmul(robotNc,G)))

    stateCo =0
    if (phaseChange == True) and (phaseChange1 == False):
        stateCo = 1
        rate = quinticSpline(time, double2Single_pre, double2Single, 1, 0, 0, 0, 0, 0)
        TorqueContact = contactRedistributionWalking(G_r, 0.0, rate, foot_step[current_step_num, 6])                
    elif (phaseChange == False) and (phaseChange1 == True):
        rate = quinticSpline(time, single2Double_pre, single2Double, 0, 0, 0, 1, 0, 0)
        stateCo =2 
        if (current_step_num < foot_step_number - 1):
            TorqueContact = contactRedistributionWalking(G_r, 0.0, rate, foot_step[current_step_num, 6]) 
        else:
            TorqueContact = contactRedistributionWalking(G_r, 0.0, rate, foot_step[foot_step_number - 1, 6])  
    else:
        rate = 1.0
        stateCo =3
        if (current_step_num < foot_step_number - 1):
            TorqueContact = contactRedistributionWalking(G_r, 0.0, rate, foot_step[current_step_num + 1, 6])
        else:
            TorqueContact = contactRedistributionWalking(G_r, 0.0, rate, foot_step[foot_step_number - 2, 6])

    robotTorque[6:model.nq] = robotTorque_r + TorqueContact
    robotTorque_r = robotTorque_r + TorqueContact

    if contactState == 1:
        robotContactForce = np.subtract(np.subtract(np.matmul(robotJcinvT[0:12,6:6+model.nq],robotTorque_r), robotPc),robotmuc)
    elif contactState == 2:
        robotContactForce[0:6] = np.zeros(6)
        robotContactForce[6:12] = np.subtract(np.subtract(np.matmul(robotJcinvT[0:6,6:6+model.nq],robotTorque_r), robotPc),robotmuc)
    else:   
        robotContactForce[0:6] = np.subtract(np.subtract(np.matmul(robotJcinvT[0:6,6:6+model.nq],robotTorque_r), robotPc),robotmuc)
        robotContactForce[6:12] = np.zeros(6)  
    
    contact_dof_ = int(np.shape(robotJac)[0])
    
    if(time == 1501 or time == 1500 or time == 4000 or time == 5001 or time == 1001):
        print(time)
    else:
    '''

    robotContactForce = pinocchio.utils.zero(12)
    robotTorque = np.matmul(scipy.linalg.pinv2(robotNc),np.subtract(np.add(np.add(np.matmul(M,qddot_desired),b),G),robothc))
    
    if (phaseChange == True and phaseChange1 == False):
        rate = quinticSpline(time, double2Single_pre, double2Single, 1, 0, 0, 0, 0, 0)
        TorqueContact = contactRedistributionWalking(robotTorque, 0.0, rate, foot_step[current_step_num, 6])                
    elif (phaseChange == False and phaseChange1 == True):
        rate = quinticSpline(time, single2Double_pre, single2Double, 0, 0, 0, 1, 0, 0)
        if (current_step_num < foot_step_number - 1):
            TorqueContact = contactRedistributionWalking(robotTorque, 0.0, rate, foot_step[current_step_num, 6]) 
        else:
            TorqueContact = contactRedistributionWalking(robotTorque, 0.0, rate, foot_step[foot_step_number - 1, 6])  
    else:
        rate = 1.0
        if (current_step_num < foot_step_number - 1):
            TorqueContact = contactRedistributionWalking(robotTorque, 0.0, rate, foot_step[current_step_num + 1, 6])
        else:
            TorqueContact = contactRedistributionWalking(robotTorque, 0.0, rate, foot_step[foot_step_number - 2, 6])
        if(current_step_num != 0):
            TorqueContact = np.zeros(33)
    
    robotTorque[6:model.nq] = robotTorque[6:model.nq]# + TorqueContact
    
    if contactState == 1:
        robotContactForce = np.subtract(np.subtract(np.matmul(robotJcinvT,robotTorque), robotPc),robotmuc)
    elif contactState == 2:
        robotContactForce[0:6] = np.zeros(6)
        robotContactForce[6:12] = np.subtract(np.subtract(np.matmul(robotJcinvT,robotTorque), robotPc),robotmuc)
    else:   
        robotContactForce[0:6] = np.subtract(np.subtract(np.matmul(robotJcinvT,robotTorque), robotPc),robotmuc)
        robotContactForce[6:12] = np.zeros(6)  

    #print(np.shape(robotdJac))
    #robotPc1 = np.matmul(robotJcinvT[0:6,12:18],G[12:18])
    #robotmuc1 = np.matmul(robotLambdac,np.subtract(np.matmul(np.matmul(robotJac[0:6,12:18],Minv[12:18,12:18]),b[12:18]),np.matmul(robotdJac[0:6,12:18],qdot_command[12:18])))
 
    #robotContactForce1 = np.subtract(np.subtract(np.matmul(robotJcinvT[0:6,12:18],robotTorque[12:18]), robotPc1),robotmuc1)
    #robotContactForce1 = np.matmul(robotJcinvT[0:12,12:18],robotTorque[12:18])
    #f2.write('%f %f %f %f %f %f %f %f' % (time, contactState, robotContactForce1[0], robotContactForce1[1], robotContactForce1[2], robotContactForce1[3], robotContactForce1[4], robotContactForce1[5]))
    #f2.write("\n")
       
def calMargin():
    global contact_margin
    #print("CalMargin")
    contact_margin = np.zeros(10)
    mu_s = 0.3
    lx = 0.26
    ly = 0.1

    if contactState == 1:
        contact_margin[0] = robotContactForce[2]
        contact_margin[1] = np.sqrt(robotContactForce[0] * robotContactForce[0] + robotContactForce[1] * robotContactForce[1]) - mu_s * abs(robotContactForce[2])
        contact_margin[2] = np.abs(robotContactForce[5]) - mu_s * np.abs(robotContactForce[2])
        contact_margin[3] = np.abs(robotContactForce[3]) - ly/2 * np.abs(robotContactForce[2])
        contact_margin[4] = np.abs(robotContactForce[4]) - lx/2 * np.abs(robotContactForce[2])
        contact_margin[5] = robotContactForce[8]
        contact_margin[6] = np.sqrt(robotContactForce[6] * robotContactForce[6] + robotContactForce[7] * robotContactForce[7]) - mu_s * abs(robotContactForce[8])
        contact_margin[7] = np.abs(robotContactForce[11]) - mu_s * np.abs(robotContactForce[8])
        contact_margin[8] = np.abs(robotContactForce[9]) - ly/2 * np.abs(robotContactForce[8])
        contact_margin[9] = np.abs(robotContactForce[10]) - lx/2 * np.abs(robotContactForce[8])
    elif contactState == 2:
        contact_margin[5] = robotContactForce[8]
        contact_margin[6] = np.sqrt(robotContactForce[6] * robotContactForce[6] + robotContactForce[7] * robotContactForce[7]) - mu_s * abs(robotContactForce[8])
        contact_margin[7] = np.abs(robotContactForce[11]) - mu_s * np.abs(robotContactForce[8])
        contact_margin[8] = np.abs(robotContactForce[9]) - ly/2 * np.abs(robotContactForce[8])
        contact_margin[9] = np.abs(robotContactForce[10]) - lx/2 * np.abs(robotContactForce[8])
    else:   
        contact_margin[0] = robotContactForce[2]
        contact_margin[1] = np.sqrt(robotContactForce[0] * robotContactForce[0] + robotContactForce[1] * robotContactForce[1]) - mu_s * abs(robotContactForce[2])
        contact_margin[2] = np.abs(robotContactForce[5]) - mu_s * np.abs(robotContactForce[2])
        contact_margin[3] = np.abs(robotContactForce[3]) - ly/2 * np.abs(robotContactForce[2])
        contact_margin[4] = np.abs(robotContactForce[4]) - lx/2 * np.abs(robotContactForce[2])      

def talker():
    modelInitialize()
    walkingSetup()
    footStep()
    cpGenerator()
    comGenerator()
    swingFootGenerator()

    global contactState, f1, f2

    f = open("newfile.txt", 'w')
    f1 = open("newfile1.txt", 'w')
    f2 = open("newfile2.txt", 'w')

    global phaseChange, phaseChange1
    for i in range(0, int(total_tick)):
        LF_tran[0] = lfoot[i,0]
        LF_tran[1] = lfoot[i,1]
        LF_tran[2] = lfoot[i,2]
        RF_tran[0] = rfoot[i,0]
        RF_tran[1] = rfoot[i,1]
        RF_tran[2] = rfoot[i,2]
        PELV_tran[0] = com_refx[i]
        PELV_tran[1] = com_refy[i]
        PELV_tran[2] = PELV_tran_init[2]

        inverseKinematics(i, LF_rot, RF_rot, PELV_rot, LF_tran, RF_tran, PELV_tran, HRR_tran_init, HLR_tran_init, HRR_rot_init, HLR_rot_init, PELV_tran_init, PELV_rot_init, CPELV_tran_init)
       
    for i in range(0, int(total_tick)): #range(int(t_total + t_temp + 2), int(t_totfal * (foot_step_number - 2) + t_temp)): 
        contactState = phase_variable[i]
        jointUpdateLocal(i)
        modelUpdate(q_command,qdot_command,qddot_command)
        comJacobianinverseKinematics(i, LF_rot, RF_rot, PELV_rot, LF_tran, RF_tran, PELV_tran, HRR_tran_init, HLR_tran_init, HRR_rot_init, HLR_rot_init, PELV_tran_init, PELV_rot_init, CPELV_tran_init)
        
        f2.write('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' % (i, contactState, com_refx[i], com_refy[i], PELV_tran_cur[0], RF_tran_cur[0], V_tran[0], V_tran[1], rfoot[i,2], lfoot[i,2], LF_tran_cur[0], LF_tran_cur[2],leg_qs[i,0], leg_qs[i,1], leg_qs[i,2], leg_qs[i,3], leg_qs[i,4], leg_qs[i,5]))
        f2.write("\n")
    
        jointUpdateGlobal(i)
        modelUpdate(q_command,qdot_command,qddot_command)
        '''
        print(q_command)
        print(RF_tran_cur - TranVR[0:3,3])
        print(TranVR[0:3,3])
        print(TranVRi[0:3,3])
        print(t_last)
        '''
        if(i == 0):
            robotCAMI = np.zeros(3)

        phaseUpdate(i)
        #if(i>2999 and i <4001):
        estimateContactForce(i, qddot_command)
        calMargin()  
        robotCAMI = np.add(robotCAMI, robotCAM * 0.002)
        f.write('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ' % (i, current_step_num, com_refx[i],zmp_refx[i], com_refy[i],zmp_refy[i],com_refdy[i], com_refddy[i],leg_qs[i,0], leg_qs[i,1], leg_qs[i,2], leg_qs[i,3], leg_qs[i,4], leg_qs[i,5], leg_qdots[i,0], leg_qdots[i,1], leg_qdots[i,2], leg_qdots[i,3], leg_qdots[i,4], leg_qdots[i,5], leg_qddots[i,0], leg_qddots[i,1], leg_qddots[i,2], leg_qddots[i,3], leg_qddots[i,4], leg_qddots[i,5]))
        f.write("\n")
        #f2.write('%f %f %f %f %f %f %f %f %f %f %f %f %f' % (i, robotCAMI[0], robotCAMI[1], robotCAMI[2],robotCAM[0], robotCAM[1], robotCAM[2], robotContactForce[0], robotContactForce[1], robotContactForce[2], robotContactForce[3], robotContactForce[4], robotContactForce[5]))            
        #f2.write("\n")
        f1.write('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' % (i, com_refx[i], com_refy[i], COM_tran_cur[0], COM_tran_cur[1], V_tran[0] - V_tranG[0], V_tranG[0], lfoot[i,0], lfoot[i,2], LF_tran_cur[0], LF_tran_cur[2],leg_qs[i,0], leg_qs[i,1], leg_qs[i,2], leg_qs[i,3], leg_qs[i,4], leg_qs[i,5]))
        f1.write("\n")
    
    f.close()
    f1.close()
    f2.close()

if __name__=='__main__':
    talker()
