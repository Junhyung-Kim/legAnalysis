#!/usr/bin/env python 
from __future__ import print_function
import pinocchio
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from sys import argv
from scipy.interpolate import CubicSpline
from os.path import dirname, join, abspath

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
    global x_direction, y_direction, yaw_direction, step_length, hz, total_tick, t_total_t, t_start_real, t_temp_t, t_double, t_rest_1, t_rest_2, t_start, t_total, t_temp, t_last, t_double_1, t_double_2
    global zc, wn, current_step_num, ref_zmp, ref_com, walking_tick, total_tick, phase_variable, lfoot, rfoot, foot_height, foot_step_dir
    hz = 100
    x_direction = 1.00
    y_direction = 0.00
    yaw_direction = 0.00
    step_length = 0.10
    
    t_total_t = 1.0
    t_temp_t = 2.0
    t_double = 0.1
    
    t_total = t_total_t * hz
    t_temp = t_temp_t * hz
    t_start = t_temp + 1
    t_last = t_temp + t_total
    t_double_1 = 0.1 * hz
    t_double_2 = 0.1 * hz
    t_rest_1 = 0.1 * hz
    t_rest_2 = 0.1 * hz
    t_start_real = t_start + t_rest_1

    foot_step_dir = 1

    foot_height = 0.05

    current_step_num = 0
    zc = 0.727822
    wn = np.sqrt(9.81/zc)

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
                numberOfFootstep = numberOfFootstep + 1

    if (middle_foot_step_number != 0 or np.abs(middle_residual_length) >= 0.0001):
        if (middle_foot_step_number % 2 == 0):
            numberOfFootstep = numberOfFootstep + 2
        else:
            if (np.abs(middle_residual_length) >= 0.0001):
                numberOfFootstep = numberOfFootstep + 3
            else:
                numberOfFootstep = numberOfFootstep + 1

    if (final_foot_step_number != 0 or abs(final_residual_angle) >= 0.0001):
        if (abs(final_residual_angle) >= 0.0001):
            numberOfFootstep = numberOfFootstep + 2
        else:
            numberOfFootstep = numberOfFootstep + 1

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
    global zmp_refx, zmp_refy, com_refx, com_refy, walking_tick, capturePoint_refx, capturePoint_refy, com_refdx, com_refdy, com_refddx, com_refddy, total_tick
    capturePoint_ox = np.zeros(foot_step_number + 3)
    capturePoint_oy = np.zeros(foot_step_number + 3)
    zmp_dx = np.zeros(foot_step_number + 2)
    zmp_dy = np.zeros(foot_step_number + 2)
    b_offset = np.zeros(foot_step_number + 2)
    capturePoint_offsetx = np.zeros(foot_step_number + 3)
    capturePoint_offsety = np.zeros(foot_step_number + 3)

    total_tick = t_total * (foot_step_number + 1) + t_temp - 1

    walking_tick = np.zeros(int(total_tick))

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
        walking_tick[i] = i
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
            zmp_refx[0] = 0.0
            zmp_refy[0] = 0.0
        else:
            zmp_refx[i] = (capturePoint_refx[i - 1]) - (capturePoint_refx[i] - capturePoint_refx[i - 1]) * hz / (wn)
            zmp_refy[i] = (capturePoint_refy[i - 1]) - (capturePoint_refy[i] - capturePoint_refy[i - 1]) * hz / (wn)

def comGenerator():
    for i in range(0, int(t_total * (foot_step_number + 1) + t_temp - 1)):
        if (i >= t_temp - t_total):
            com_refx[i] = wn / hz * capturePoint_refx[i] + (1 - wn / hz) * com_refx[i - 1]
            com_refy[i] = wn / hz * capturePoint_refy[i] + (1 - wn / hz) * com_refy[i - 1]
            com_refdx[i] = (com_refx[i] - com_refx[i - 1]) * hz
            com_refdy[i] = (com_refy[i] - com_refy[i - 1]) * hz
            com_refddx[i] = (com_refdx[i] - com_refdx[i - 1]) * hz
            com_refddy[i] = (com_refdy[i] - com_refdy[i - 1]) * hz
        else:
            com_refx[i] = PELV_tran_init[0]
            com_refy[i] = PELV_tran_init[1]
            com_refdx[i] = 0.0
            com_refdy[i] = 0.0
            com_refddx[i] = 0.0
            com_refddy[i] = 0.0
                

def swingFootGenerator():
    global lfoot, rfoot, phase_variable
    phase_variable = np.zeros(int(total_tick))
    lfoot = np.zeros((int(total_tick),3))
    rfoot = np.zeros((int(total_tick),3))

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
                        if (i < t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0):
                            rfoot[i,2] = quinticSpline(i, t_start_real + t_total + t_double_1, t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2) / 2, RF_tran[2], 0.0, 0.0, RF_tran[2] + foot_height, 0.0, 0.0)               
                        else:
                            rfoot[i,2] = quinticSpline(i, t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0, t_start + t_total + t_total - t_rest_2 - t_double_2, RF_tran[2] + foot_height, 0.0, 0.0, RF_tran[2], 0.0, 0.0)
                        if(i >= t_start_real + t_total + t_double_1) and ( i <= t_start + t_total + t_total - t_rest_2 - t_double_2):
                            phase_variable[i] = 3
                    else:
                        rfoot[i,0] = foot_step[j - 1, 0]
                        rfoot[i,1] = foot_step[j - 1, 1]
                        rfoot[i,2] = RF_tran[2]

                        lfoot[i,1] = LF_tran[1]
                        lfoot[i,0] = quinticSpline(i, t_start_real + t_total + t_double_1, t_start + t_total * 2 - t_rest_2 - t_double_2, LF_tran[0], 0.0, 0.0, foot_step[j, 0], 0.0, 0.0)
                        
                        if (i < t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0):
                            lfoot[i,2] = quinticSpline(i, t_start_real + t_total + t_double_1, t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2, LF_tran[2], 0.0, 0.0, LF_tran[2] + foot_height, 0.0, 0.0)
                        else:
                            lfoot[i,2] = quinticSpline(i, t_start_real + t_total + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0, t_start + t_total + t_total - t_rest_2 - t_double_2, LF_tran[2] + foot_height, 0.0, 0.0, LF_tran[2], 0.0, 0.0)
            
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
                        if (i < t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0):
                            rfoot[i,2] = quinticSpline(i, t_start_real + t_total * j + t_double_1, t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2, LF_tran[2], 0.0, 0.0, RF_tran[2] + foot_height, 0.0, 0.0)
                        else:
                            rfoot[i,2] = quinticSpline(i, t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0, t_start + t_total * j + t_total - t_rest_2 - t_double_2, RF_tran[2] + foot_height, 0.0, 0.0, RF_tran[2], 0.0, 0.0)
                        if(i >= t_start_real + t_total * j + t_double_1) and ( i <=  t_start + t_total * j + t_total - t_rest_2 - t_double_2):
                            phase_variable[i] = 3
                    else:
                        lfoot[i,1] = foot_step[j - 2, 1]
                        rfoot[i,1] = foot_step[j - 1, 1]
                        rfoot[i,0] = foot_step[j - 1, 0]
                        rfoot[i,2] = RF_tran[2]
                        lfoot[i,0] = quinticSpline(i, t_start_real + t_total * j + t_double_1, t_start + t_total * (j + 1) - t_rest_2 - t_double_2, foot_step[j - 2, 0], 0.0, 0.0, foot_step[j, 0], 0.0, 0.0)
                        if (i < t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0):
                            lfoot[i,2] = quinticSpline(i, t_start_real + t_total * j + t_double_1, t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2, LF_tran[2], 0.0, 0.0, LF_tran[2] + foot_height, 0.0, 0.0)
                        else:
                            lfoot[i,2] = quinticSpline(i, t_start_real + t_total * j + t_double_1 + (t_total - t_rest_1 - t_rest_2 - t_double_1 - t_double_2 ) / 2.0, t_start + t_total * j + t_total - t_rest_2 - t_double_2, LF_tran[2] + foot_height, 0.0, 0.0, LF_tran[2], 0.0, 0.0)
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
                           
def inverseKinematics(time, LF_rot_c, RF_rot_c, PELV_rot_c, LF_tran_c, RF_tran_c, PELV_tran_c, HRR_tran_init_c, HLR_tran_init_c, HRR_rot_init_c, HLR_rot_init_c, PELV_tran_init_c, PELV_rot_init_c, CPELV_tran_init_c):
    global leg_q, leg_qdot, leg_qddot
    M_PI = 3.14159265358979323846
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
    lp = np.matmul(np.transpose(LF_rot_c), np.transpose(lpt))
    rp = np.matmul(np.transpose(RF_rot_c), np.transpose(rpt))
    
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

    r_tl2 = np.matmul(np.matmul(np.matmul(PELV_rot, np.transpose(r_l4l5)),np.transpose(r_l3l4)),np.transpose(r_l2l3))
    leg_q[1] = np.arcsin(r_tl2[2, 1])

    c_lq5 = np.divide(-r_tl2[0, 1], np.cos(leg_q[1]))

    if c_lq5 > 1.0:
        c_lq5 = 1.0
    elif c_lq5 < -1.0:
        c_lq5 = -1.0
    else:
        c_lq5 = c_lq5

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

    r_tr2 = np.matmul(np.matmul(np.matmul(PELV_rot, np.transpose(r_r4r5)),np.transpose(r_r3r4)),np.transpose(r_r2r3))
    leg_q[7] = np.arcsin(r_tr2[2,1])
    c_rq5 = -r_tr2[0, 1] / np.cos(leg_q[7])

    if c_rq5 > 1.0:
        c_rq5 = 1.0
    elif c_rq5 < -1.0:
        c_rq5 = -1.0
    else:
        c_rq5 = c_rq5  
    
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
        leg_qddots[time,:] =np.subtract(leg_qdots[time,:], leg_qdots[time-1,:]) * hz

    leg_qdot = leg_qdots[time,:]
    leg_qddot = leg_qddots[time,:]
        
def modelInitialize():
    global model, foot_distance, data, LFframe_id, RFframe_id, PELVjoint_id, LHjoint_id, RHjoint_id, LFjoint_id, RFjoint_id, LFcframe_id, RFcframe_id, q, qdot, qddot, LF_tran, RF_tran, PELV_tran, LF_rot, RF_rot, PELV_rot, qdot_z, qddot_z, HRR_rot_init, HLR_rot_init, HRR_tran_init, HLR_tran_init, PELV_tran_init, PELV_rot_init, CPELV_tran_init, q_command, qdot_command, qddot_command
    model = pinocchio.buildModelFromUrdf("/home/jhk/legAnalysis/dyros_tocabi_with_redhands.urdf",pinocchio.JointModelFreeFlyer())      
   
    LFframe_id = model.getFrameId("L_Foot_Link")
    RFframe_id = model.getFrameId("R_Foot_Link")

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
    q_init = [0.0, 0.0, 0.814175, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.56, 1.3, -0.73, 0.0, 0.0, 0.0, -0.56, 1.3, -0.73, 0.0, 0.0, 0.0, 0.0, 0.2, 0.6, 1.5, -1.47, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -0.2, -0.6, -1.5, 1.47, 1.0, 0.0, 1.0, 0.0]

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

    LF_tran = data.oMi[LFjoint_id].translation
    RF_tran = data.oMi[RFjoint_id].translation
    LF_rot = data.oMi[LFjoint_id].rotation
    RF_rot = data.oMi[RFjoint_id].rotation

    PELV_tran = np.add(data.oMi[PELVjoint_id].translation, model.inertias[PELVjoint_id].lever)
    PELV_rot = data.oMi[PELVjoint_id].rotation

    LF_tran_init = data.oMi[LFjoint_id].translation
    RF_tran_init = data.oMi[RFjoint_id].translation
    HLR_tran_init = data.oMi[LHjoint_id].translation
    HRR_tran_init = data.oMi[RHjoint_id].translation
    LF_rot_init = data.oMi[LFjoint_id].rotation
    RF_rot_init = data.oMi[RFjoint_id].rotation
    HLR_rot_init = data.oMi[LHjoint_id].rotation
    HRR_rot_init = data.oMi[RHjoint_id].rotation

    PELV_tran_init = np.add(data.oMi[PELVjoint_id].translation, model.inertias[PELVjoint_id].lever)
    CPELV_tran_init = data.oMi[PELVjoint_id].translation 
    PELV_rot_init = data.oMi[PELVjoint_id].rotation

    foot_distance = LF_tran_init - RF_tran_init

def modelUpdate(q_desired, qdot_desired, qddot_desired):
    global contactnum, M, G, COR, Minv, b, robotJac, robotdJac, robotIc, LF_j, RF_j, LF_cj, RF_cj, LF_cdj, RF_cdj, robotLambdac, robotJcinvT, robotNc, robotPc, robotmuc, robothc, LF_tran_cur, RF_tran_cur, PELV_tran_cur
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

    contactnum = 0

    if contactState == 1:
        contactnum = 2
        robotJac = np.zeros((2 * 6, model.nv))
        robotdJac = np.zeros((2 * 6, model.nv))
        robotIc = np.zeros((2 * 6, 2 * 6))
    elif contactState == 2:
        contactnum = 1
        robotJac = np.zeros((1 * 6, model.nv))
        robotdJac = np.zeros((1 * 6, model.nv))
        robotIc = np.zeros((1 * 6, 1 * 6))
    else:
        contactnum = 1
        robotJac = np.zeros((1 * 6, model.nv))
        robotdJac = np.zeros((1 * 6, model.nv))
        robotIc = np.zeros((1 * 6, 1 * 6))

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
            robotJac[0:6,0:model.nv] = LF_cj
            robotdJac[0:6,0:model.nv] = LF_cdj
        elif i == 1:
            robotJac[6:12,0:model.nv] = RF_cj
            robotdJac[6:12,0:model.nv] = RF_cdj

    robotLambdac = np.linalg.inv(np.matmul(np.matmul(robotJac,Minv),np.transpose(robotJac)))
    robotJcinvT = np.matmul(np.matmul(robotLambdac, robotJac),Minv)
    robotNc = np.subtract(np.identity(model.nv),np.matmul(np.transpose(robotJac),robotJcinvT))
    robotPc = np.matmul(robotJcinvT,G)

    robotmuc = np.matmul(robotLambdac,np.subtract(np.matmul(np.matmul(robotJac,Minv),b),np.matmul(robotdJac,qdot_desired)))
    robothc = np.matmul(np.transpose(robotJac),np.add(robotmuc, robotPc))

def jointUpdate(time):
    q_command[0] = com_refx[time]
    qdot_command[0] = com_refdx[time]
    qddot_command[0] = com_refddx[time]
    q_command[1] = com_refy[time]
    qdot_command[1] = com_refdy[time]
    qddot_command[1] = com_refddy[time]

    for i in range(7, 19):
        q_command[i] = leg_q[i-7]

    for i in range(6, 18):
        qdot_command[i] = leg_qdot[i-6]
        qddot_command[i] = leg_qddot[i-6]

    print(qdot_command)

    print("a")
    print(leg_qdot)

def estimateContactForce(qddot_desired):
    global robotContactForce
    robotContactForce = pinocchio.utils.zero(12)

    robotTorque = np.matmul(np.linalg.pinv(robotNc),np.subtract(np.add(np.add(np.matmul(M,qddot_desired),b),G),robothc))

    if contactState == 1:
        robotContactForce = np.subtract(np.subtract(np.matmul(robotJcinvT,robotTorque), robotPc),robotmuc)
    elif contactState == 2:
        robotContactForce[0:6] = np.zeros(6)
        robotContactForce[6:12] = np.subtract(np.subtract(np.matmul(robotJcinvT,robotTorque), robotPc),robotmuc)
    else:   
        robotContactForce[0:6] = np.subtract(np.subtract(np.matmul(robotJcinvT,robotTorque), robotPc),robotmuc)
        robotContactForce[6:12] = np.zeros(6)
    #print(robotTorque)
    #print(robotContactForce)

def calMargin():
    print("CalMargin")
    contact_margin = np.zeros(10)
    mu_s = 0.3
    lx = 1.0
    ly = 0.5

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

    global contactState

    f = open("newfile.txt", 'w')
    f1 = open("newfile1.txt", 'w')
    
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
        
        contactState = phase_variable[i]
        inverseKinematics(i, LF_rot, RF_rot, PELV_rot, LF_tran, RF_tran, PELV_tran, HRR_tran_init, HLR_tran_init, HRR_rot_init, HLR_rot_init, PELV_tran_init, PELV_rot_init, CPELV_tran_init)
        jointUpdate(i)
        modelUpdate(q_command,qdot_command,qddot_command)
        estimateContactForce(qddot_command)
        print(i)
        
        #if(np.abs(robotContactForce[2])<2000):
        f.write('%f %f %f' % (robotContactForce[2], contactState, com_refy[i]))
        f.write("\n")
        
        #calMargin()
    
    f.close()
    f1.close()

if __name__=='__main__':
    talker()
