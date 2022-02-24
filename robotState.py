#!/usr/bin/env python 
import pinocchio
import numpy as np
from sys import argv
from os.path import dirname, join, abspath


def talker():
    model = pinocchio.buildModelFromUrdf("/home/jhk/legAnalysis/dyros_tocabi_sim.urdf",pinocchio.JointModelFreeFlyer())      
    data = model.createData()

    q = pinocchio.randomConfiguration(model)
    qdot = pinocchio.utils.zero(model.nv)
    qdot_t = pinocchio.utils.zero(model.nv)
    qddot_t = pinocchio.utils.zero(model.nv)


if __name__=='__main__':
    talker()