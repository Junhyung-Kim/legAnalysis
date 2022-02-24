%load humanoid model
robot = importrobot("dyros_tocabi_sim.urdf");
config = homeConfiguration(robot);
q_init = [0.0; 0.0; -0.55; 1.26; -0.71; 0.0; 0.0; 0.0; -0.55; 1.26; -0.71; 0.0; 0.0; 0.0; 0.0; 0.2; 0.6; 1.5; -1.47; -1.0; 0.0; -1.0; 0.0; 0.0; 0.0; -0.2; -0.6; -1.5; 1.47; 1.0; 0.0; 1.0; 0.0];
jointnum = size(q_init,1);
for i = 1:jointnum
    config(i).JointPosition = q_init(i,1);
end

PELV = getTransform(robot,config,'Pelvis_Link','Pelvis_Link');
RF = getTransform(robot,config,'R_Foot_Link','Pelvis_Link');
LF = getTransform(robot,config,'L_Foot_Link','Pelvis_Link');




