%load humanoid model
robot = importrobot("dyros_tocabi_sim.urdf");
config = homeConfiguration(robot);
q_init = [0.0; 0.0; -0.55; 1.26; -0.71; 0.0; 0.0; 0.0; -0.55; 1.26; -0.71; 0.0; 0.0; 0.0; 0.0; 0.2; 0.6; 1.5; -1.47; -1.0; 0.0; -1.0; 0.0; 0.0; 0.0; -0.2; -0.6; -1.5; 1.47; 1.0; 0.0; 1.0; 0.0];
qd_init = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0];
jointnum = size(q_init,1);
for i = 1:jointnum
    config(i).JointPosition = q_init(i,1);
end

PELV = getTransform(robot,config,'Pelvis_Link','Pelvis_Link');
T_RF = getTransform(robot,config,'R_Foot_Link','Pelvis_Link');
T_LF = getTransform(robot,config,'L_Foot_Link','Pelvis_Link');

robot.DataFormat = 'column'; 
robot.Gravity = [0 0 -9.81];

G = gravityTorque(robot,q_init);
Cor = velocityProduct(robot,q_init,q_init);
M = massMatrix(robot,q_init);

J_RF = geometricJacobian(robot,q_init,'R_Foot_Link');
J_LF = geometricJacobian(robot,q_init,'L_Foot_Link');
