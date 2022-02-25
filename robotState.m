clear all
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

%%%%%%%%foot step input%%%%%%%%%%
X_direction = 1.0;
step_length = 0.3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

foot_step_number = fix(X_direction / step_length);
final_step_length = rem(X_direction, step_length);
if (rem(X_direction,step_length) ~= 0)
    foot_step_number = foot_step_number + 1;    
end

foot_step(1,1) = T_LF(1,4);
foot_step(1,2) = T_LF(2,4);

for i = 2:foot_step_number
    foot_step(i,1) = step_length + foot_step(i-1,1);
    foot_step(i,2) = foot_step(i-1,2) * (-1);
end
if (final_step_length == 0)
    foot_step(foot_step_number+1,1) = step_length + foot_step(foot_step_number,1);
    foot_step(foot_step_number+1,2) = foot_step(foot_step_number,2) * (-1);
else
    foot_step(foot_step_number+1,1) = final_step_length + foot_step(foot_step_number,1);
    foot_step(foot_step_number+1,2) = foot_step(foot_step_number,2) * (-1);
end