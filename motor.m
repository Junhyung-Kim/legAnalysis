clear all
close all
clc

for i = 1:856
   hipyaw(i,1) = 0.001 * i; 
   
end

for i = 1:1220
   hiproll(i,1) = 0.001 * i; 
end

for i = 1:1540
   hippitch(i,1) = 0.001 * i; 
   knee(i,1) = 0.001 * i; 
   anklepitch(i,1) = 0.001 * i; 
end

for i = 1:429
   ankleroll(i,1) = 0.001 * i; 
end

hipyaw_res = 0.753;
hiproll_res = 0.949;
hippitch_res = 0.458;
knee_res = 0.458;
anklepitch_res = 0.458;
ankleroll_res = 0.615;

hipyaw_ki = 0.121;
hiproll_ki = 0.184;
hippitch_ki = 0.157;
knee_ki = 0.157;
anklepitch_ki = 0.157;
ankleroll_ki = 0.061;

for i = 1:856
   hipyaw(i,2) = hipyaw(i,1)/hipyaw_ki*hipyaw_res;
end

for i = 1:1220
   hiproll(i,2) = hiproll(i,1)/hiproll_ki*hiproll_res; 
end

for i = 1:1540
   hippitch(i,2) = hippitch(i,1)/hippitch_ki*hippitch_res;
   knee(i,2) = knee(i,1)/knee_ki*knee_res; 
   anklepitch(i,2) = anklepitch(i,1)/anklepitch_ki*anklepitch_res;
end

for i = 1:429
   ankleroll(i,2) = ankleroll(i,1)/ankleroll_ki*ankleroll_res; 
end

Vcc = 48;

for i = 1:856
   hipyaw(i,3) = Vcc - hipyaw(i,2);
end

for i = 1:1220
   hiproll(i,3) = Vcc - hiproll(i,2);
end

for i = 1:1540
   hippitch(i,3) = Vcc - hippitch(i,2);
   knee(i,3) = Vcc - knee(i,2); 
   anklepitch(i,3) = Vcc - anklepitch(i,2);
end

for i = 1:429
   ankleroll(i,3) = Vcc - ankleroll(i,2);
end

for i = 1:856
   hipyaw(i,4) = hipyaw(i,3)/hipyaw_ki;
end

for i = 1:1220
   hiproll(i,4) = hiproll(i,3)/hiproll_ki;
end

for i = 1:1540
   hippitch(i,4) = hippitch(i,3)/hippitch_ki;
   knee(i,4) = knee(i,3)/knee_ki; 
   anklepitch(i,4) = anklepitch(i,3)/anklepitch_ki;
end

for i = 1:429
   ankleroll(i,4) = ankleroll(i,3)/ankleroll_ki;
end

for i = 1:856
   hipyaw(i,5) = hipyaw(i,4)*60/(2*pi);
end

for i = 1:1220
   hiproll(i,5) = hiproll(i,4)*60/(2*pi);
end

for i = 1:1540
   hippitch(i,5) = hippitch(i,4)*60/(2*pi);
   knee(i,5) = knee(i,4)*60/(2*pi);
   anklepitch(i,5) = anklepitch(i,4)*60/(2*pi);
end
for i = 1:429
   ankleroll(i,5) = ankleroll(i,4)*60/(2*pi);
end

hipyaw_angv(1,1) = 40;
hipyaw_torque(1,1) = 0.0;
hiproll_angv(1,1) = 26;
hiproll_torque(1,1) = 0.0;
hippitch_angv(1,1) = 30;
hippitch_torque(1,1) = 0.0;
knee_angv(1,1) = 30;
knee_torque(1,1) = 0.0;
anklepitch_angv(1,1) = 30;
anklepitch_torque(1,1) = 0.0;
ankleroll_angv(1,1) = 80;
ankleroll_torque(1,1) = 0.0;

for i = 1:856
   hipyaw_angv(i+1,1) = hipyaw(i,5)/100;
   hipyaw_torque(i+1,1) = hipyaw(i,1)*100*0.8;
end

for i = 1:1220
   hiproll_angv(i+1,1) =hiproll(i,5)/100;
   hiproll_torque(i+1,1) =hiproll(i,1)*100*0.8;
end

for i = 1:1540
   hippitch_angv(i+1,1) =hippitch(i,5)/100;
   knee_angv(i+1,1) =knee(i,5)/100;
   anklepitch_angv(i+1,1) =anklepitch(i,5)/100;
   hippitch_torque(i+1,1) =hippitch(i,1)*100*0.8;
   knee_torque(i+1,1) =knee(i,1)*100*0.8;
   anklepitch_torque(i+1,1) =anklepitch(i,1)*100*0.8;
end
for i = 1:429
   ankleroll_angv(i+1,1) =ankleroll(i,5)/100;
   ankleroll_torque(i+1,1) =ankleroll(i,1)*100*0.8;
end

   hipyaw_angv(858,1) = 0;
   hipyaw_torque(858,1) =   hipyaw_torque(857,1);
   hiproll_angv(1222,1) = 0;
   hiproll_torque(1222,1) = hiproll_torque(1221,1);
   hippitch_angv(1542,1) = 0;
   hippitch_torque(1542,1) = hippitch_torque(1541,1);
   knee_angv(1542,1) = 0;
   knee_torque(1542,1) = knee_torque(1541,1);
   anklepitch_angv(1542,1) = 0;
   anklepitch_torque(1542,1) = anklepitch_torque(1541,1);
   ankleroll_angv(431,1) = 0;
   ankleroll_torque(431,1) = ankleroll_torque(430,1);
   walking1020 = importdata("0_tocabi_-8.txt");
for i = 1:9999
   walking1020(i,2+4) =  walking1020(i,2+4)*60/(2*pi);
   walking1020(i,4+4) =  walking1020(i,4+4)*60/(2*pi);
   walking1020(i,6+4) =  walking1020(i,6+4)*60/(2*pi);
   walking1020(i,8+4) =  walking1020(i,8+4)*60/(2*pi);
   walking1020(i,10+4) =  walking1020(i,10+4)*60/(2*pi); 
   walking1020(i,12+4) =  walking1020(i,12+4)*60/(2*pi);
end

figure(2);
sgtitle('Pose2, Step time : 1s, Step length: 0.2m');
subplot(2,3,1);
plot(hipyaw_angv(:,1),hipyaw_torque(:,1));
hold on
title('Hip Yaw Joint')
plot(hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),hipyaw_torque(:,1));
plot(walking1020(:,2+4),walking1020(:,1+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 


subplot(2,3,2);
plot(hiproll_angv(:,1),hiproll_torque(:,1));
hold on
title('Hip Roll Joint')
plot(hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),hiproll_torque(:,1));
plot(walking1020(:,4+4),walking1020(:,3+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 


subplot(2,3,3);
plot(hippitch_angv(:,1),hippitch_torque(:,1));
hold on
title('Hip Pitch Joint')
plot(hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),hippitch_torque(:,1));
plot(walking1020(:,6+4),walking1020(:,5+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 


subplot(2,3,4);
hold on
title('Knee Pitch Joint')
plot(knee_angv(:,1),knee_torque(:,1));
plot(knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),knee_torque(:,1));
plot(walking1020(:,8+4),walking1020(:,7+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 


subplot(2,3,5);
hold on
title('Ankle Pitch Joint')
plot(anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(walking1020(:,10+4),walking1020(:,9+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 


subplot(2,3,6);
hold on
title('Ankle Roll Joint')
plot(ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(walking1020(:,12+4),walking1020(:,11+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 
saveas(gcf,'Pose2, Step time : 1s, Step length: 0.2m.png')

walking1030 = importdata("0_tocabi_-12.txt");
for i = 1:9999
   walking1030(i,2+4) =  walking1030(i,2+4)*60/(2*pi);
   walking1030(i,4+4) =  walking1030(i,4+4)*60/(2*pi);
   walking1030(i,6+4) =  walking1030(i,6+4)*60/(2*pi);
   walking1030(i,8+4) =  walking1030(i,8+4)*60/(2*pi);
   walking1030(i,10+4) =  walking1030(i,10+4)*60/(2*pi); 
   walking1030(i,12+4) =  walking1030(i,12+4)*60/(2*pi);
end

figure(6);
subplot(2,3,1);
sgtitle('Pose2, Step time : 1s, Step length: 0.3m');
plot(hipyaw_angv(:,1),hipyaw_torque(:,1));
hold on
title('Hip Yaw Joint')
plot(hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),hipyaw_torque(:,1));
plot(walking1030(:,2+4),walking1030(:,1+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 


subplot(2,3,2);
plot(hiproll_angv(:,1),hiproll_torque(:,1));
hold on
title('Hip Roll Joint')
plot(hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),hiproll_torque(:,1));
plot(walking1030(:,4+4),walking1030(:,3+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 


subplot(2,3,3);
plot(hippitch_angv(:,1),hippitch_torque(:,1));
hold on
title('Hip Pitch Joint')
plot(hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),hippitch_torque(:,1));
plot(walking1030(:,6+4),walking1030(:,5+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 


subplot(2,3,4);
hold on
title('Knee Pitch Joint')
plot(knee_angv(:,1),knee_torque(:,1));
plot(knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),knee_torque(:,1));
plot(walking1030(:,8+4),walking1030(:,7+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 


subplot(2,3,5);
hold on
title('Ankle Pitch Joint')
plot(anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(walking1030(:,10+4),walking1030(:,9+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 


subplot(2,3,6);
hold on
title('Ankle Roll Joint')
plot(ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(walking1030(:,12+4),walking1030(:,11+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 
saveas(gcf,'Pose2, Step time : 1s, Step length: 0.3m.png')

walking1040 = importdata("0_tocabi_-9.txt");
for i = 1:9999
   walking1040(i,2+4) =  walking1040(i,2+4)*60/(2*pi);
   walking1040(i,4+4) =  walking1040(i,4+4)*60/(2*pi);
   walking1040(i,6+4) =  walking1040(i,6+4)*60/(2*pi);
   walking1040(i,8+4) =  walking1040(i,8+4)*60/(2*pi);
   walking1040(i,10+4) =  walking1040(i,10+4)*60/(2*pi); 
   walking1040(i,12+4) =  walking1040(i,12+4)*60/(2*pi);
end

figure(3);
subplot(2,3,1);
sgtitle('Pose2, Step time : 1s, Step length: 0.4m');
plot(hipyaw_angv(:,1),hipyaw_torque(:,1));
hold on
title('Hip Yaw Joint')
plot(hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),hipyaw_torque(:,1));
plot(walking1040(:,2+4),walking1040(:,1+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 


subplot(2,3,2);
plot(hiproll_angv(:,1),hiproll_torque(:,1));
hold on
title('Hip Roll Joint')
plot(hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),hiproll_torque(:,1));
plot(walking1040(:,4+4),walking1040(:,3+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 


subplot(2,3,3);
plot(hippitch_angv(:,1),hippitch_torque(:,1));
hold on
title('Hip Pitch Joint')
plot(hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),hippitch_torque(:,1));
plot(walking1040(:,6+4),walking1040(:,5+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 


subplot(2,3,4);
hold on
title('Knee Pitch Joint')
plot(knee_angv(:,1),knee_torque(:,1));
plot(knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),knee_torque(:,1));
plot(walking1040(:,8+4),walking1040(:,7+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,5);
hold on
title('Ankle Pitch Joint')
plot(anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(walking1040(:,10+4),walking1040(:,9+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 


subplot(2,3,6);
hold on
title('Ankle Roll Joint')
plot(ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(walking1040(:,12+4),walking1040(:,11+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

saveas(gcf,'Pose2, Step time : 1s, Step length: 0.4m.png')
walking1050 = importdata("0_tocabi_-10.txt");
for i = 1:9999
   walking1050(i,2+4) =  walking1050(i,2+4)*60/(2*pi);
   walking1050(i,4+4) =  walking1050(i,4+4)*60/(2*pi);
   walking1050(i,6+4) =  walking1050(i,6+4)*60/(2*pi);
   walking1050(i,8+4) =  walking1050(i,8+4)*60/(2*pi);
   walking1050(i,10+4) =  walking1050(i,10+4)*60/(2*pi); 
   walking1050(i,12+4) =  walking1050(i,12+4)*60/(2*pi);
end

figure(4);
sgtitle('Pose2, Step time : 1s, Step length: 0.5m');
subplot(2,3,1);
plot(hipyaw_angv(:,1),hipyaw_torque(:,1));
hold on
title('Hip Yaw Joint')
plot(hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),hipyaw_torque(:,1));
plot(walking1050(:,2+4),walking1050(:,1+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,2);
plot(hiproll_angv(:,1),hiproll_torque(:,1));
hold on
title('Hip Roll Joint')
plot(hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),hiproll_torque(:,1));
plot(walking1050(:,4+4),walking1050(:,3+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,3);
plot(hippitch_angv(:,1),hippitch_torque(:,1));
hold on
title('Hip Pitch Joint')
plot(hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),hippitch_torque(:,1));
plot(walking1050(:,6+4),walking1050(:,5+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,4);
hold on
title('Knee Pitch Joint')
plot(knee_angv(:,1),knee_torque(:,1));
plot(knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),knee_torque(:,1));
plot(walking1050(:,8+4),walking1050(:,7+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,5);
hold on
title('Ankle Pitch Joint')
plot(anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(walking1050(:,10+4),walking1050(:,9+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,6);
hold on
title('Ankle Roll Joint')
plot(ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(walking1050(:,12+4),walking1050(:,11+4));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 
saveas(gcf,'Pose2, Step time : 1s, Step length: 0.5m.png')

walking0810 = importdata("0_tocabi_-13.txt");
for i = 1:8399
   walking0810(i,2) =  walking0810(i,2)*60/(2*pi);
   walking0810(i,4) =  walking0810(i,4)*60/(2*pi);
   walking0810(i,6) =  walking0810(i,6)*60/(2*pi);
   walking0810(i,8) =  walking0810(i,8)*60/(2*pi);
   walking0810(i,10) =  walking0810(i,10)*60/(2*pi); 
   walking0810(i,12) =  walking0810(i,12)*60/(2*pi);
end

figure(7);
sgtitle('Pose2, Step time : 0.8s, Step length: 0.1m');
subplot(2,3,1);
plot(hipyaw_angv(:,1),hipyaw_torque(:,1));
hold on
title('Hip Yaw Joint')
plot(hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),hipyaw_torque(:,1));
plot(walking0810(:,2),walking0810(:,1));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,2);
plot(hiproll_angv(:,1),hiproll_torque(:,1));
hold on
title('Hip Roll Joint')
plot(hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),hiproll_torque(:,1));
plot(walking0810(:,4),walking0810(:,3));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,3);
plot(hippitch_angv(:,1),hippitch_torque(:,1));
hold on
title('Hip Pitch Joint')
plot(hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),hippitch_torque(:,1));
plot(walking0810(:,6),walking0810(:,5));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,4);
hold on
title('Knee Pitch Joint')
plot(knee_angv(:,1),knee_torque(:,1));
plot(knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),knee_torque(:,1));
plot(walking0810(:,8),walking0810(:,7));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,5);
hold on
title('Ankle Pitch Joint')
plot(anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(walking0810(:,10),walking0810(:,9));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,6);
hold on
title('Ankle Roll Joint')
plot(ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(walking0810(:,12),walking0810(:,11));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 
saveas(gcf,'Pose2, Step time : 0.8s, Step length: 0.1m.png')
walking0820 = importdata("0_tocabi_-14.txt");
for i = 1:8399
   walking0820(i,2) =  walking0820(i,2)*60/(2*pi);
   walking0820(i,4) =  walking0820(i,4)*60/(2*pi);
   walking0820(i,6) =  walking0820(i,6)*60/(2*pi);
   walking0820(i,8) =  walking0820(i,8)*60/(2*pi);
   walking0820(i,10) =  walking0820(i,10)*60/(2*pi); 
   walking0820(i,12) =  walking0820(i,12)*60/(2*pi);
end

figure(9);
sgtitle('Pose2, Step time : 0.8s, Step length: 0.2m');
subplot(2,3,1);
plot(hipyaw_angv(:,1),hipyaw_torque(:,1));
hold on
title('Hip Yaw Joint')
plot(hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),hipyaw_torque(:,1));
plot(walking0820(:,2),walking0820(:,1));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,2);
plot(hiproll_angv(:,1),hiproll_torque(:,1));
hold on
title('Hip Roll Joint')
plot(hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),hiproll_torque(:,1));
plot(walking0820(:,4),walking0820(:,3));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,3);
plot(hippitch_angv(:,1),hippitch_torque(:,1));
hold on
title('Hip Pitch Joint')
plot(hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),hippitch_torque(:,1));
plot(walking0820(:,6),walking0820(:,5));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,4);
hold on
title('Knee Pitch Joint')
plot(knee_angv(:,1),knee_torque(:,1));
plot(knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),knee_torque(:,1));
plot(walking0820(:,8),walking0820(:,7));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,5);
hold on
title('Ankle Pitch Joint')
plot(anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(walking0820(:,10),walking0820(:,9));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,6);
hold on
title('Ankle Roll Joint')
plot(ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(walking0820(:,12),walking0820(:,11));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 
saveas(gcf,'Pose2, Step time : 0.8s, Step length: 0.2m.png')
walking0830 = importdata("0_tocabi_-18.txt");
for i = 1:8399
   walking0830(i,2) =  walking0830(i,2)*60/(2*pi);
   walking0830(i,4) =  walking0830(i,4)*60/(2*pi);
   walking0830(i,6) =  walking0830(i,6)*60/(2*pi);
   walking0830(i,8) =  walking0830(i,8)*60/(2*pi);
   walking0830(i,10) =  walking0830(i,10)*60/(2*pi); 
   walking0830(i,12) =  walking0830(i,12)*60/(2*pi);
end

figure(10);
sgtitle('Pose2, Step time : 0.8s, Step length: 0.3m');
subplot(2,3,1);
plot(hipyaw_angv(:,1),hipyaw_torque(:,1));
hold on
title('Hip Yaw Joint')
plot(hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),hipyaw_torque(:,1));
plot(walking0830(:,2),walking0830(:,1));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,2);
plot(hiproll_angv(:,1),hiproll_torque(:,1));
hold on
title('Hip Roll Joint')
plot(hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),hiproll_torque(:,1));
plot(walking0830(:,4),walking0830(:,3));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,3);
plot(hippitch_angv(:,1),hippitch_torque(:,1));
hold on
title('Hip Pitch Joint')
plot(hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),hippitch_torque(:,1));
plot(walking0830(:,6),walking0830(:,5));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,4);
hold on
title('Knee Pitch Joint')
plot(knee_angv(:,1),knee_torque(:,1));
plot(knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),knee_torque(:,1));
plot(walking0830(:,8),walking0830(:,7));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,5);
hold on
title('Ankle Pitch Joint')
plot(anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(walking0830(:,10),walking0830(:,9));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,6);
hold on
title('Ankle Roll Joint')
plot(ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(walking0830(:,12),walking0830(:,11));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 
saveas(gcf,'Pose2, Step time : 0.8s, Step length: 0.3m.png')
walking0840 = importdata("0_tocabi_-16.txt");
for i = 1:8399
   walking0840(i,2) =  walking0840(i,2)*60/(2*pi);
   walking0840(i,4) =  walking0840(i,4)*60/(2*pi);
   walking0840(i,6) =  walking0840(i,6)*60/(2*pi);
   walking0840(i,8) =  walking0840(i,8)*60/(2*pi);
   walking0840(i,10) =  walking0840(i,10)*60/(2*pi); 
   walking0840(i,12) =  walking0840(i,12)*60/(2*pi);
end

figure(11);
sgtitle('Pose2, Step time : 0.8s, Step length: 0.4m');
subplot(2,3,1);
plot(hipyaw_angv(:,1),hipyaw_torque(:,1));
hold on
title('Hip Yaw Joint')
plot(hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),hipyaw_torque(:,1));
plot(walking0840(:,2),walking0840(:,1));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,2);
plot(hiproll_angv(:,1),hiproll_torque(:,1));
hold on
title('Hip Roll Joint')
plot(hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),hiproll_torque(:,1));
plot(walking0840(:,4),walking0840(:,3));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,3);
plot(hippitch_angv(:,1),hippitch_torque(:,1));
hold on
title('Hip Pitch Joint')
plot(hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),hippitch_torque(:,1));
plot(walking0840(:,6),walking0840(:,5));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,4);
hold on
title('Knee Pitch Joint')
plot(knee_angv(:,1),knee_torque(:,1));
plot(knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),knee_torque(:,1));
plot(walking0840(:,8),walking0840(:,7));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,5);
hold on
title('Ankle Pitch Joint')
plot(anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(walking0840(:,10),walking0840(:,9));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,6);
hold on
title('Ankle Roll Joint')
plot(ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(walking0840(:,12),walking0840(:,11));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 
saveas(gcf,'Pose2, Step time : 0.8s, Step length: 0.4m.png')
walking0850 = importdata("0_tocabi_-17.txt");
for i = 1:8399
   walking0850(i,2) =  walking0850(i,2)*60/(2*pi);
   walking0850(i,4) =  walking0850(i,4)*60/(2*pi);
   walking0850(i,6) =  walking0850(i,6)*60/(2*pi);
   walking0850(i,8) =  walking0850(i,8)*60/(2*pi);
   walking0850(i,10) =  walking0850(i,10)*60/(2*pi); 
   walking0850(i,12) =  walking0850(i,12)*60/(2*pi);
end

figure(12);
sgtitle('Pose2, Step time : 0.8s, Step length: 0.5m');
subplot(2,3,1);
plot(hipyaw_angv(:,1),hipyaw_torque(:,1));
hold on
title('Hip Yaw Joint')
plot(hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),hipyaw_torque(:,1));
plot(walking0850(:,2),walking0850(:,1));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,2);
plot(hiproll_angv(:,1),hiproll_torque(:,1));
hold on
title('Hip Roll Joint')
plot(hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),hiproll_torque(:,1));
plot(walking0850(:,4),walking0850(:,3));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,3);
plot(hippitch_angv(:,1),hippitch_torque(:,1));
hold on
title('Hip Pitch Joint')
plot(hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),hippitch_torque(:,1));
plot(walking0850(:,6),walking0850(:,5));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,4);
hold on
title('Knee Pitch Joint')
plot(knee_angv(:,1),knee_torque(:,1));
plot(knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),knee_torque(:,1));
plot(walking0850(:,8),walking0850(:,7));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,5);
hold on
title('Ankle Pitch Joint')
plot(anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(walking0850(:,10),walking0850(:,9));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,6);
hold on
title('Ankle Roll Joint')
plot(ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(walking0850(:,12),walking0850(:,11));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 
saveas(gcf,'Pose2, Step time : 0.8s, Step length: 0.6m.png')

walking0610 = importdata("0_tocabi_-19.txt");
for i = 1:6799
   walking0610(i,2) =  walking0610(i,2)*60/(2*pi);
   walking0610(i,4) =  walking0610(i,4)*60/(2*pi);
   walking0610(i,6) =  walking0610(i,6)*60/(2*pi);
   walking0610(i,8) =  walking0610(i,8)*60/(2*pi);
   walking0610(i,10) =  walking0610(i,10)*60/(2*pi); 
   walking0610(i,12) =  walking0610(i,12)*60/(2*pi);
end

figure(13);
sgtitle('Pose2, Step time : 0.6s, Step length: 0.1m');
subplot(2,3,1);
plot(hipyaw_angv(:,1),hipyaw_torque(:,1));
hold on
title('Hip Yaw Joint')
plot(hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),hipyaw_torque(:,1));
plot(walking0610(:,2),walking0610(:,1));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,2);
plot(hiproll_angv(:,1),hiproll_torque(:,1));
hold on
title('Hip Roll Joint')
plot(hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),hiproll_torque(:,1));
plot(walking0610(:,4),walking0610(:,3));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,3);
plot(hippitch_angv(:,1),hippitch_torque(:,1));
hold on
title('Hip Pitch Joint')
plot(hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),hippitch_torque(:,1));
plot(walking0610(:,6),walking0610(:,5));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,4);
hold on
title('Knee Pitch Joint')
plot(knee_angv(:,1),knee_torque(:,1));
plot(knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),knee_torque(:,1));
plot(walking0610(:,8),walking0610(:,7));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,5);
hold on
title('Ankle Pitch Joint')
plot(anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(walking0610(:,10),walking0610(:,9));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,6);
hold on
title('Ankle Roll Joint')
plot(ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(walking0610(:,12),walking0610(:,11));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 
saveas(gcf,'Pose2, Step time : 0.6s, Step length: 0.1m.png')

walking0620 = importdata("0_tocabi_-20.txt");
for i = 1:6799
   walking0620(i,2) =  walking0620(i,2)*60/(2*pi);
   walking0620(i,4) =  walking0620(i,4)*60/(2*pi);
   walking0620(i,6) =  walking0620(i,6)*60/(2*pi);
   walking0620(i,8) =  walking0620(i,8)*60/(2*pi);
   walking0620(i,10) =  walking0620(i,10)*60/(2*pi); 
   walking0620(i,12) =  walking0620(i,12)*60/(2*pi);
end

figure(14);
sgtitle('Pose2, Step time : 0.6s, Step length: 0.2m');
subplot(2,3,1);
plot(hipyaw_angv(:,1),hipyaw_torque(:,1));
hold on
title('Hip Yaw Joint')
plot(hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),hipyaw_torque(:,1));
plot(walking0620(:,2),walking0620(:,1));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,2);
plot(hiproll_angv(:,1),hiproll_torque(:,1));
hold on
title('Hip Roll Joint')
plot(hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),hiproll_torque(:,1));
plot(walking0620(:,4),walking0620(:,3));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,3);
plot(hippitch_angv(:,1),hippitch_torque(:,1));
hold on
title('Hip Pitch Joint')
plot(hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),hippitch_torque(:,1));
plot(walking0620(:,6),walking0620(:,5));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,4);
hold on
title('Knee Pitch Joint')
plot(knee_angv(:,1),knee_torque(:,1));
plot(knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),knee_torque(:,1));
plot(walking0620(:,8),walking0620(:,7));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,5);
hold on
plot(anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(walking0620(:,10),walking0620(:,9));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,6);
hold on
title('Ankle Roll Joint')
plot(ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(walking0620(:,12),walking0620(:,11));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 
saveas(gcf,'Pose2, Step time : 0.6s, Step length: 0.2m.png')
walking0630 = importdata("0_tocabi_-21.txt");
for i = 1:6799
   walking0630(i,2) =  walking0630(i,2)*60/(2*pi);
   walking0630(i,4) =  walking0630(i,4)*60/(2*pi);
   walking0630(i,6) =  walking0630(i,6)*60/(2*pi);
   walking0630(i,8) =  walking0630(i,8)*60/(2*pi);
   walking0630(i,10) =  walking0630(i,10)*60/(2*pi); 
   walking0630(i,12) =  walking0630(i,12)*60/(2*pi);
end

figure(15);
sgtitle('Pose2, Step time : 0.6s, Step length: 0.3m');
subplot(2,3,1);
plot(hipyaw_angv(:,1),hipyaw_torque(:,1));
hold on
title('Hip Yaw Joint')
plot(hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),hipyaw_torque(:,1));
plot(walking0630(:,2),walking0630(:,1));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,2);
plot(hiproll_angv(:,1),hiproll_torque(:,1));
hold on
title('Hip Roll Joint')
plot(hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),hiproll_torque(:,1));
plot(walking0630(:,4),walking0630(:,3));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,3);
plot(hippitch_angv(:,1),hippitch_torque(:,1));
hold on
title('Hip Pitch Joint')
plot(hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),hippitch_torque(:,1));
plot(walking0630(:,6),walking0630(:,5));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,4);
hold on
title('Knee Pitch Joint')
plot(knee_angv(:,1),knee_torque(:,1));
plot(knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),knee_torque(:,1));
plot(walking0630(:,8),walking0630(:,7));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,5);
hold on
title('Ankle Pitch Joint')
plot(anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(walking0630(:,10),walking0630(:,9));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,6);
hold on
title('Ankle Roll Joint')
plot(ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(walking0630(:,12),walking0630(:,11));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 
saveas(gcf,'Pose2, Step time : 0.6s, Step length: 0.3m.png')
walking0640 = importdata("0_tocabi_-22.txt");
for i = 1:6799
   walking0640(i,2) =  walking0640(i,2)*60/(2*pi);
   walking0640(i,4) =  walking0640(i,4)*60/(2*pi);
   walking0640(i,6) =  walking0640(i,6)*60/(2*pi);
   walking0640(i,8) =  walking0640(i,8)*60/(2*pi);
   walking0640(i,10) =  walking0640(i,10)*60/(2*pi); 
   walking0640(i,12) =  walking0640(i,12)*60/(2*pi);
end

figure(16);
sgtitle('Pose2, Step time : 0.6s, Step length: 0.4m');
subplot(2,3,1);
plot(hipyaw_angv(:,1),hipyaw_torque(:,1));
hold on
plot(hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),hipyaw_torque(:,1));
plot(walking0640(:,2),walking0640(:,1));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,2);
plot(hiproll_angv(:,1),hiproll_torque(:,1));
hold on
title('Hip Roll Joint')
plot(hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),hiproll_torque(:,1));
plot(walking0640(:,4),walking0640(:,3));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,3);
plot(hippitch_angv(:,1),hippitch_torque(:,1));
hold on
title('Hip Pitch Joint')
plot(hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),hippitch_torque(:,1));
plot(walking0640(:,6),walking0640(:,5));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,4);
hold on
title('Knee Pitch Joint')
plot(knee_angv(:,1),knee_torque(:,1));
plot(knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),knee_torque(:,1));
plot(walking0640(:,8),walking0640(:,7));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,5);
hold on
title('Ankle Pitch Joint')
plot(anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(walking0640(:,10),walking0640(:,9));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,6);
hold on
title('Ankle Roll Joint')
plot(ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(walking0640(:,12),walking0640(:,11));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 
saveas(gcf,'Pose2, Step time : 0.6s, Step length: 0.4m.png')
walking0650 = importdata("0_tocabi_-23.txt");
for i = 1:6799
   walking0650(i,2) =  walking0650(i,2)*60/(2*pi);
   walking0650(i,4) =  walking0650(i,4)*60/(2*pi);
   walking0650(i,6) =  walking0650(i,6)*60/(2*pi);
   walking0650(i,8) =  walking0650(i,8)*60/(2*pi);
   walking0650(i,10) =  walking0650(i,10)*60/(2*pi); 
   walking0650(i,12) =  walking0650(i,12)*60/(2*pi);
end

figure(17);
sgtitle('Pose2, Step time : 0.6s, Step length: 0.5m');
subplot(2,3,1);
plot(hipyaw_angv(:,1),hipyaw_torque(:,1));
hold on
title('Hip Yaw Joint')
plot(hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),-1*hipyaw_torque(:,1));
plot(-1*hipyaw_angv(:,1),hipyaw_torque(:,1));
plot(walking0650(:,2),walking0650(:,1));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,2);
plot(hiproll_angv(:,1),hiproll_torque(:,1));
hold on
title('Hip Roll Joint')
plot(hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),-1*hiproll_torque(:,1));
plot(-1*hiproll_angv(:,1),hiproll_torque(:,1));
plot(walking0650(:,4),walking0650(:,3));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,3);
plot(hippitch_angv(:,1),hippitch_torque(:,1));
hold on
title('Hip Pitch Joint')
plot(hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),-1*hippitch_torque(:,1));
plot(-1*hippitch_angv(:,1),hippitch_torque(:,1));
plot(walking0650(:,6),walking0650(:,5));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,4);
hold on
title('Knee Pitch Joint')
plot(knee_angv(:,1),knee_torque(:,1));
plot(knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),-1*knee_torque(:,1));
plot(-1*knee_angv(:,1),knee_torque(:,1));
plot(walking0650(:,8),walking0650(:,7));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,5);
hold on
title('Ankle Pitch Joint')
plot(anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),-1*anklepitch_torque(:,1));
plot(-1*anklepitch_angv(:,1),anklepitch_torque(:,1));
plot(walking0650(:,10),walking0650(:,9));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 

subplot(2,3,6);
hold on
title('Ankle Roll Joint')
plot(ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),-1*ankleroll_torque(:,1));
plot(-1*ankleroll_angv(:,1),ankleroll_torque(:,1));
plot(walking0650(:,12),walking0650(:,11));
xlabel('Joint Velocity (rpm)') 
ylabel('Joint Torque (Nm)') 
saveas(gcf,'Pose2, Step time : 0.6s, Step length: 0.5m.png')