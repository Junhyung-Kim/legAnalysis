clear all;
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
step_length = 0.2;
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

%% ZMP generator

%%%%%%%%%%%%%%input%%%%%%%%%%%%%%%%%%%
hz_ = 1000;
% walking_tick_ = 1/hz;
total_tick = 10000;

t_total_t = 1.2;
t_temp_t = 1.0;
t_double = 0.1;

t_total_ = t_total_t * hz_;
t_temp_ = t_temp_t * hz_;
t_start_ = t_temp_ + 1;
t_last_ = t_temp_ + t_total_;
t_double_1 = 0.1 * hz_;
t_double_2 = 0.1 * hz_;

ZMP_cutoff = 0.005;
ref_zmp_ = zeros(total_tick,2);
walking_tick = zeros(1,total_tick);

current_step_num_ = 0;
% planning_step_num = 3;

zc_ = 0.727822;
wn_ = sqrt(9.81/zc_);

foot_height_ = 0.5;

PHASE_variable = zeros(total_tick);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for walking_tick_ = 0:total_tick-1
    walking_tick(walking_tick_+1) = walking_tick_;
    if walking_tick_ <= t_temp_
        
        if walking_tick_ < 0.5 * hz_
            ref_zmp_(walking_tick_+1,1) = T_LF(1,4);
            ref_zmp_(walking_tick_+1,2) = 0;
        elseif walking_tick_ < 1.5 * hz_
            del_x = walking_tick_ - 0.5 * hz_;
            ref_zmp_(walking_tick_+1,1) = T_LF(1,4);
            ref_zmp_(walking_tick_+1,2) = 0;
        else
            ref_zmp_(walking_tick_+1,1) = T_LF(1,4);
            ref_zmp_(walking_tick_+1,2) = 0;
        end
    elseif walking_tick_ <= t_last_

        if current_step_num_ == 0
            A_ = foot_step(current_step_num_+1,2);
            B_ = (foot_step(current_step_num_+1,1) + foot_step(current_step_num_+2,1)) / 2;
            Kx_ = (B_ * t_double * wn_) / (t_double * wn_ + tanh(wn_ * (t_total_t/2 - t_double)));
            Ky_ = A_ * t_double * wn_ * tanh(wn_ * (t_total_t/2 - t_double)) / (1 + t_double * wn_ * tanh(wn_ * (t_total_t/2 - t_double)));
            if walking_tick_ < t_start_ + t_double_1
                ref_zmp_(walking_tick_+1,1) = foot_step(current_step_num_+1,1);
                ref_zmp_(walking_tick_+1,2) = Ky_ / t_double_1 * (walking_tick_- t_start_);
            elseif (walking_tick_ < t_last_ - t_double_2) && (walking_tick_ >= t_start_ + t_double_1)
                ref_zmp_(walking_tick_+1,1) = foot_step(current_step_num_+1,1);
                ref_zmp_(walking_tick_+1,2) = A_;
            elseif (walking_tick_ <= t_last_) && (walking_tick_ >= t_last_ - t_double_2)
                ref_zmp_(walking_tick_+1,1) = (B_ - Kx_) + (Kx_ / t_double_2) * (walking_tick_ - t_start_ - (t_total_ - t_double_2)); 
                ref_zmp_(walking_tick_+1,2) = (Ky_ / t_double_2) * (t_total_ - (walking_tick_ - t_start_));
            end
            if walking_tick_ == t_last_
                current_step_num_ = current_step_num_ + 1;
                t_start_ = t_start_ + t_total_;
                t_last_ = t_last_ + t_total_;
            end

        elseif current_step_num_ < foot_step_number 
            A_ = foot_step(current_step_num_+1,2);
            B_ = (foot_step(current_step_num_,1) + foot_step(current_step_num_+1,1)) / 2 - foot_step(current_step_num_,1);
            Kx_ = (B_ * t_double * wn_) / (t_double * wn_ + tanh(wn_ * (t_total_t/2 - t_double)));
            Ky_ = A_ * t_double * wn_ * tanh(wn_ * (t_total_t/2 - t_double)) / (1 + t_double * wn_ * tanh(wn_ * (t_total_t/2 - t_double)));
            if walking_tick_ < t_start_ + t_double_1
                ref_zmp_(walking_tick_+1,1) = (foot_step(current_step_num_,1) + foot_step(current_step_num_+1,1)) / 2 + (Kx_ / t_double_1) * (walking_tick_ - t_start_);
                ref_zmp_(walking_tick_+1,2) = Ky_ / t_double_1 * (walking_tick_- t_start_);
            elseif walking_tick_ < t_last_ - t_double_2
                ref_zmp_(walking_tick_+1,1) = foot_step(current_step_num_+1,1);
                ref_zmp_(walking_tick_+1,2) = A_;
            else
                ref_zmp_(walking_tick_+1,1) = ((foot_step(current_step_num_+1,1) + foot_step(current_step_num_+2,1)) / 2 - Kx_) + (Kx_ / t_double_2) * (walking_tick_ - t_start_ - (t_total_ - t_double_2)); 
                ref_zmp_(walking_tick_+1,2) = (Ky_ / t_double_2) * (t_total_ - (walking_tick_ - t_start_));
            end
            if walking_tick_ == t_last_
                current_step_num_ = current_step_num_ + 1;
                t_start_ = t_start_ + t_total_;
                t_last_ = t_last_ + t_total_;
            end
        elseif current_step_num_ == foot_step_number
            if final_step_length == 0
                A_ = foot_step(current_step_num_+1,2);
                B_ = (foot_step(current_step_num_,1) + foot_step(current_step_num_+1,1)) / 2 - foot_step(current_step_num_,1);
                Kx_ = (B_ * t_double * wn_) / (t_double * wn_ + tanh(wn_ * (t_total_t/2 - t_double)));
                Ky_ = A_ * t_double * wn_ * tanh(wn_ * (t_total_t/2 - t_double)) / (1 + t_double * wn_ * tanh(wn_ * (t_total_t/2 - t_double)));
                if walking_tick_ < t_start_ + t_double_1
                    ref_zmp_(walking_tick_+1,1) = (foot_step(current_step_num_,1) + foot_step(current_step_num_+1,1)) / 2 + (Kx_ / t_double_1) * (walking_tick_ - t_start_);
                    ref_zmp_(walking_tick_+1,2) = Ky_ / t_double_1 * (walking_tick_- t_start_);
                elseif walking_tick_ < t_last_ - t_double_2
                    ref_zmp_(walking_tick_+1,1) = foot_step(current_step_num_+1,1);
                    ref_zmp_(walking_tick_+1,2) = A_;
                elseif (walking_tick_ <= t_last_) && (walking_tick_ >= t_last_ - t_double_2)
                    ref_zmp_(walking_tick_+1,1) = foot_step(current_step_num_+1,1); 
                    ref_zmp_(walking_tick_+1,2) = (Ky_ / t_double_2) * (t_total_ - (walking_tick_ - t_start_));   
                end
            else % final_step_length ~= 0
            end
        end
    elseif walking_tick_ > t_last_
        ref_zmp_(walking_tick_+1,1) = foot_step(current_step_num_+1,1);
        ref_zmp_(walking_tick_+1,2) = 0;   
    end
end


%% COM trajectory
t_start_ = t_temp_ + 1;
t_last_ = t_temp_ + t_total_;
current_step_num_ = 0;

for walking_tick_ = 0:total_tick-1
    walking_tick(walking_tick_+1) = walking_tick_;
    if walking_tick_ <= t_temp_
        if walking_tick_ < 0.5 * hz_
            ref_com_(walking_tick_+1,1) = T_LF(1,4);
            ref_com_(walking_tick_+1,2) = 0;
        elseif walking_tick_ < 1.5 * hz_
            del_x = walking_tick_ - 0.5 * hz_;
            ref_com_(walking_tick_+1,1) = T_LF(1,4);
            ref_com_(walking_tick_+1,2) = 0;
        else
            ref_com_(walking_tick_+1,1) = T_LF(1,4);
            ref_com_(walking_tick_+1,2) = 0;
        end
    elseif walking_tick_ <= t_last_
      
        if current_step_num_ == 0
            A_ = foot_step(current_step_num_+1,2);
            B_ = (foot_step(current_step_num_+1,1) + foot_step(current_step_num_+2,1)) / 2;
            Kx_ = (B_ * t_double * wn_) / (t_double * wn_ + tanh(wn_ * (t_total_t/2 - t_double)));
            Ky_ = A_ * t_double * wn_ * tanh(wn_ * (t_total_t/2 - t_double)) / (1 + t_double * wn_ * tanh(wn_ * (t_total_t/2 - t_double)));
            Cx1_ = Kx_ - B_;
            Cx2_ = Kx_ / (wn_ * t_double);
            Cy1_ = Ky_ - A_;
            Cy2_ = Ky_ / (wn_ * t_double);
            if walking_tick_ < t_start_ + t_double_1
                ref_com_(walking_tick_+1,1) = foot_step(current_step_num_+1,1);
                ref_com_(walking_tick_+1,2) = Ky_ / t_double_1 * (walking_tick_- t_start_);
            elseif (walking_tick_ < t_last_ - t_double_2) && (walking_tick_ >= t_start_ + t_double_1)
                ref_com_(walking_tick_+1,1) = foot_step(current_step_num_+1,1) + Cx1_ * cosh(wn_*((walking_tick_ - t_start_)/hz_ - t_double)) + Cx2_ * sinh(wn_*((walking_tick_ - t_start_)/hz_ - t_double));
                ref_com_(walking_tick_+1,2) = Cy1_ * cosh(wn_*((walking_tick_ - t_start_)/hz_ - t_double)) + Cy2_ * sinh(wn_*((walking_tick_ - t_start_)/hz_ - t_double)) + A_;
            elseif (walking_tick_ <= t_last_) && (walking_tick_ >= t_last_ - t_double_2)
                ref_com_(walking_tick_+1,1) = (B_ - Kx_) + (Kx_ / t_double_2) * (walking_tick_ - t_start_ - (t_total_ - t_double_2)); 
                ref_com_(walking_tick_+1,2) = (Ky_ / t_double_2) * (t_total_ - (walking_tick_ - t_start_));
            end
            if walking_tick_ == t_last_
                current_step_num_ = current_step_num_ + 1;
                t_start_ = t_start_ + t_total_;
                t_last_ = t_last_ + t_total_;
            end

        elseif current_step_num_ < foot_step_number 
            A_ = foot_step(current_step_num_+1,2);
            B_ = (foot_step(current_step_num_,1) + foot_step(current_step_num_+1,1)) / 2 - foot_step(current_step_num_,1);
            Kx_ = (B_ * t_double * wn_) / (t_double * wn_ + tanh(wn_ * (t_total_t/2 - t_double)));
            Ky_ = A_ * t_double * wn_ * tanh(wn_ * (t_total_t/2 - t_double)) / (1 + t_double * wn_ * tanh(wn_ * (t_total_t/2 - t_double)));
            Cx1_ = Kx_ - B_;
            Cx2_ = Kx_ / (wn_ * t_double);
            Cy1_ = Ky_ - A_;
            Cy2_ = Ky_ / (wn_ * t_double);

            if walking_tick_ < t_start_ + t_double_1
                ref_com_(walking_tick_+1,1) = (foot_step(current_step_num_,1) + foot_step(current_step_num_+1,1)) / 2 + (Kx_ / t_double_1) * (walking_tick_ - t_start_);
                ref_com_(walking_tick_+1,2) = Ky_ / t_double_1 * (walking_tick_ - t_start_);
            elseif walking_tick_ < t_last_ - t_double_2
                ref_com_(walking_tick_+1,1) = (foot_step(current_step_num_,1) + foot_step(current_step_num_+1,1)) / 2 + Cx1_ * cosh(wn_*((walking_tick_ - t_start_)/hz_ - t_double)) + Cx2_ * sinh(wn_*((walking_tick_ - t_start_)/hz_ - t_double)) + B_;
                ref_com_(walking_tick_+1,2) = Cy1_ * cosh(wn_*((walking_tick_ - t_start_)/hz_ - t_double)) + Cy2_ * sinh(wn_*((walking_tick_ - t_start_)/hz_ - t_double)) + A_;
            else
                ref_com_(walking_tick_+1,1) = ((foot_step(current_step_num_+1,1) + foot_step(current_step_num_+2,1)) / 2 - Kx_) + (Kx_ / t_double_2) * (walking_tick_ - t_start_ - (t_total_ - t_double_2)); 
                ref_com_(walking_tick_+1,2) = (Ky_ / t_double_2) * (t_total_ - (walking_tick_ - t_start_));
            end
            if walking_tick_ == t_last_
                current_step_num_ = current_step_num_ + 1;
                t_start_ = t_start_ + t_total_;
                t_last_ = t_last_ + t_total_;
            end
        elseif current_step_num_ == foot_step_number
            if final_step_length == 0
                A_ = foot_step(current_step_num_+1,2);
                B_ = (foot_step(current_step_num_,1) + foot_step(current_step_num_+1,1)) / 2 - foot_step(current_step_num_,1);
                Kx_ = (B_ * t_double * wn_) / (t_double * wn_ + tanh(wn_ * (t_total_t/2 - t_double)));
                Ky_ = A_ * t_double * wn_ * tanh(wn_ * (t_total_t/2 - t_double)) / (1 + t_double * wn_ * tanh(wn_ * (t_total_t/2 - t_double)));
                Cx1_ = Kx_ - B_;
                Cx2_ = Kx_ / (wn_ * t_double);
                Cy1_ = Ky_ - A_;
                Cy2_ = Ky_ / (wn_ * t_double);                
                if walking_tick_ < t_start_ + t_double_1
                    ref_com_(walking_tick_+1,1) = (foot_step(current_step_num_,1) + foot_step(current_step_num_+1,1)) / 2 + (Kx_ / t_double_1) * (walking_tick_ - t_start_);
                    ref_com_(walking_tick_+1,2) = Ky_ / t_double_1 * (walking_tick_- t_start_);
                elseif walking_tick_ < t_last_ - t_double_2
                    ref_com_(walking_tick_+1,1) = foot_step(current_step_num_+1,1);
                    ref_com_(walking_tick_+1,2) = Cy1_ * cosh(wn_*((walking_tick_ - t_start_)/hz_ - t_double)) + Cy2_ * sinh(wn_*((walking_tick_ - t_start_)/hz_ - t_double)) + A_;
                elseif (walking_tick_ <= t_last_) && (walking_tick_ >= t_last_ - t_double_2)
                    ref_com_(walking_tick_+1,1) = foot_step(current_step_num_+1,1); 
                    ref_com_(walking_tick_+1,2) = (Ky_ / t_double_2) * (t_total_ - (walking_tick_ - t_start_));   
                end
            else % final_step_length ~= 0
            end
        end
    elseif walking_tick_ > t_last_
        ref_com_(walking_tick_+1,1) = foot_step(current_step_num_+1,1);
        ref_com_(walking_tick_+1,2) = 0;   
    end
end


%% Foot trajectory


y = [0 0.2];

cs = spline([-4 4],[0 y 0]);
xx = linspace(-4,4,101);
yy = ppval(cs,xx);
plot(xx,yy);

t_start_ = t_temp_ + 1;
t_last_ = t_temp_ + t_total_;
current_step_num_ = 0;

rfoot = zeros(total_tick,3);
lfoot = zeros(total_tick,3);


for walking_tick_ = 0:total_tick-1
    walking_tick(walking_tick_+1) = walking_tick_;
    if walking_tick_ <= t_temp_

        PHASE_variable(walking_tick_+1) = 1; %DSP

        if walking_tick_ < 0.5 * hz_
            lfoot(walking_tick_+1,1) = T_LF(1,4);
            lfoot(walking_tick_+1,2) = T_LF(2,4);
            lfoot(walking_tick_+1,3) = T_LF(3,4);

            rfoot(walking_tick_+1,1) = T_RF(1,4);
            rfoot(walking_tick_+1,2) = T_RF(2,4);
            rfoot(walking_tick_+1,3) = T_RF(3,4);
        elseif walking_tick_ < 1.5 * hz_
            del_x = walking_tick_ - 0.5 * hz_;
            lfoot(walking_tick_+1,1) = T_LF(1,4);
            lfoot(walking_tick_+1,2) = T_LF(2,4);
            lfoot(walking_tick_+1,3) = T_LF(3,4);
            
            rfoot(walking_tick_+1,1) = T_RF(1,4);
            rfoot(walking_tick_+1,2) = T_RF(2,4);
            rfoot(walking_tick_+1,3) = T_RF(3,4);
        else
            lfoot(walking_tick_+1,1) = T_LF(1,4);
            lfoot(walking_tick_+1,2) = T_LF(2,4);
            lfoot(walking_tick_+1,3) = T_LF(3,4);

            rfoot(walking_tick_+1,1) = T_RF(1,4);
            rfoot(walking_tick_+1,2) = T_RF(2,4);
            rfoot(walking_tick_+1,3) = T_RF(3,4);
        end
    elseif walking_tick_ <= t_last_
        if current_step_num_ == 0
            
            if walking_tick_ < t_start_ + t_double_1
                PHASE_variable(walking_tick_+1) = 1; %DSP
            elseif (walking_tick_ < t_last_ - t_double_2) && (walking_tick_ >= t_start_ + t_double_1)
                PHASE_variable(walking_tick_+1) = 3; %SSP
            elseif (walking_tick_ <= t_last_) && (walking_tick_ >= t_last_ - t_double_2)
                PHASE_variable(walking_tick_+1) = 1; %DSP
            end
            lfoot(walking_tick_+1,1) = T_LF(1,4);
            lfoot(walking_tick_+1,2) = T_LF(2,4);
            lfoot(walking_tick_+1,3) = T_LF(3,4);

            RF = [foot_step(current_step_num_+1,1) foot_step(current_step_num_+2,1)];
            rf_y = spline([t_start_ t_last_],[0 RF 0]);
            rf_x = linspace(t_start_,t_last_,t_last_-t_start_+1);
            prf = ppval(rf_y,rf_x);         

            RF1 = [T_RF(3,4) T_RF(3,4)+foot_height_ T_RF(3,4)];
            rf_y1 = spline([t_start_ (t_start_+t_last_)/2 t_last_],[0 RF1 0]);
            rf_x1 = linspace(t_start_,t_last_,t_last_-t_start_+1);
            prf1 = ppval(rf_y1,rf_x1);

            rfoot(walking_tick_+1,1) = prf(walking_tick_+1-t_start_);
            rfoot(walking_tick_+1,2) = T_RF(2,4);
            rfoot(walking_tick_+1,3) = prf1(walking_tick_+1-t_start_);
            
            if walking_tick_ == t_last_
                current_step_num_ = current_step_num_ + 1;
                t_start_ = t_start_ + t_total_;
                t_last_ = t_last_ + t_total_;
            end
        elseif current_step_num_ < foot_step_number
            if mod(current_step_num_,2) == 1
                if walking_tick_ < t_start_ + t_double_1
                    PHASE_variable(walking_tick_+1) = 1; %DSP
                elseif (walking_tick_ < t_last_ - t_double_2) && (walking_tick_ >= t_start_ + t_double_1)
                    PHASE_variable(walking_tick_+1) = 2; %SSP
                elseif (walking_tick_ <= t_last_) && (walking_tick_ >= t_last_ - t_double_2)
                    PHASE_variable(walking_tick_+1) = 1; %DSP
                end                
                LF = [foot_step(current_step_num_,1) foot_step(current_step_num_+2,1)];
                lf_y = spline([t_start_ t_last_],[0 LF 0]);
                lf_x = linspace(t_start_,t_last_,t_last_-t_start_+1);
                plf = ppval(lf_y,lf_x);

                LF1 = [T_LF(3,4) T_LF(3,4)+foot_height_ T_LF(3,4)];
                lf_y1 = spline([t_start_ (t_start_+t_last_)/2 t_last_],[0 LF1 0]);
                lf_x1 = linspace(t_start_,t_last_,t_last_-t_start_+1);
                plf1 = ppval(lf_y1,lf_x1);
                
                lfoot(walking_tick_+1,1) = plf(walking_tick_+1-t_start_);
    
                rfoot(walking_tick_+1,1) = foot_step(current_step_num_+1,1);

                lfoot(walking_tick_+1,2) = T_LF(2,4);
                rfoot(walking_tick_+1,2) = T_RF(2,4);

                lfoot(walking_tick_+1,3) = plf1(walking_tick_+1-t_start_);
                rfoot(walking_tick_+1,3) = T_RF(3,4);
    
                if walking_tick_ == t_last_
                    current_step_num_ = current_step_num_ + 1;
                    t_start_ = t_start_ + t_total_;
                    t_last_ = t_last_ + t_total_;
                end
    
            elseif mod(current_step_num_,2) == 0 && current_step_num_ > 0
                if walking_tick_ < t_start_ + t_double_1
                    PHASE_variable(walking_tick_+1) = 1; %DSP
                elseif (walking_tick_ < t_last_ - t_double_2) && (walking_tick_ >= t_start_ + t_double_1)
                    PHASE_variable(walking_tick_+1) = 3; %SSP
                elseif (walking_tick_ <= t_last_) && (walking_tick_ >= t_last_ - t_double_2)
                    PHASE_variable(walking_tick_+1) = 1; %DSP
                end    
                RF = [foot_step(current_step_num_,1) foot_step(current_step_num_+2,1)];
                rf_y = spline([t_start_ t_last_],[0 RF 0]);
                rf_x = linspace(t_start_,t_last_,t_last_-t_start_+1);
                prf = ppval(rf_y,rf_x);
                
                RF1 = [T_RF(3,4) T_RF(3,4)+foot_height_ T_RF(3,4)];
                rf_y1 = spline([t_start_ (t_start_+t_last_)/2 t_last_],[0 RF1 0]);
                rf_x1 = linspace(t_start_,t_last_,t_last_-t_start_+1);
                prf1 = ppval(rf_y1,rf_x1);

                rfoot(walking_tick_+1,1) = prf(walking_tick_+1-t_start_);
    
                lfoot(walking_tick_+1,1) = foot_step(current_step_num_+1,1);

                lfoot(walking_tick_+1,2) = T_LF(2,4);
                rfoot(walking_tick_+1,2) = T_RF(2,4);

                rfoot(walking_tick_+1,3) = prf1(walking_tick_+1-t_start_);
                lfoot(walking_tick_+1,3) = T_LF(3,4);
    
                if walking_tick_ == t_last_
                    current_step_num_ = current_step_num_ + 1;
                    t_start_ = t_start_ + t_total_;
                    t_last_ = t_last_ + t_total_;
                end
            end
        elseif current_step_num_ == foot_step_number
            if final_step_length == 0                       
                if mod(current_step_num_,2) == 1
                    if walking_tick_ < t_start_ + t_double_1
                        PHASE_variable(walking_tick_+1) = 1; %DSP
                    elseif (walking_tick_ < t_last_ - t_double_2) && (walking_tick_ >= t_start_ + t_double_1)
                        PHASE_variable(walking_tick_+1) = 2; %SSP
                    elseif (walking_tick_ <= t_last_) && (walking_tick_ >= t_last_ - t_double_2)
                        PHASE_variable(walking_tick_+1) = 1; %DSP
                    end    
                    LF = [foot_step(current_step_num_,1) foot_step(current_step_num_+1,1)];
                    lf_y = spline([t_start_ t_last_],[0 LF 0]);
                    lf_x = linspace(t_start_,t_last_,t_last_-t_start_+1);
                    plf = ppval(lf_y,lf_x);

                    LF1 = [T_LF(3,4) T_LF(3,4)+foot_height_ T_LF(3,4)];
                    lf_y1 = spline([t_start_ (t_start_+t_last_)/2 t_last_],[0 LF1 0]);
                    lf_x1 = linspace(t_start_,t_last_,t_last_-t_start_+1);
                    plf1 = ppval(lf_y1,lf_x1);
                    
                    lfoot(walking_tick_+1,1) = plf(walking_tick_+1-t_start_);
                    rfoot(walking_tick_+1,1) = foot_step(current_step_num_+1,1);

                    lfoot(walking_tick_+1,2) = T_LF(2,4);
                    rfoot(walking_tick_+1,2) = T_RF(2,4);
        
                    lfoot(walking_tick_+1,3) = plf1(walking_tick_+1-t_start_);
                    rfoot(walking_tick_+1,3) = T_RF(3,4);
                elseif mod(current_step_num_,2) == 0 && current_step_num_ > 0
                    if walking_tick_ < t_start_ + t_double_1
                        PHASE_variable(walking_tick_+1) = 1; %DSP
                    elseif (walking_tick_ < t_last_ - t_double_2) && (walking_tick_ >= t_start_ + t_double_1)
                        PHASE_variable(walking_tick_+1) = 3; %SSP
                    elseif (walking_tick_ <= t_last_) && (walking_tick_ >= t_last_ - t_double_2)
                        PHASE_variable(walking_tick_+1) = 1; %DSP
                    end        
                    RF = [foot_step(current_step_num_,1) foot_step(current_step_num_+1,1)];
                    rf_y = spline([t_start_ t_last_],[0 RF 0]);
                    rf_x = linspace(t_start_,t_last_,t_last_-t_start_+1);
                    prf = ppval(rf_y,rf_x);

                    RF1 = [T_RF(3,4) T_RF(3,4)+foot_height_ T_RF(3,4)];
                    rf_y1 = spline([t_start_ (t_start_+t_last_)/2 t_last_],[0 RF1 0]);
                    rf_x1 = linspace(t_start_,t_last_,t_last_-t_start_+1);
                    prf1 = ppval(rf_y1,rf_x1);
                    
                    rfoot(walking_tick_+1,1) = prf(walking_tick_+1-t_start_);
                    lfoot(walking_tick_+1,1) = foot_step(current_step_num_+1,1);

                    lfoot(walking_tick_+1,2) = T_LF(2,4);
                    rfoot(walking_tick_+1,2) = T_RF(2,4);

                    rfoot(walking_tick_+1,3) = prf1(walking_tick_+1-t_start_);
                    lfoot(walking_tick_+1,3) = T_LF(3,4);
                end
            else % final_step_length ~= 0
            end
        end
    elseif walking_tick_ > t_last_
        PHASE_variable(walking_tick_+1) = 1; %DSP
        lfoot(walking_tick_+1,1) = foot_step(current_step_num_+1,1);
        lfoot(walking_tick_+1,2) = T_LF(2,4);
        lfoot(walking_tick_+1,3) = T_LF(3,4);

        rfoot(walking_tick_+1,1) = foot_step(current_step_num_+1,1);
        rfoot(walking_tick_+1,2) = T_RF(2,4);
        rfoot(walking_tick_+1,3) = T_RF(3,4);
    end
end


subplot(3,1,1)
plot(walking_tick,ref_zmp_(:,1),walking_tick,ref_com_(:,1),walking_tick,lfoot(:,1),walking_tick,rfoot(:,1))
legend('zmp','com','left','right')
subplot(3,1,2)
plot(walking_tick,ref_zmp_(:,2),walking_tick,ref_com_(:,2),walking_tick,lfoot(:,2),walking_tick,rfoot(:,2))
subplot(3,1,3)
plot(walking_tick,lfoot(:,3),walking_tick,rfoot(:,3),walking_tick,PHASE_variable)
legend('left','right')

