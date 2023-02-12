%% Homework I
% Siyu Wang (swang333)
%% Question 1.4 R-C circuit on capcitor simulation 
[TimeAllList, VoutAllList, iAllList, VinAllList, VtList] = RC_circuit_Q1(1000,0.00025);
figure();
plot(TimeAllList,VoutAllList);
hold on 
plot(TimeAllList,VinAllList);
hold on 
plot(0.35, VtList(1), 'xr', 'MarkerSize',20)
hold on
plot(1.05, VtList(2), 'xr', 'MarkerSize',20)
hold on
title("Voltage accross capacitor in a RC circuit"); % set figure title
legend("Voltage out","Voltage in", "voltage at rising time constant", "voltage at falling time constant"); % set legend 
xlabel("Simulation Time (s)"); % set x label 
ylabel("Voltage (V)"); % set y label

figure();
plot(TimeAllList,iAllList);
title("Current in a RC circuit"); % set figure title
xlabel("Simulation Time (s)"); % set x label 
ylabel("Current (A)"); % set y label

%% Question 2.2 R-C circuit on resistor simulation
[TimeAllList, VoutAllList, VinAllList, iAllList, VtList] = RC_circuit_Q2(1000,0.00025);
figure();
plot(TimeAllList,VoutAllList);
hold on 
plot(TimeAllList,VinAllList);
hold on 
plot(0.35, VtList(1), 'xr', 'MarkerSize',20)
hold on
plot(1.85, VtList(2), 'xr', 'MarkerSize',20)
hold on
title("Voltage accross Resistor in a RC circuit"); % set figure title
legend("Voltage out","Voltage in", "voltage at rising time constant", "voltage at falling time constant"); % set legend
xlabel("Simulation Time (s)"); % set x label 
ylabel("Voltage (V)"); % set y label

figure();
plot(TimeAllList,iAllList);
title("Current in a RC circuit in Question 2"); % set figure title
xlabel("Simulation Time (s)"); % set x label 
ylabel("Current (A)"); % set y label


%% Question 4.3 RCL circuit Band-stop filter 
L = 0.0101;
C = 1e-9;
R = 1273.24;
[f_list,angle_list,mag_list] = RCL_circuit_Q4(R,C,L);
figure();
semilogx(f_list,mag_list);
title("Bold Plot"); % set figure title
legend("Magnitude"); % set legend
xlabel("Frequency(Hz)");
ylabel("|H|");
grid on;
figure();
semilogx(f_list,angle_list);
title("Bold Plot"); % set figure title
legend("Angle"); % set legend
xlabel("Frequency(Hz)");
ylabel("<H (Degrees)");
grid on;

%% Question 4.5 RCL circuit freq @ 40000 Hz
figure();
semilogx(f_list,mag_list);
hold on
semilogx(40000,0.749,'xr', 'MarkerSize',20);
hold on
title("Bold Plot"); % set figure title
legend("Magnitude"); % set legend
xlabel("Frequency(Hz)");
ylabel("|H|");
grid on;
figure();
semilogx(f_list,angle_list);
hold on
semilogx(40000,-42.35,'xr', 'MarkerSize',20);
hold on
title("Bold Plot"); % set figure title
legend("Angle"); % set legend
xlabel("Frequency(Hz)");
ylabel("<H (Degrees)");
grid on;

%%
%-------------------- Question 1 helper function --------------
function [Time_all_list,Vout_all_list, i_all_list, Vin_all_list, Vt_list]= RC_circuit_Q1(R,C)
    % 0v at t = 0s
    % steps up to 1v at t = 0.1s
    % steps backdown to 0v at 1.6s
    % ends at 3.1s 
    Time_all_list = [];
    Vout_all_list = [];
    Vin_all_list = [];
    i_all_list = [];
    Time_current_list = [];
    Vout_current_list = [];
    Vt_list = [];

    Vin = 0;
    Vout_prev = 0;
    Time_prev=0;
    [Time_current_list, i_current_list, Vin_current_list, Time_prev, Vout_current_list, Vout_prev] = RC_step_C(R,C,Vout_prev,Vin, Time_prev,0.1);
    Time_all_list = [Time_all_list,Time_current_list];
    Vout_all_list = [Vout_all_list, Vout_current_list];
    i_all_list = [i_all_list,i_current_list];
    Vin_all_list = [Vin_all_list,Vin_current_list];
    
    Vin = 1;
    [Time_current_list, i_current_list, Vin_current_list, Time_prev, Vout_current_list, Vout_prev, Voltage_at_TimeConstant] = RC_step_C(R,C,Vout_prev,Vin, Time_prev,0.7);
    Time_all_list = [Time_all_list,Time_current_list];
    Vout_all_list = [Vout_all_list, Vout_current_list];
    i_all_list = [i_all_list,i_current_list];
    Vin_all_list = [Vin_all_list,Vin_current_list];
    Vt_list = [Vt_list,Voltage_at_TimeConstant];

    Vin = 0;
    [Time_current_list, i_current_list, Vin_current_list, Time_prev, Vout_current_list, Vout_prev, Voltage_at_TimeConstant] = RC_step_C(R,C,Vout_prev,Vin, Time_prev,0.7);
    Time_all_list = [Time_all_list,Time_current_list];
    Vout_all_list = [Vout_all_list, Vout_current_list];
    i_all_list = [i_all_list,i_current_list];
    Vin_all_list = [Vin_all_list,Vin_current_list];
    Vt_list = [Vt_list,Voltage_at_TimeConstant];
end 

% This function models the Vout in a step voltage
function [Time_list, i_list, Vin_list, Time_last, Vout_list, Vout, Voltage_at_TimeConstant]=RC_step_C(R,C,Vout_prev, Vin, Time_prev, duration)
    dt = 0.001;    
    Vout = Vout_prev; %initilazation
    Time_list = [];
    Vout_list = [];
    Vin_list = [];
    i_list = [];
    for t=0:dt:duration
        dVout = dt*((Vin-Vout)/(R*C));
        i = C*(dVout/dt);
        Vout = Vout + dVout;
        Time_list = [Time_list,t+Time_prev];
        Time_last = t+Time_prev;
        Vout_list = [Vout_list, Vout];
        i_list = [i_list,i];
        Vin_list = [Vin_list,Vin];
        if t == 0.25
            Voltage_at_TimeConstant = Vout;
        end 
    end 
end 

%-------------------- Question 2 helper function --------------
function [Time_all_list,Vout_all_list, Vin_all_list, i_all_list, Vt_list]= RC_circuit_Q2(R,C)
    % 0v at t = 0s
    % steps up to 1v at t = 0.1s
    % steps backdown to 0v at 0.8s
    % ends at 1.5s 
    Time_all_list = [];
    Vout_all_list = [];
    Vin_all_list = [];
    i_all_list = [];
    Time_current_list = [];
    Vout_current_list = [];
    Vin_current_list = []
    Vt_list = [];

    Vin = 0;
    Vout_prev = 0;
    Vc_prev = 0;
    Time_prev=0;
    [Time_current_list, i_current_list, Vin_current_list, Time_prev, Vout_current_list, Vout_prev, Vc_prev] = RC_step_R(R,C,Vout_prev,Vc_prev, Vin,0, Time_prev,0.1);
    Time_all_list = [Time_all_list,Time_current_list];
    Vout_all_list = [Vout_all_list, Vout_current_list];
    Vin_all_list = [Vin_all_list,Vin_current_list];
    i_all_list = [i_all_list,i_current_list];
    
    Vin = 1;
    [Time_current_list, i_current_list, Vin_current_list, Time_prev, Vout_current_list, Vout_prev, Vc_prev, Voltage_at_TimeConstant] = RC_step_R(R,C,Vout_prev,Vc_prev, Vin,1, Time_prev,1.5);
    Time_all_list = [Time_all_list,Time_current_list];
    Vout_all_list = [Vout_all_list, Vout_current_list];
    Vin_all_list = [Vin_all_list,Vin_current_list];
    i_all_list = [i_all_list,i_current_list];
    Vt_list = [Vt_list,Voltage_at_TimeConstant];

    Vin = 0;
    [Time_current_list, i_current_list, Vin_current_list, Time_prev, Vout_current_list, Vout_prev, Vc_prev, Voltage_at_TimeConstant] = RC_step_R(R,C,Vout_prev,Vc_prev, Vin,-1, Time_prev,1.5);
    Time_all_list = [Time_all_list,Time_current_list];
    Vout_all_list = [Vout_all_list, Vout_current_list];
    Vin_all_list = [Vin_all_list,Vin_current_list];
    i_all_list = [i_all_list,i_current_list];
    Vt_list = [Vt_list,Voltage_at_TimeConstant];
end 

% This function models the Vout in a step voltage
function [Time_list, i_list, Vin_list, Time_last, Vout_list, Vout, Vc, Voltage_at_TimeConstant]=RC_step_R(R,C,Vout_prev, Vc_prev, Vin, dVin, Time_prev, duration)
    dt = 0.001;    
    Vout = Vout_prev; %initilazation
    Vc = Vc_prev;
    Time_list = [];
    Vout_list = [];
    Vin_list = [];
    i_list = [];
    for t=0:dt:duration
        if t==0
            dVout = dVin-dt*(Vout/(R*C));
        else
            dVout = -dt*(Vout/(R*C));
        end
        dVc = dt*((Vin-Vc)/(R*C));
        i = C*(dVc/dt);
        Vout = Vout + dVout;
        Vc = Vc + dVc;
        Time_list = [Time_list,t+Time_prev];
        Time_last = t+Time_prev;
        Vout_list = [Vout_list, Vout];
        i_list = [i_list,i];
        Vin_list = [Vin_list,Vin];
        if t == 0.25
            Voltage_at_TimeConstant = Vout;
        end 
    end 
end 

%-------------------- Question 4 helper function --------------
function [f_list,angle_list, mag_list]= RCL_circuit_Q4(R,C,L)
    mag_list = [];
    angle_list = [];
    f_list = [];
    for f_exp=2:0.001:8
        f = 10^f_exp;
        w = 2*pi*f;
        H = (1-L*C*w^2)/(complex(1-L*C*w^2,w*C*R));
        mag = abs(H);
        phase = angle(H)*180/pi;
        f_list = [f_list,f];
        mag_list = [mag_list, mag];
        angle_list = [angle_list,phase];
    end
end 
