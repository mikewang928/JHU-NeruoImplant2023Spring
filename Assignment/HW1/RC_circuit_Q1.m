%%
[TimeAllList, VoutAllList] = RC_circuit(1000,0.00025);
plot(TimeAllList,VoutAllList);

%%
function [Time_all_list,Vout_all_list]= RC_circuit(R,C)
    % 0v at t = 0s
    % steps up to 1v at t = 0.1s
    % steps backdown to 0v at 0.8s
    % ends at 1.5s 
    Time_all_list = [];
    Vout_all_list = [];
    Time_current_list = [];
    Vout_current_list = [];

    Vin = 0;
    Vout_prev = 0;
    Time_prev=0;
    [Time_current_list, Time_prev, Vout_current_list, Vout_prev] = RC_step(R,C,Vout_prev,Vin, Time_prev,0.1);
    Time_all_list = [Time_all_list,Time_current_list];
    Vout_all_list = [Vout_all_list, Vout_current_list];
    
    Vin = 1;
    [Time_current_list, Time_prev, Vout_current_list, Vout_prev] = RC_step(R,C,Vout_prev,Vin, Time_prev,0.7);
    Time_all_list = [Time_all_list,Time_current_list];
    Vout_all_list = [Vout_all_list, Vout_current_list];
    
    Vin = 0;
    [Time_current_list, Time_prev, Vout_current_list, Vout_prev] = RC_step(R,C,Vout_prev,Vin, Time_prev,0.7);
    Time_all_list = [Time_all_list,Time_current_list];
    Vout_all_list = [Vout_all_list, Vout_current_list];
end 

% This function models the Vout in a step down step voltage
function [Time_list, Time_last, Vout_list, Vout]=RC_step(R,C,Vout_prev, Vin, Time_prev, duration)
    dt = 0.001;    
    Vout = Vout_prev; %initilazation
    Time_list = [];
    Vout_list = [];
    for t=0:dt:duration
        dVout = dt*((Vin-Vout)/(R*C));
        Vout = Vout + dVout;
        Time_list = [Time_list,t+Time_prev];
        Time_last = t+Time_prev;
        Vout_list = [Vout_list, Vout];
    end 
end 



