%% NI&I HW0 2023”
% Siyu Wang
%% load data 
data_load = load("HW0 NI&I 2023 DATA.mat");
H = data_load.H;
L = data_load.L;
numDinY = data_load.numDinY;
numDinM = data_load.numDinM;
numMinY  = data_load.numMinY;


%% Problem I: plot the H and L temperture information 
DaysInAYear = 1:1:numDinY; % constructing each day in the year

% declear a figure
figure();
plot(DaysInAYear,H,"r"); % plot H data
hold on
plot(DaysInAYear,L,"b"); % plot L data
hold off 

xlim([1,365]);
ylim([0,100]); % set y axis limitations
title("Average Daily High and Low Temperatures in Baltimore"); % set figure title
legend("High","Low"); % set legend 
xlabel("Day of the Year"); % set x label 
ylabel("Temperature (deg F)"); % set y label


%% Problem II: Mean average daily low temp for each month
current_date = 1;
for i = 1:1:numMinY
    meanOfMonthL = mean(L(current_date:current_date+numDinM(i)-1))
    current_date = current_date + numDinM(i);
end 


%% Problem III: Finding values in an Array
% 1)
dayHAbove = find(H>=50);
numOfMonthHAbove = length(dayHAbove);
disp("Q3 1) There are " + numOfMonthHAbove +" days in a year with an average daily high temperature greater than or equal to 50 degrees F")

% 2)
maxTemp = max(H);
dayHMax = find(H==maxTemp);
numOfDaysMaxHTemp = length(dayHMax);
disp("Q3 2) There are " + numOfDaysMaxHTemp +" day(s) of the year does the hottest/highest average daily high temperature occur, they are: ")
disp(dayHMax)

