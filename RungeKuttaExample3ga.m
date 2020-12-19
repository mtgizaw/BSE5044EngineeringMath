%RungeKuttaExample3ga.m

%Written by Ryan Senger, 6/10/2020

%this script simulates the following coupled ODE's:
%d[X]/dt = u[X]
%u=(umax*[S])/(Ks+[S])
%d[S]/dt = -1/Yxs*d[X]/dt
 
%this script also compares simulation results to raw data.

%this script was built from RungeKuttaExample3mse.m and does NOT clear
%variables, define umax Ks Yxs, plot, or print to the screen.  

%this version is called by GeneticAlgorith.m

umax = chrom(m,1);
Ks = chrom(m,2);
 
%clear
%clc
 
%specify the raw data here
data(:,1)=[0; 2; 4; 8; 10; 12; 14; 16; 18]; %time
data(:,2)=[0.2; 0.211; 0.305; 0.98; 1.77; 3.2; 5.6; 6.15; 6.2]; %[X]
data(:,3)=[9.23; 9.21; 9.07; 8.03; 6.8; 4.6; 0.92; 0.077; 0]; %[S]
 
%specify the kinetic rate constants
% umax=0.394; %guesses from llsa
% Ks=1.83;
Yxs=0.65;
 
%define the initial values and time-step
h=0.1;
t0=0;
x0=0.2;
s0=9.23;
 
%define the length of the simulation
tfinal=18;
 
%calculate the time vector
t=(t0:h:tfinal);

%calculate the number of time steps
steps=size(t,2);
 
%initiate answer vectors
x=zeros(1,size(t,2));
s=zeros(1,size(t,2));
 
%specify initial values
x(1,1)=x0;
s(1,1)=s0;
 
%start the Runge-Kutta loop
for i=1:steps-1
    
    %the k1 value
    T=t(1,i);
    X=x(1,i);
    S=s(1,i);
    
    k1x=(umax*S*X)/(Ks+S);
    k1s=-1/Yxs*k1x;
    
    %the k2 value
    T=t(1,i)+0.5*h;
    X=x(1,i)+0.5*h*k1x;
    S=s(1,i)+0.5*h*k1s;
 
    k2x=(umax*S*X)/(Ks+S);
    k2s=-1/Yxs*k2x;
    
    %the k3 value
    T=t(1,i)+0.5*h;
    X=x(1,i)+0.5*h*k2x;
    S=s(1,i)+0.5*h*k2s;
 
    k3x=(umax*S*X)/(Ks+S);
    k3s=-1/Yxs*k3x;
    
    %the k4 value
    T=t(1,i)+h;
    X=x(1,i)+h*k3x;
    S=s(1,i)+h*k3s;
 
    k4x=(umax*S*X)/(Ks+S);
    k4s=-1/Yxs*k4x;
    
    %calculate the new values
    x(1,i+1)=x(1,i)+h/6*(k1x+2*k2x+2*k3x+k4x);
    s(1,i+1)=s(1,i)+h/6*(k1s+2*k2s+2*k3s+k4s);
    
end

%calculate the mse between raw data and simulation values
data(:,4) = 0;  %this column will be the simulation value of x
data(:,5) = 0;  %this column will be the simulation value of s
data(:,6) = 0;  %this is the squared error in x
data(:,7) = 0;  %this is the squared error in s
data(:,8) = 0;  %this is the total error in x and s

%searching for all raw time points of the data matrix in the "t" vector
for i = 1:size(data,1)
    
    %finding the column in t
    id = find(t==data(i,1));
    
    %extracting x and s values
    if isempty(id) ~= 1
        
        data(i,4) = x(1,id);    %the value of x
        data(i,5) = s(1,id);    %the value of s
        
    else
        
        fprintf(['Simulation data points do not exist for all raw data points.' '\n']);
        return
        
    end
    
end

%calculate the squared error values
data(:,6) = (data(:,2)-data(:,4)).^2;
data(:,7) = (data(:,3)-data(:,5)).^2;
data(:,8) = data(:,6) + data(:,7);

%calculate the mse
MeanSquaredError = sum(data(:,8))./size(data,1);
result = MeanSquaredError;

% %plotting the raw data points as circles
% subplot(2,1,1)
% plot(data(:,1),data(:,2),'o');
% xlabel(['Time']);
% ylabel(['[X]']);
% hold on
% subplot(2,1,2)
% plot(data(:,1),data(:,3),'o');
% xlabel(['Time']);
% ylabel(['[S]']);
% hold on
%  
% %plotting simulation results
% subplot(2,1,1)
% plot(t,x,'b')           %plots as blue line
% % plot(t,x,'r')     %plots red line
% 
% subplot(2,1,2)
% plot(t,s,'b')
% % plot(t,s,'r')
% 
% %print the value of mse to the command window
% fprintf(['Mean-squared error = ' num2str(MeanSquaredError) '\n']);
% fprintf(['Done!' '\n']);




