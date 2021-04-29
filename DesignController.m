clear all
close all
clc

% Routine to compute the gain of regulator
% by pole placement and optimal control LQG
% Olivier Verlinden and Hoai Nam Huynh (2020)
% Olivier.Verlinden@umons.ac.be
% HoaiNam.Huynh@ubc.ca

%%%%%%%%%%
% Inputs %
%%%%%%%%%%

% Load the system matrices generated from EasyDyn

% Name of system matrices
MatrixName = 'segway';

% Matrix C: Measured states
Cp=[1 0]; % Position [x theta]
Cv=[0 0]; % Velocity [x_d theta_d]

% Matrix D: influence of input on measure
D=0;

% Optimal control parameters
R = 0.0001; % State and input weight

% Discrete time step
h=0.01; 

% Pole placement parameters (conjugate poles)
zeta=[0.9 0.9]; % Desired damping
freq=[1.1 1.2]; % Desired frequency [Hz]

% Frequency of extended state
freq_extended_state = [1.15];

%----------------------------------------------------------%

%%%%%%%%%%%
% Program %
%%%%%%%%%%%
% Load the system matrices generated from EasyDyn
MM = load(strcat(MatrixName,'.mm')); % Mass Matrix
CC = load(strcat(MatrixName,'.cc')); % Damping Matrix
KK = load(strcat(MatrixName,'.kk')); % Stiffness Matrix
FF = load(strcat(MatrixName,'.ff')); % Input Matrix

% Construction of the equivalent state-space form
I = eye(size(MM));
Z = zeros(size(MM));

A=[Z I; -inv(MM)*KK -inv(MM)*CC];
B=[zeros(size(FF)); inv(MM)*FF];
C=[Cp Cv];
[nm,ns] = size(C); % Retrieve dimension of C for new state

% Construction of the augmented system with a new state
% corresponding to the integral of the error
A2=[A zeros(ns,nm);-C zeros(nm,nm)];
B2=[B; -D];
C2=[C zeros(nm,nm)];
D2=D;

% Continuous model
states= {'x' 'alpha' 'xd' 'alphad' 'interr'};
inputs= 'Vmot';
outputs= 'x';
sys_roul=ss(A2,B2,C2,D2,'statename', states,'inputname', inputs,...
            'outputname', outputs);
        
% Discrete-time model
sys_roulDT=c2d(sys_roul,h);

%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal control: LQR %
%%%%%%%%%%%%%%%%%%%%%%%%

% Designing a state-feedback controller by LQR
KLQGDT=lqr(sys_roulDT,diag(ones(1,(2*size(MM,1))+1)),R)

%%%%%%%%%%%%%%%%%%
% Pole placement %
%%%%%%%%%%%%%%%%%%

% Designing a state feedback controller by pole placement
% The roots are chosen at 1.1 Hz and 1.2 Hz with a damping of 90 % 
w=2*pi; % 1 Hz
DesiredPoles=zeros(1,2*size(MM,1));
DesiredPolesDiscrete=zeros(1,2*size(MM,1));
for i_pole = 1 : size(MM,1)
    DesiredPoles(1,1+2*(i_pole-1))=freq(1,i_pole)*(-zeta(1,i_pole)*w+1i*sqrt(1-zeta(1,i_pole)*zeta(1,i_pole))*w);
    DesiredPolesDiscrete(1,1+2*(i_pole-1))= exp( DesiredPoles(1,1+2*(i_pole-1))*h  );
    DesiredPoles(1,2*i_pole)=freq(1,i_pole)*(-zeta(1,i_pole)*w-1i*sqrt(1-zeta(1,i_pole)*zeta(1,i_pole))*w);
    DesiredPolesDiscrete(1,2*i_pole)=  exp( DesiredPoles(1,2*i_pole) *h  );
end
% Extended state
PoleExtendedState = zeros(1,size(freq_extended_state,2));
for iExtPole=1 : size(freq_extended_state,2)
    PoleExtendedState(1,iExtPole)=exp(-freq_extended_state(1,iExtPole)*w*h);
end
KPPDT=place(sys_roulDT.a,sys_roulDT.b,[DesiredPolesDiscrete(1,:) PoleExtendedState(1,:)])
