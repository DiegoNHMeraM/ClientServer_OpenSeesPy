% created by Diego Mera Muñoz
% modified by Diego Mera Muñoz (diego.mera@sansano.usm.cl) 11/2021

dtsim = 0.02;
dtfast = 1/4000;

%%
%%%%%%%%%%%%%%%%%% Substructuring and non-linear model %%%%%%%%%%%%%%%%%%
% Mass, damping and initial stiffness of the experimental sub. 
me = 0.0001;             % [lb/in/s^2]
ke = 2.9;                % [lb/in]
ce = 2*0.05*sqrt(me*ke); % [lb/in/s]

adaptive = 0; % 0 = Not adaptive; 1 = Adaptive
learningrate = 0.01;

%% Description of the parameters for non-linear models
%%%%%%%%%%%%%%%%%%%%%%%% BoucWen parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% asiv = 0 pure hysteretic ... 1 linear
% xy = Yield displacement
% nsiv = regulates curve shape
% n1siv = regulates initial form of discharge
% n2siv = 1-n1siv

%%%%%%%%%%%% Degradation parameters (mostagel models) %%%%%%%%%%%%%%%%%%%%%
% stiffdeg = Stiffness degradation
% strdegHmos = Degradation of dissipated hysteretic energy

%%%%%%%%%%%%%%%%%%%%%%%% gap-closing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kgap = kgap*ke*Ngap
% Ngap = (x-xgap)^Ngap-1
% xgap = There is only force if x> xgap

%%%%%%%%%%%%%%%%%%%%%%%% Cubic spring stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%
% khard -> Fhard=khard*x^3

%% Lineal (Case 0)
%
% BoucWen parameters
asiv=0.9999;         
xy=100;       
nsiv=40;          
n1siv=0.5;       
n2siv=1-n1siv;  

% Degradation parameters (mostagel models)
stiffdeg=0;    
strdegHmos=0;  

% gap-closing parameters
kgap=4*0;         
Ngap=1.9;       
xgap=2;   

% Cubic spring stiffness
khard=0;    
%}

%% Non-Linear (Case 1)
%{
% BoucWen parameters
asiv=0.9;          
xy=0.2;       
nsiv=40;          
n1siv=0.1;       
n2siv=1-n1siv;   

% Degradation parameters (mostagel models)
stiffdeg=0;    
strdegHmos=0;  

% gap-closing parameters
kgap=4*0;         
Ngap=1.9;       
xgap=2;   

% Cubic spring stiffness
khard=5e-3;     
%}

%% Non-Linear (Case 2)
%{
% BoucWen parameters
asiv=0.9;          
xy=0.2;      
nsiv=10;          
n1siv=0.5;       
n2siv=1-n1siv;   

% Degradation parameters (mostagel models)
stiffdeg=7e-10;    
strdegHmos=7e-10;

% gap-closing parameters
kgap=4*0;         
Ngap=1.9;      
xgap=2;   

% Cubic spring stiffness
khard=5e-3;     
%}

%% Non-Linear (Case 3)
%{
% BoucWen parameters
asiv=0.9;         
xy=0.2;       
nsiv=10;         
n1siv=0.5;       
n2siv=1-n1siv;   

% Degradation parameters (mostagel models)
stiffdeg=5e-1;    
strdegHmos=5e-1;  

% gap-closing parameters
kgap=4*0;         
Ngap=1.9;       
xgap=2;   

% Cubic spring stiffness
khard=5e-3;     
%}

%% Non-Linear (Case 4)
%{
% BoucWen parameters
asiv=0.5;          
xy=0.2;       
nsiv=1.5;          
n1siv=0.25;       
n2siv=1-n1siv;   

% Degradation parameters (mostagel models)
stiffdeg=2e-1;    
strdegHmos=2e-1;  

% gap-closing parameters
kgap=4*0;         
Ngap=1.9;       
xgap=0.3;   

% Cubic spring stiffness
khard=5e-3;     
%}

Fyi=(1-asiv)*ke*xy; % Initial yield force

%SS exp. sub. (lineal section)
ssAe=[0,1;-asiv*ke/me,-ce/me];
ssBe=[0;1/me];
ssCe=[1,0;0,1];
ssDe=[0;0];


%%
%%%%%%%%%%%%%%%%%% Carrion & Spencer Actuator Parameters %%%%%%%%%%%%%%%%%%
% Controller
Kp = 3.0;       % mA/in
% Servovalve
tauv = 0.00332; %s
kv = 1;
Kq = 23.01;     % in^3/s/mA
Kc = 1.36e-5;   %in^3/s/psi
%Actuator

Area = 0.751;   % in^2
Cl = 5.89e-6;   % in^3/s/psi
Vt = 48.66;     % in^3
Be = 95958;     % psi


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPENSATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivatives
u = [1 0 0 0];
dudt = [1 -1 0 0]/dtsim;
dudt2 = [1 -2 1 0]/dtsim^2;
dudt3 = [1 -3 3 -1]/dtsim^3;
Derivadas = [u;dudt;dudt2;dudt3];

% Transfer system and model for initial conditions
A_amb_i = [1 20/1000 8e-5 2e-7];
G_act = tf(1,flip(A_amb_i));

% Feedforward controller
% Filter
n = 4;
fc = 20;
fs = 1/dtsim;
[numfilter,denfilter] = butter(n,fc/(fs/2));

% Initial controller
dAMB_i = A_amb_i*Derivadas;
fir_coef = length(dAMB_i);

% RLS design parameters
forgfact = 1;
P_corr_i = eye(4)*1e10;


%% 
%%%%%%%%%%%%%%%%%%%%%%%%% FORCE & DISP. SENSOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%
seed = rng(123);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Sensor RMS noise (Displacement transducer)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rms_noise_DT = 6.25e-13;        % 
Seed_snr = rng(seed);           % Random number generator for noise
Gsnsr_DT = 1;                   % Sensor Gain
rmsD=sqrt(rms_noise_DT/(dtsim));   %rms asociado en N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Sensor RMS noise (Force transducer)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rms_noise_FT = 1.16e-13;            % 
Seed_snr = rng(seed);           % Random number generator for noise
Gsnsr_FT = 1;                   % Sensor Gain
rmsF=sqrt(rms_noise_FT/(dtsim));   %rms asociado en N

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DDMRT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WIENER FILTER
% dt = 0.02
w_80 = [4.15196215127984;-8.61304476276613;11.5668233445886;-11.3369579306043;8.58283182617035;-4.94233392101430;1.94188198247825;-0.394838977645662];

% Interpolation (Monomial Basis)
N1=80;
N=4;
x0=linspace(0,1,N)';
V=x0.^(0:N-1);
V1=inv(V);   %faster than \
x02=linspace(0,1,(N-1)*N1+1);
x02=x02(end-N1+1:end);

















