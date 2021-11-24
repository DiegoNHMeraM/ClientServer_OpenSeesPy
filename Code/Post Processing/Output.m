%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post-Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Adaptive references %%%%%%%%%%%%%%%%%
% Referencia = load('Ref_Linear.mat');
% Referencia = load('Ref_NonLinear_1.mat');
% Referencia = load('Ref_NonLinear_2.mat');
% Referencia = load('Ref_NonLinear_3.mat');
% Referencia = load('Ref_NonLinear_4.mat');

%%%%%%%%%%%%% Fixed references %%%%%%%%%%%%%%%%%
% Referencia = load('Ref_Linear_Fixed.mat');
% Referencia = load('Ref_NonLinear_1_Fixed.mat');
% Referencia = load('Ref_NonLinear_2_Fixed.mat');
% Referencia = load('Ref_NonLinear_3_Fixed.mat');
% Referencia = load('Ref_NonLinear_4_Fixed.mat');


WIENER = out.x_t.Data;
x_t = zeros(length(WIENER),1);
for i = 1:length(x_t)
    x_t(i) = WIENER(:,:,i);
end
t_t = out.x_t.Time;
x_m = out.x_m.Data;
t_m = out.x_m.Time;
totaltime = t_t(end);

% Hysteresis
% Displacement results
Disp = out.Desp.Data;       % Measured
t = out.Desp.Time;
% Strength results
rtn = out.rt.Data(:,1);     % Non-linear restoring force (without inertial or viscous)
Kcurr = out.Kcur.Data(:,1); % Rigidity over time
Fyi = out.Fy.Data(:,1);     % Fluence over time
F_hard = out.Fhard.Data;    % Cubic Spring Force
F_m = out.Fm.Data;          % Measured experimental force


% Outcome indicators
[Amptotal,phitotal,feqtotal,delta] = Freq_Resp_Tong(Referencia.x_t,x_m,1/(1/4000));
aux = abs(round(delta,2)*4000)+1;

J2=rms(x_t-x_m)/rms(x_t)*100;
J4_Disp=rms(Referencia.x_t(aux:end)-x_m(1:end-(aux-1)))/rms(Referencia.x_t(aux:end))*100;
J4_Force=rms(Referencia.rtn(aux:end,1)-rtn(1:end-(aux-1)))/rms(Referencia.rtn(aux:end,1))*100;
[Amptotal,phitotal,feqtotal,delaytotal] = Freq_Resp_Tong(x_t,x_m,1/(1/4000));


Ji = ["J2 [%]";"J4 Disp [%]";"J4 Force [%]";"delay [ms]"];
Ji = cellstr(Ji);
results=[J2;J4_Disp;J4_Force;delaytotal*1000];
table(results,'VariableNames',{'Results'},'RowNames',Ji)

%%%%%%%%%%%%%%%%%%%%%%%%%%% TRACKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,4,[1,2,3])
plot(t_t,x_t,'k','LineWidth',1.5)
hold on
plot(t_m,x_m,'r--','LineWidth',1.5)
legend('x_{Target}','x_{Measured}','Orientation','Horizontal','Location','best','FontSize',12)
xlabel('Time [s]','FontSize',15)
ylabel('Disp. [in]','FontSize',15)
xlim([0 totaltime])
grid on

subplot(2,4,[5,6,7])
plot(t_t,abs(x_t-x_m),'k','LineWidth',1.5)
hold on
plot(t_t(1,1),abs(x_t(1,1)-x_m(1,1)),'k','LineWidth',1.5)
legend(['J_2 = ',num2str(J2),' %'],['Delay = ',num2str(abs(delaytotal)*1000), ' ms'],'FontSize',12)  
xlabel('Time [sec]','FontSize',15)
ylabel('|error| [in]','FontSize',15)
xlim([0 totaltime])
grid on

subplot(2,4,[4,8])
plot(t_t,x_t,'k','LineWidth',1.5)
hold on
plot(t_m,x_m,'r--','LineWidth',1.5)
legend('x_{Target}','x_{Measured}','Orientation','Vertical','Location','best','FontSize',12)
xlabel('Time [s]','FontSize',15)
ylabel('Disp. [in]','FontSize',15)
xlim([11.5 13])
grid on



%%%%%%%%%%%%%%%%%%%%%%%%%%%% REFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,4,[1,2,3])
plot(Referencia.t_t+delta,Referencia.x_t,'k','LineWidth',1.5)
hold on
plot(t_m,x_m,'r--','LineWidth',1.5)
legend('x_{Reference}','x_{Measured}','Orientation','Horizontal','Location','best','FontSize',12)
xlabel('Time [s]','FontSize',15)
ylabel('Disp. [in]','FontSize',15)
xlim([0 totaltime])
grid on

subplot(2,4,[5,6,7])
plot(t_t(1:end-(aux-1)),abs(Referencia.x_t(aux:end)-x_m(1:end-(aux-1))),'k','LineWidth',1.5)
legend(['J_4 = ',num2str(J4_Disp),' %'],'FontSize',12)  
xlabel('Time [sec]','FontSize',15)
ylabel('|error| [in]','FontSize',15)
xlim([0 totaltime])
grid on

subplot(2,4,[4,8])
plot(Referencia.t_t+delta,Referencia.x_t,'k','LineWidth',1.5)
hold on
plot(t_m,x_m,'r--','LineWidth',1.5)
legend('x_{Reference}','x_{Measured}','Orientation','Vertical','Location','best','FontSize',12)
xlabel('Time [s]','FontSize',15)
ylabel('Disp. [in]','FontSize',15)
xlim([11.5 13])
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%% W PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_adap = out.adaptive_parameters.Data;
w = zeros(length(w_adap),8);
for i = 1:length(w)
    w(i,:) = w_adap(:,:,i);
end


figure; hold on; grid on
plot(out.adaptive_parameters.Time,w,'LineWidth',1.5)
legend('w_1','w_2','w_3','w_4','w_5','w_6','w_7','w_8','Location','EastOutside','FontSize',12)
xlabel('Time [s]','FontSize',15); ylabel('Coefficients','FontSize',15);
title('FIR Filter Coefficients')
xlim([0 out.adaptive_parameters.Time(end)])


figure()
 for plotId = 1 : 8
     subplot(2, 4, plotId) ;
     
     x=out.adaptive_parameters.Time;
     y = w(:,plotId);
     [p, S, mu] = polyfit(x,y,10);

     % Plots
     hold on
     plot(x,y,'ro','markersize',2,'markerfacecolor','white')
     z=@(x) polyval(p,x,S,mu);
     fplot(z,[x(1),x(end)],'k','Linewidth',1.5)
     legend('Data',['w',num2str(plotId)],'Location','best','FontSize',12)
     xlabel('Time [s]','FontSize',15); ylabel('Coefficient','FontSize',15);
     xlim([0 out.adaptive_parameters.Time(end)])
     grid on
 end
 sgtitle('FIR Filter Coefficients')

%%%%%%%%%%%%%%%%%%%%%%%%%%%% HYSTERESIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on
plot(Disp,rtn,'r-','Linewidth',1.0)
plot(Referencia.x, Referencia.rtn,'k--','Linewidth',1.0)
legend('Hysteresis','Reference hysteresis','Location','best','FontSize',10)
xlabel('Measured Disp [in]','FontSize',12)
ylabel('Force [Kips]','FontSize',12)
grid on

%%% Only in Local Simulations 
%{
time = out.MT.Time;
ticks = out.MT.Data;

figure
plot(time,ticks,'k','LineWidth',1.5)
hold on
plot(time(find(ticks==max(abs(ticks)))),max(ticks),'*r','LineWidth',1.5)
xlabel('Time [sec]','FontSize',15)
ylabel('Missed Ticks','FontSize',15)
legend('Missed Ticks',['Max_{MT} = ',num2str(max(ticks))],'FontSize',12)
grid on
%}


function [A,phi,feq,d] = Freq_Resp_Tong(InSig,OutSig,fsamp)
%        A : Amplitude error (Output / Input)
%       phi: Phase error
%       feq: Equivalent Frequency
%        d : Delay
%     InSig: Input Signal
%    OutSig: Output Signal
%     fsamp: Sampling frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l=2;    % Exponent ||fft[I(t)]||^l

% Signal FFT
FI0=fft(InSig);         
FO0=fft(OutSig); 

% Input Signal Process
L=length(InSig);                              % Signal size
FI=FI0(1:ceil(L/2));                          % One Sided Spectrum
FI(2:end-1)=2*FI(2:end-1);
fi=fsamp*(0:ceil(L/2)-1)/L;                   % Frequencies

feq=sum((abs(FI).^l).*fi')/sum(abs(FI).^l);   % Equivalent Frequency

% Output Signal Process
L=length(InSig);                              % Signal size
FO=FO0(1:ceil(L/2));                          % One Sided Spectrum
FO(2:end-1)=2*FO(2:end-1);
fo=fsamp*(0:ceil(L/2)-1)/L;                   % Frequencies      


%FEI
sumFI=sum(abs(FI).^l);          
FEI=FO.*abs(FI).^l./FI/sumFI;   
FEI=sum(FEI);

% Amplitude and Phase
A=abs(FEI);
phi=atan(imag(FEI)/real(FEI));

% Time delay
d=-phi/(2*pi*feq); 

end