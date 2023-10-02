%Array processing course basic code
clear
clc
close all
format shortG
%+++++ BEAMFORMING ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----- Scenario -----
%Number of elements in the array
N = 10;
%Inter-element spacing (in wavelength)
d = 0.5;
pos = d * (0:N-1)'; %positions of the antennas
%Mainlobe width
theta_3dB = 0.9/(N*d);
%White noise
sigma2 = 1;	%white noise power
%Interference
thetaj = [-20;15]/180*pi;	%angles of arrival	
INR = [20;20];			%interference to noise ratio (dB)
Pj = sigma2 * 10.^(INR/10);		%interference power
J = length(thetaj);
%Interference + noise covariance matrix
Aj = exp(1i*2*pi*pos*sin(thetaj'));	%interference steering matrix N|J
C = Aj*diag(Pj)*Aj' + sigma2*eye(N);	%interference + noise covariance matrix
%Signal of interest
thetas = 0/180*pi;	%angle of arrival
SNR = 0;            %signal to noise ratio (dB)
Ps = sigma2 * 10^(SNR/10);			%signal power
as = exp(1i*2*pi*pos*sin(thetas));	%steering vector
%Total covariance matrix (signal + interference + noise)
R = Ps*(as*as') + C;

i = 1;
shift_theta = linspace(-30,30,31);
[opt_arr,cbf_arr,mvdr_arr,mpdr_arr] = deal(zeros(size(shift_theta)));
[opt_arr1,cbf_arr1,mvdr_arr1,mpdr_arr1] = deal(zeros(size(shift_theta)));

while i <= length(shift_theta)

%----- Natural beampattern -----
% %Weight vector (all weights equal to 1/N)
w = ones(N,1); 
w = w/(ones(1,N)*w);  
%Diagram
tab_theta = (-90:0.5:90)/180*pi;        %Angles where to evaluate beampattern
A = exp(1i*2*pi*pos*sin(tab_theta));    %Steering matrix: each column is a(theta)
G = 20*log10(abs(w'*A));        %beampattern (power in dB)


%----- CONVENTIONAL AND OPTIMAL BEAMFORMERS -----
%Looked direction
theta0 = (shift_theta(i)*theta_3dB)/180*pi;
a0 = exp(1i*2*pi*pos*sin(theta0));
%Conventional beamformer
w_CBF = a0; 
w_CBF = w_CBF/(a0'*w_CBF);
G_CBF = 20*log10(abs(w_CBF'*A));        %beampattern (power in dB)
SINR_CBF = Ps*(abs(w_CBF'*as)^2)/(abs(w_CBF'*C*w_CBF)); %SINR
cbf_arr(i) = SINR_CBF;
A_WN_CBF = 1/(norm(w_CBF)^2);   %White noise array gain
%Optimal beamformer
w_opt = (C\as); 
w_opt = w_opt/(as'*w_opt);
G_opt = 20*log10(abs(w_opt'*A));
SINR_opt = Ps*(abs(w_opt'*as)^2)/(abs(w_opt'*C*w_opt));
opt_arr(i) = SINR_opt;

A_WN_opt = 1/(norm(w_opt)^2);

%----- ADAPTIVE BEAMFORMING WITH ESTIMATED COVARIANCE MATRICES -----
%Number of snapshots
K = 1000;
%Signal
S = sqrt(Ps/2) * as * (randn(1,K)+1i*randn(1,K));
%Interference + noise
IN = Aj * diag(sqrt(Pj/2)) * (randn(J,K)+1i*randn(J,K));
NOISE = sqrt(sigma2/2)*(randn(N,K)+1i*randn(N,K));
%MVDR-SMI
Y_MVDR = IN + NOISE;
C_hat = (Y_MVDR*Y_MVDR')/K;
w_MVDR_SMI = (C_hat\a0);
w_MVDR_SMI = w_MVDR_SMI / (a0'*w_MVDR_SMI);
G_MVDR_SMI = 20*log10(abs(w_MVDR_SMI'*A));
SINR_MVDR_SMI = Ps*(abs(w_MVDR_SMI'*as)^2)/(abs(w_MVDR_SMI'*C*w_MVDR_SMI));
mvdr_arr(i) = SINR_MVDR_SMI;
A_WN_MVDR_SMI = 1 / (norm(w_MVDR_SMI)^2);
%MPDR-SMI
Y_MPDR = S + IN + NOISE;
R_hat = (Y_MPDR*Y_MPDR')/K;
mu = N;
w_MPDR_SMI = ((R_hat+mu*eye(length(R_hat)))\a0);
w_MPDR_SMI = w_MPDR_SMI / (a0'*w_MPDR_SMI);
G_MPDR_SMI = 20*log10(abs(w_MPDR_SMI'*A));
SINR_MPDR_SMI = Ps*(abs(w_MPDR_SMI'*as)^2)/(abs(w_MPDR_SMI'*C*w_MPDR_SMI));
mpdr_arr(i) = SINR_MPDR_SMI;
A_WN_MPDR_SMI = 1 / (norm(w_MPDR_SMI)^2);
A_WN_CBF = 1 / (norm(w_CBF)^2);
A_WN_opt = 1 / (norm(w_opt)^2);
opt_arr1(i) = A_WN_opt;
cbf_arr1(i) = A_WN_CBF;
mvdr_arr1(i) = A_WN_MVDR_SMI;
mpdr_arr1(i) = A_WN_MPDR_SMI;


i = i + 1;
end


figure
plot(shift_theta*theta_3dB,10*log10(opt_arr),'k-^','LineWidth',0.7)
%hold on
%plot(shift_theta,cbf_arr)
hold on
plot(shift_theta*theta_3dB,10*log10(mvdr_arr),'k--o','LineWidth',0.7)
hold on
plot(shift_theta*theta_3dB,10*log10(mpdr_arr),'k--x','LineWidth',0.7)
legend('Optimal','MVDR','MPDR')
xlabel('\Delta\theta (°)')
ylabel('SINR (dB)')
grid on
title('MVDR vs MPDR vs Optimal : SINR Comparison in function of DOA error')

figure
plot(shift_theta*theta_3dB,10*log10(opt_arr1),'k-^','LineWidth',0.7)
%hold on
%plot(shift_theta,cbf_arr)
hold on
plot(shift_theta*theta_3dB,10*log10(mvdr_arr1),'k--o','LineWidth',0.7)
hold on
plot(shift_theta*theta_3dB,10*log10(mpdr_arr1),'k--x','LineWidth',0.7)
legend('Optimal','MVDR','MPDR')
xlabel('\Delta\theta (°)')
ylabel('SINR (dB)')
grid on
title('MVDR vs MPDR vs Optimal : SINR Comparison in function of DOA error')

