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
NoI = 5;
thetaj = [linspace(-40,-20,floor(NoI/2))';linspace(20,40,ceil(NoI/2))']/180*pi;	%angles of arrival	
INR = 20*ones(NoI,1);			%interference to noise ratio (dB)			%interference to noise ratio (dB)
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
theta0 = 3/180*pi;
a0 = exp(1i*2*pi*pos*sin(theta0));
%Conventional beamformer
w_CBF = a0; 
w_CBF = w_CBF/(a0'*w_CBF);
G_CBF = 20*log10(abs(w_CBF'*A));        %beampattern (power in dB)
SINR_CBF = Ps*(abs(w_CBF'*as)^2)/(abs(w_CBF'*C*w_CBF)); %SINR
A_WN_CBF = 1/(norm(w_CBF)^2);   %White noise array gain
%Optimal beamformer
w_opt = (C\as); 
w_opt = w_opt/(as'*w_opt);
G_opt = 20*log10(abs(w_opt'*A));
SINR_opt = Ps*(abs(w_opt'*as)^2)/(abs(w_opt'*C*w_opt));
A_WN_opt = 1/(norm(w_opt)^2);


%----- ADAPTIVE BEAMFORMING WITH ESTIMATED COVARIANCE MATRICES -----
%Number of snapshots
K = 100;
%Signal
S = sqrt(Ps/2) * as * (randn(1,K)+1i*randn(1,K));
%Interference + noise
IN = Aj * diag(sqrt(Pj/2)) * (randn(J,K)+1i*randn(J,K));
NOISE = sqrt(sigma2/2)*(randn(N,K)+1i*randn(N,K));
%MVDR-SMI
Y_MVDR = IN + NOISE;
C_hat = (Y_MVDR*Y_MVDR')/K;
kr = 1;
[SINR_MVDR_SMI, A_WN_MVDR_SMI] = deal(zeros(1,N));
while kr <= N
    w_MVDR_SMI = conjugate_gradient_method(zeros(length(w_opt),1),C_hat,a0,1e-20,kr);
    w_MVDR_SMI = w_MVDR_SMI / (a0'*w_MVDR_SMI);
    G_MVDR_SMI = 20*log10(abs(w_MVDR_SMI'*A));
    SINR_MVDR_SMI(kr) = Ps*(abs(w_MVDR_SMI'*as)^2)/(abs(w_MVDR_SMI'*C*w_MVDR_SMI));
    A_WN_MVDR_SMI(kr) = 1 / (norm(w_MVDR_SMI)^2);
    kr = kr + 1;
end

w_MVDR_SMI = C_hat\a0;
w_MVDR_SMI = w_MVDR_SMI / (a0'*w_MVDR_SMI);
G_MVDR_SMI = 20*log10(abs(w_MVDR_SMI'*A));
SINR_MVDR_SMI_noob = Ps*(abs(w_MVDR_SMI'*as)^2)/(abs(w_MVDR_SMI'*C*w_MVDR_SMI));
A_WN_MVDR_SMI_noob = 1 / (norm(w_MVDR_SMI)^2);
kr = kr + 1;
%MPDR-SMI
% Y_MPDR = S + IN + NOISE;
% R_hat = (Y_MPDR*Y_MPDR')/K;
% w_MPDR_SMI = (R_hat\a0);
% w_MPDR_SMI = w_MPDR_SMI / (a0'*w_MPDR_SMI);
% G_MPDR_SMI = 20*log10(abs(w_MPDR_SMI'*A));
% SINR_MPDR_SMI = Ps*(abs(w_MPDR_SMI'*as)^2)/(abs(w_MPDR_SMI'*C*w_MPDR_SMI));
% A_WN_MPDR_SMI = 1 / (norm(w_MPDR_SMI)^2);
A_WN_CBF = 1 / (norm(w_CBF)^2);
A_WN_opt = 1 / (norm(w_opt)^2);


figure;
plot(1:1:N,ones(N,1)*10*log10(SINR_opt),'k-','LineWidth',1)
hold on
plot(1:1:N, 10*log10(SINR_MVDR_SMI_noob)*ones(N,1),'--','LineWidth',1);
hold on
plot(1:1:N, 10*log10(SINR_MVDR_SMI),'k--','LineWidth',1)
legend('Optimal','MVDR','CG-MVDR')
xlabel('Krylov Subspace Dimension')
ylabel('SINR (dB)')
grid on
