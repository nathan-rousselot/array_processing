%Array processing course basic code
clear
clc
close all
format shortG
rng(42)
%+++++ BEAMFORMING ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----- Scenario -----
%Number of elements in the array
N = 10;
Ns = 1; % Monte Carlo Samples
roots_mpdr = zeros(Ns,N-1);

% k = logspace(1,2,30);
k = 1000;
[opt_arr,cbf_arr,mvdr_arr,mpdr_arr] = deal(zeros(Ns,length(k)));
[opt_arr1,cbf_arr1,mvdr_arr1,mpdr_arr1] = deal(zeros(Ns,length(k)));

i = 1;
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

while i <= length(k)
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
theta0 = 0/180*pi;
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
K = round(k(i));
sample = 1;
while sample <= Ns
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
mvdr_arr(sample,i) = SINR_MVDR_SMI;
A_WN_MVDR_SMI = 1 / (norm(w_MVDR_SMI)^2);
%MPDR-SMI
Y_MPDR = S + IN + NOISE;
R_hat = (Y_MPDR*Y_MPDR')/K;
w_MPDR_SMI = (R_hat\a0);
w_MPDR_SMI = w_MPDR_SMI / (a0'*w_MPDR_SMI);

%%% UNIT CIRCLE RECTIFIED MPDR
xi_n = roots(w_MPDR_SMI);
roots_mpdr(sample,:) = xi_n;

% Initialize variables
arg_xi_n = angle(xi_n);
omega_n = zeros(size(xi_n));

for i = 1:length(xi_n)
    if abs(arg_xi_n(i) - pi*theta0) < 2*pi/N
        omega_n(i) = sign(arg_xi_n(i))*2*pi/N;
    else
        omega_n(i) = angle(exp(1j*arg_xi_n(i)));
    end
    xi_n(i) = exp(1j*omega_n(i));
end
c_n = poly(xi_n).';
w_MPDR_SMI_UC = c_n/abs(c_n'*a0);


G_MPDR_SMI = 20*log10(abs(w_MPDR_SMI'*A));
SINR_MPDR_SMI = Ps*(abs(w_MPDR_SMI'*as)^2)/(abs(w_MPDR_SMI'*C*w_MPDR_SMI));
mpdr_arr(sample,i) = SINR_MPDR_SMI;
A_WN_MPDR_SMI = 1 / (norm(w_MPDR_SMI)^2);
A_WN_CBF = 1 / (norm(w_CBF)^2);
A_WN_opt = 1 / (norm(w_opt)^2);

mvdr_arr1(sample,i) = A_WN_MVDR_SMI;
opt_arr1(sample,i) = A_WN_opt;
cbf_arr1(sample,i) = A_WN_CBF;
mpdr_arr1(sample,i) = A_WN_MPDR_SMI;
sample=sample+1;
end
i = i + 1;
end

figure
plot(roots_mpdr,'d','LineWidth',1)
hold on
plot(xi_n,'ro','LineWidth',1.5)
hold on
fplot(@(t) sin(t), @(t) cos(t),'k');
grid on
axis equal
xlim([-1.5 1.5])
ylim([-1.5 1.5])
legend('Original Zeroes','Unit Circle Rectified')
title('Unit Circle and zeros of array polynomial')



%+++++ BEAMFORMING ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----- Scenario -----
%Number of elements in the array
N = 10;
Ns = 1000; % Monte Carlo Samples
roots_mpdr = zeros(Ns,N-1);

% k = logspace(1,2,30);
k = 1000;
[opt_arr,cbf_arr,mvdr_arr,mpdr_arr] = deal(zeros(Ns,length(k)));
[opt_arr1,cbf_arr1,mvdr_arr1,mpdr_arr1] = deal(zeros(Ns,length(k)));

i = 1;
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

while i <= length(k)
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
theta0 = 0/180*pi;
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
K = round(k(i));
sample = 1;
while sample <= Ns
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
mvdr_arr(sample,i) = SINR_MVDR_SMI;
A_WN_MVDR_SMI = 1 / (norm(w_MVDR_SMI)^2);
%MPDR-SMI
Y_MPDR = S + IN + NOISE;
R_hat = (Y_MPDR*Y_MPDR')/K;
w_MPDR_SMI = (R_hat\a0);
w_MPDR_SMI = w_MPDR_SMI / (a0'*w_MPDR_SMI);

%%% UNIT CIRCLE RECTIFIED MPDR
xi_n = roots(w_MPDR_SMI);
roots_mpdr(sample,:) = xi_n;



G_MPDR_SMI = 20*log10(abs(w_MPDR_SMI'*A));
SINR_MPDR_SMI = Ps*(abs(w_MPDR_SMI'*as)^2)/(abs(w_MPDR_SMI'*C*w_MPDR_SMI));
mpdr_arr(sample,i) = SINR_MPDR_SMI;
A_WN_MPDR_SMI = 1 / (norm(w_MPDR_SMI)^2);
A_WN_CBF = 1 / (norm(w_CBF)^2);
A_WN_opt = 1 / (norm(w_opt)^2);

mvdr_arr1(sample,i) = A_WN_MVDR_SMI;
opt_arr1(sample,i) = A_WN_opt;
cbf_arr1(sample,i) = A_WN_CBF;
mpdr_arr1(sample,i) = A_WN_MPDR_SMI;
sample=sample+1;
end
i = i + 1;
end

figure
plot(roots_mpdr,'r.','LineWidth',1)
hold on
fplot(@(t) sin(t), @(t) cos(t),'k');
grid on
axis equal
xlim([-1.5 1.5])
ylim([-1.5 1.5])
legend('Sampled Zeroes')
title('Unit Circle and zeros of array polynomial')