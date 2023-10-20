%Array processing course basic code
clear
clc
close all
format shortG
%+++++ BEAMFORMING ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----- Scenario -----
%Number of elements in the array
N = 10;
k = 3:2:N+90;
Ns = 100;
% SINAR_partial = zeros(length(Rk),1);
% AWN_partial = zeros(length(Rk),1);
let_through_level = 0e-10;
[mvdr_gsc,mpdr_gsc,SINR_MVDR_SMI,AWN_partial] = deal(zeros(Ns,length(k)));
[SINAR_partial,A_WN_MVDR_SMI] = deal(zeros(Ns,length(k),4));
%Inter-element spacing (in wavelength)
d = 0.5;
pos = d * (0:N-1)'; %positions of the antennas
%Mainlobe width
theta_3dB = 0.9/(N*d);
errors = [-10 -5 0 5 10]*theta_3dB;
%White noise
sigma2 = 1;	%white noise power
%Interference
NoI = 5;
thetaj = [linspace(-40,-20,floor(NoI/2))';linspace(20,40,ceil(NoI/2))']/180*pi;	%angles of arrival	
INR = 20*ones(NoI,1);			%interference to noise ratio (dB)
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
theta0 = 0/180*pi;
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


i = 1;
while (i <= length(k))
error = 1;
while (error <= length(errors))
sample = 1;
% while (sample <= Ns)
%----- ADAPTIVE BEAMFORMING WITH ESTIMATED COVARIANCE MATRICES -----
%Number of snapshots
K = k(i);
while (sample <= Ns)
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
SINR_MVDR_SMI(sample,i) = Ps*(abs(w_MVDR_SMI'*as)^2)/(abs(w_MVDR_SMI'*C*w_MVDR_SMI));
A_WN_MVDR_SMI(sample,i) = 1 / (norm(w_MVDR_SMI)^2);

%----- GSC implementation of MVDR beamformer -----
%Matrix B
B = null(a0')+let_through_level*a0;
%MVDR
Y = Y_MVDR;
%Data in the main and auxilliary channels
d = w_CBF' * Y;     %signal in main channel 1|K
Z = B' * Y;         %signal in auxilliary channels N-1|K
Rz = (Z*Z')/K;      %estimate of R_z
% R = kr;
% [V, D] = eig(Rz); % sorted eigen value the vector with the highest eigen values are the one containing information on interferences
% U = V(:,end-R+1:end);
err = errors(error);
thetaj_err = thetaj+err;
Aj_err = exp(1i*2*pi*pos*sin(thetaj_err'));
U = B'*Aj_err;
Z_t = U'*Z;
Rz_t = (Z_t*Z_t')/K;
rdz = Z*d'/K;
Wa_t = (U'*Rz*U)\U'*rdz;

%\U'*B'*Rz*W_CBF;
w_MVDR_SMI_GSC = w_CBF-B*U*Wa_t;
%Estimate covariance matrix and cross correlation
G_MVDR_SMI = 20*log10(abs(w_MVDR_SMI_GSC'*A));
SINAR_partial(sample,i,error) = Ps*(abs(w_MVDR_SMI_GSC'*as)^2)/(abs(w_MVDR_SMI_GSC'*C*w_MVDR_SMI_GSC));
AWN_partial(sample,i,error) = 1 / (norm(w_MVDR_SMI_GSC)^2);
% legend('boxoff')
sample = sample + 1;
end
error = error + 1;
end
i = i + 1;
end

figure
plot(k,10*log10(mean(SINAR_partial(:,:,3))),'k-x')
hold on
plot(k,10*log10(mean(SINAR_partial(:,:,2))),'k--')
hold on
plot(k,10*log10(mean(SINAR_partial(:,:,4))),'k-o')
hold on
plot(k,10*log10(mean(SINAR_partial(:,:,1))),'k-.')
hold on
plot(k,10*log10(mean(SINAR_partial(:,:,5))),'k-^')
hold on
plot(k,10*log10(SINR_opt)*ones(length(k),1),'k-','LineWidth',0.7)
xlabel('Number of snapshots')
ylabel('SINR (dB)')
legend('Optimal no error','Optimal -0.9°','Optimal +0.9°','Optimal -1.8°','Optimal +1.8°','w_{opt}','Location','southeast')
grid on

