%Array processing course basic code
clear
clc
close all
format shortG
%+++++ BEAMFORMING ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----- Scenario -----
%Number of elements in the array
N = 10;
k = 10:1:100;
let_through_level = 0e-10;
Ns = 10;
[mvdr_gsc,mpdr_gsc] = deal(zeros(Ns,length(k)));
%Inter-element spacing (in wavelength)
d = 0.5;
pos = d * (0:N-1)'; %positions of the antennas
%Mainlobe width
theta_3dB = 0.9/(N*d);
%White noise
sigma2 = 1;	%white noise power
%Interference
thetaj = [-20;15;35]/180*pi;	%angles of arrival	
INR = [20;20;20];			%interference to noise ratio (dB)
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
%----- ADAPTIVE BEAMFORMING WITH ESTIMATED COVARIANCE MATRICES -----
%Number of snapshots
K = k(i);
sample = 1;
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
SINR_MVDR_SMI = Ps*(abs(w_MVDR_SMI'*as)^2)/(abs(w_MVDR_SMI'*C*w_MVDR_SMI));
A_WN_MVDR_SMI = 1 / (norm(w_MVDR_SMI)^2);

%----- GSC implementation of MVDR beamformer -----
%Matrix B
B = null(a0')+let_through_level*a0;
%MVDR
Y = Y_MVDR;
%Data in the main and auxilliary channels
d = w_CBF' * Y;     %signal in main channel 1|K
Z = B' * Y;         %signal in auxilliary channels N-1|K
Rz = (Z*Z')/K;      %estimate of R_z
R = 2;
[V, D] = eig(Rz); % sorted eigen value the vector with the highest eigen values are the one containing information on interferences
U = V(:,end-R+1:end);
Z_t = Z' * U;
Wa_t = zeros(R,1);
for m=R:-1:0
    disp(m);
    Wa_t(m+1) = D(end-m,end-m);
end
%\U'*B'*Rz*W_CBF;
w_MVDR_SMI_GSC = B*U*Wa_t;
%Estimate covariance matrix and cross correlation

% disp(['VARIANTE MVDR: ||w_{gsc}-w_{df}||=',num2str(norm(w_MVDR_SMI-w_MVDR_SMI_GSC))])
mvdr_gsc(sample,i) = norm(w_MVDR_SMI-w_MVDR_SMI_GSC,2);

sample = sample + 1;
end
i = i + 1;
end

figure;
semilogy(k,mvdr_gsc,'k.')
hold on

xlabel('K')
ylabel('||w_{gsc}-w_{df}||_2')
grid on
% legend('boxoff')

