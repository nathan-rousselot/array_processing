%Array processing course basic code
clear
clc
close all
format shortG
rng(42)
figure;
%+++++ BEAMFORMING ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----- Scenario -----
%Number of elements in the array
N = 10;
Ns = 100; % Monte Carlo Samples
mu_arr = logspace(0,2.5,20);
SINR_mu_arr = zeros(size(mu_arr));
%%%% Hyperparams
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

%----- Natural beampattern -----
% %Weight vector (all weights equal to 1/N)
w = ones(N,1); 
w = w/(ones(1,N)*w);  
%Diagram
tab_theta = (-90:0.5:90)/180*pi;        %Angles where to evaluate beampattern
A = exp(1i*2*pi*pos*sin(tab_theta));    %Steering matrix: each column is a(theta)
G = 20*log10(abs(w'*A));        %beampattern (power in dB)
w_opt = (C\as); 
w_opt = w_opt/(as'*w_opt);
G_opt = 20*log10(abs(w_opt'*A));
SINR_opt = Ps*(abs(w_opt'*as)^2)/(abs(w_opt'*C*w_opt));

j = 1;
while (j <= length(mu_arr))
mu = mu_arr(j);
i = 1;
shift_theta = linspace(-30,30,31);
[opt_arr,cbf_arr,mvdr_arr,mpdr_arr] = deal(zeros(1,length(shift_theta)));
[opt_arr1,cbf_arr1,mvdr_arr1,mpdr_arr1] = deal(zeros(1,length(shift_theta)));

while i <= length(shift_theta)

%----- CONVENTIONAL AND OPTIMAL BEAMFORMERS -----
%Looked direction
theta0 = (shift_theta(i)*theta_3dB)/180*pi;
a0 = exp(1i*2*pi*pos*sin(theta0));

%Optimal beamformer

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
w_MPDR_SMI = ((R_hat+mu*eye(length(R_hat)))\a0);
w_MPDR_SMI = w_MPDR_SMI / (a0'*w_MPDR_SMI);
G_MPDR_SMI = 20*log10(abs(w_MPDR_SMI'*A));
SINR_MPDR_SMI = Ps*(abs(w_MPDR_SMI'*as)^2)/(abs(w_MPDR_SMI'*C*w_MPDR_SMI));
mpdr_arr(i) = SINR_MPDR_SMI;
A_WN_MPDR_SMI = 1 / (norm(w_MPDR_SMI)^2);
A_WN_opt = 1 / (norm(w_opt)^2);
opt_arr1(i) = A_WN_opt;
mvdr_arr1(i) = A_WN_MVDR_SMI;
mpdr_arr1(i) = A_WN_MPDR_SMI;


i = i + 1;
end

SINR_mu_arr(j) = norm(SINR_MPDR_SMI-SINR_opt,2)/norm(SINR_opt,2);


%%%%%%%%%%%%%%%%%%%% K Convergence %%%%%%%%%%%%%%%%%%%%

k = logspace(1,2,30);
[opt_arr2,cbf_arr2,mvdr_arr2,mpdr_arr2] = deal(zeros(Ns,length(k)));
[opt_arr3,cbf_arr3,mvdr_arr3,mpdr_arr3] = deal(zeros(Ns,length(k)));
i = 1;
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
mvdr_arr2(sample,i) = SINR_MVDR_SMI;
A_WN_MVDR_SMI = 1 / (norm(w_MVDR_SMI)^2);
%MPDR-SMI
Y_MPDR = S + IN + NOISE;
R_hat = (Y_MPDR*Y_MPDR')/K;
w_MPDR_SMI = ((R_hat+mu*eye(length(R_hat)))\a0);
w_MPDR_SMI = w_MPDR_SMI / (a0'*w_MPDR_SMI);
G_MPDR_SMI = 20*log10(abs(w_MPDR_SMI'*A));
SINR_MPDR_SMI = Ps*(abs(w_MPDR_SMI'*as)^2)/(abs(w_MPDR_SMI'*C*w_MPDR_SMI));
mpdr_arr2(sample,i) = SINR_MPDR_SMI;
A_WN_MPDR_SMI = 1 / (norm(w_MPDR_SMI)^2);
A_WN_opt = 1 / (norm(w_opt)^2);

mvdr_arr3(sample,i) = A_WN_MVDR_SMI;
opt_arr3(sample,i) = A_WN_opt;
mpdr_arr3(sample,i) = A_WN_MPDR_SMI;
sample = sample + 1; % Monte Carlo loop
end
i = i + 1;
end

clf;

sgtitle("Constrained Method with \mu = " + mu)
subplot(221)
plot(shift_theta*theta_3dB,10*log10(opt_arr),'k-^','LineWidth',0.7)
hold on
plot(shift_theta*theta_3dB,10*log10(mvdr_arr),'k--o','LineWidth',0.7)
hold on
plot(shift_theta*theta_3dB,10*log10(mpdr_arr),'k--x','LineWidth',0.7)
legend('Optimal','MVDR','MPDR')
xlabel('\Delta\theta (°)')
ylabel('SINR (dB)')
grid on
title('SINR s.t. \Delta\theta')

subplot(222)
plot(shift_theta*theta_3dB,10*log10(opt_arr1),'k-^','LineWidth',0.7)
hold on
plot(shift_theta*theta_3dB,10*log10(mvdr_arr1),'k--o','LineWidth',0.7)
hold on
plot(shift_theta*theta_3dB,10*log10(mpdr_arr1),'k--x','LineWidth',0.7)
legend('Optimal','MVDR','MPDR')
xlabel('\Delta\theta (°)')
ylabel('SINR (dB)')
grid on
title('A_{WN} s.t. \Delta\theta')

subplot(223)
plot(k,10*log10(opt_arr2),'k-^','LineWidth',0.7)
hold on
plot(k,10*log10(mean(mvdr_arr2)),'k--o','LineWidth',0.7)
hold on
plot(k,10*log10(mean(mpdr_arr2)),'k--x','LineWidth',0.7)
legend('Optimal','MVDR','MPDR')
xlabel('Number of snapshots')
ylabel('SINR (dB)')
grid on
title('SINR Convergence (s.t. k)')

subplot(224)
plot(round(k),10*log10(mean(opt_arr3)),'k-^','LineWidth',0.7)
hold on
plot(round(k),10*log10(mean(mvdr_arr3)),'k--o','LineWidth',0.7)
hold on
plot(round(k),10*log10(mean(mpdr_arr3)),'k--x','LineWidth',0.7)
legend('Optimal','MVDR','MPDR')
xlabel('Number of snapshots')
ylabel('A_WN (dB)')
grid on
title('A_{WN} convergence (s.t. k)')
drawnow;
% pause(0.01);
j = j+1;
end
figure
plot(mu_arr,SINR_mu_arr,'k','LineWidth',1)
grid on
xlabel('\mu')
ylabel("Relative Error")
title('Robust MPDR relative error to Optimal in function of \mu')
x = [0 max(mu_arr)];
[mu_sinr_opt, minx] = min(SINR_mu_arr);
mu_opt = mu_arr(minx);
y = [mu_sinr_opt,mu_sinr_opt];
line(x,y,'Color','black','LineStyle','--');
legend('Relative Error',"Minimum reached for \mu = " + mu_opt);



%%% Optimal result
mu = mu_opt;
K = 100;
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
% mvdr_arr2(sample,i) = SINR_MVDR_SMI;
A_WN_MVDR_SMI = 1 / (norm(w_MVDR_SMI)^2);
%MPDR-SMI
Y_MPDR = S + IN + NOISE;
R_hat = (Y_MPDR*Y_MPDR')/K;
w_MPDR_SMI = ((R_hat+mu*eye(length(R_hat)))\a0);
w_MPDR_SMI = w_MPDR_SMI / (a0'*w_MPDR_SMI);
G_MPDR_SMI = 20*log10(abs(w_MPDR_SMI'*A));
SINR_MPDR_SMI = Ps*(abs(w_MPDR_SMI'*as)^2)/(abs(w_MPDR_SMI'*C*w_MPDR_SMI));
% mpdr_arr2(sample,i) = SINR_MPDR_SMI;
A_WN_MPDR_SMI = 1 / (norm(w_MPDR_SMI)^2);
A_WN_opt = 1 / (norm(w_opt)^2);

figure
plot(tab_theta*180/pi,G_opt,'-',...
    tab_theta*180/pi,G_MVDR_SMI,'-.',...
    tab_theta*180/pi,G_MPDR_SMI,':.','linewidth',2);
for k=1:length(thetaj)
    xline(thetaj(k)*180/pi,'--','color','k','linewidth',2)
end
title(['Beampatterns $K=$',num2str(K)],'interpreter','latex');
ylabel('dB','interpreter','latex');
xlabel('Angle of Arrival (degrees)','interpreter','latex');
legend('opt','MVDR-SMI','MPDR-SMI');
axis([-90 90 -70 10]);
grid on
