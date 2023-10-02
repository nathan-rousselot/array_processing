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


%----- Natural beampattern -----
% %Weight vector (all weights equal to 1/N)
w = ones(N,1); 
w = w/(ones(1,N)*w);  
%Diagram
tab_theta = (-90:0.5:90)/180*pi;        %Angles where to evaluate beampattern
A = exp(1i*2*pi*pos*sin(tab_theta));    %Steering matrix: each column is a(theta)
G = 20*log10(abs(w'*A));        %beampattern (power in dB)
%Plot
figure
plot(tab_theta*180/pi,G,'linewidth',2);
title('Natural beampattern $w_n=1/N$','interpreter','latex');
ylabel('dB','interpreter','latex');
xlabel('Angle of Arrival (degrees)','interpreter','latex');
axis([-90 90 -70 10]);
grid on


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
%Plot
figure
plot(tab_theta*180/pi,G_opt,'-',...
    tab_theta*180/pi,G_CBF,'--','linewidth',2);
for k=1:length(thetaj)
    xline(thetaj(k)*180/pi,'--','color','k','linewidth',2)
end
title('Beampatterns of CBF and optimal beamformers','interpreter','latex');
ylabel('dB','interpreter','latex');
xlabel('Angle of Arrival (degrees)','interpreter','latex');
legend('opt','CBF');
axis([-90 90 -70 10]);
grid on



%----- ADAPTIVE BEAMFORMING WITH TRUE COVARIANCE MATRICES -----
%MVDR
w_MVDR = (C\a0); 
w_MVDR = w_MVDR/(a0'*w_MVDR);
G_MVDR = 20*log10(abs(w_MVDR'*A));
SINR_MVDR = Ps*(abs(w_MVDR'*as)^2)/(abs(w_MVDR'*C*w_MVDR));
A_WN_MVDR = 1/(norm(w_MVDR)^2);
%MPDR
w_MPDR = (R\a0); 
w_MPDR = w_MPDR/(a0'*w_MPDR);
G_MPDR = 20*log10(abs(w_MPDR'*A));
SINR_MPDR = Ps*(abs(w_MPDR'*as)^2)/(abs(w_MPDR'*C*w_MPDR));
A_WN_MPDR = 1/(norm(w_MPDR)^2);
%Plot
figure
plot(tab_theta*180/pi,G_opt,'-',...
    tab_theta*180/pi,G_CBF,'--',...
    tab_theta*180/pi,G_MVDR,'-.',...
    tab_theta*180/pi,G_MPDR,':.','linewidth',2);
for k=1:length(thetaj)
    xline(thetaj(k)*180/pi,'--','color','k','linewidth',2)
end
title('Beampatterns $K=\infty$','interpreter','latex');
ylabel('dB','interpreter','latex');
xlabel('Angle of Arrival (degrees)','interpreter','latex');
legend('opt','CBF','MVDR','MPDR');
axis([-90 90 -70 10]);
grid on


disp(['SINR_CBF = ',num2str(10*log10(SINR_CBF)),' dB'])
disp(['SINR_opt = ',num2str(10*log10(SINR_opt)),' dB'])
disp(['SINR_MVDR = ',num2str(10*log10(SINR_MVDR)),' dB'])
disp(['SINR_MPDR = ',num2str(10*log10(SINR_MPDR)),' dB'])

%----- ADAPTIVE BEAMFORMING WITH ESTIMATED COVARIANCE MATRICES -----
%Number of snapshots
K = 50;
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
%MPDR-SMI
Y_MPDR = S + IN + NOISE;
R_hat = (Y_MPDR*Y_MPDR')/K;
w_MPDR_SMI = (R_hat\a0);
w_MPDR_SMI = w_MPDR_SMI / (a0'*w_MPDR_SMI);
G_MPDR_SMI = 20*log10(abs(w_MPDR_SMI'*A));
SINR_MPDR_SMI = Ps*(abs(w_MPDR_SMI'*as)^2)/(abs(w_MPDR_SMI'*C*w_MPDR_SMI));
A_WN_MPDR_SMI = 1 / (norm(w_MPDR_SMI)^2);
A_WN_CBF = 1 / (norm(w_CBF)^2);
A_WN_opt = 1 / (norm(w_opt)^2);
%Plot
figure
plot(tab_theta*180/pi,G_opt,'-',...
    tab_theta*180/pi,G_CBF,'--',...
    tab_theta*180/pi,G_MVDR_SMI,'-.',...
    tab_theta*180/pi,G_MPDR_SMI,':.','linewidth',2);
for k=1:length(thetaj)
    xline(thetaj(k)*180/pi,'--','color','k','linewidth',2)
end
title(['Beampatterns $K=$',num2str(K)],'interpreter','latex');
ylabel('dB','interpreter','latex');
xlabel('Angle of Arrival (degrees)','interpreter','latex');
legend('opt','CBF','MVDR-SMI','MPDR-SMI');
axis([-90 90 -70 10]);
grid on


disp(['SINR_MVDR_SMI = ',num2str(10*log10(SINR_MVDR_SMI)),' dB'])
disp(['SINR_MPDR_SMI = ',num2str(10*log10(SINR_MPDR_SMI)),' dB'])
%White noise array gain
disp(['A_WN_CBF = ',num2str(10*log10(A_WN_CBF)),' dB'])
disp(['A_WN_opt = ',num2str(10*log10(A_WN_opt)),' dB'])
disp(['A_WN_MVDR = ',num2str(10*log10(A_WN_MVDR)),' dB'])
disp(['A_WN_MPDR = ',num2str(10*log10(A_WN_MPDR)),' dB'])
disp(['A_WN_MVDR_SMI = ',num2str(10*log10(A_WN_MVDR_SMI)),' dB'])
disp(['A_WN_MPDR_SMI = ',num2str(10*log10(A_WN_MPDR_SMI)),' dB'])


%----- GSC implementation of MVDR beamformer -----
%Matrix B
B = null(a0');
%MVDR
Y = Y_MVDR;
%Data in the main and auxilliary channels
d = w_CBF' * Y;     %signal in main channel 1|K
Z = B' * Y;         %signal in auxilliary channels N-1|K
%Estimate covariance matrix and cross correlation
Rz = (Z*Z')/K;      %estimate of R_z
rdz = Z*d'/K;       %estimate of R_{dz}
wa = Rz\rdz;        %estimate of R_z^{-1} r_{dz}
w_MVDR_SMI_GSC = w_CBF - B * wa;
disp(['||w_{gsc}-w_{df}||=',num2str(norm(w_MVDR_SMI-w_MVDR_SMI_GSC))])

%+++++++++++ DIRECTION OF ARRIVAL ESTIMATION ++++++++++++++++++++++++++++++
%Number of elements in the array
N = 20;
%Inter-element spacing (in wavelength)
d = 0.5;
pos = d * (0:N-1)'; %positions of the antennas
%Signals impinging on the array 
%White noise
sigma2 = 1;	%white noise power
%Signals
thetas = [-20;0;10]/180*pi;	%angles of arrival	
SNR = [10;10;10];			%interference to noise ratio
Ps = sigma2 * 10.^(SNR/10);		%interference power
P = length(thetas);
As = exp(1i*2*pi*pos*sin(thetas'));
%Number of snapshots
K = N;
%Number of FFT bins
vec_doa = (-90:0.2:90)/180*pi;       
A = exp(1i*2*pi*pos*sin(vec_doa)); 
%Snapshots
Y = As * diag(sqrt(Ps/2)) * (randn(P,K)+1i*randn(P,K))...
    + sqrt(sigma2/2)*(randn(N,K)+1i*randn(N,K));

%Conventional beamforming
P_CBF = sum(abs(A'*Y).^2,2)/(N^2)/K;
P_CBF = 10*log10(P_CBF);

figure
vec_doa = vec_doa*180/pi;
plot(vec_doa,P_CBF,'-','Linewidth',2);
title('CBF spectrum','interpreter','latex');
xlabel('Angle of arrival (degrees)','interpreter','latex');
ylabel('dB','interpreter','latex');

