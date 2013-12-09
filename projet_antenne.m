%--------------------------------------------------------------------------
%   Projet Traitement d'antennes
%   Pierre Apap - Thomas Nicot
%   3 EN
%--------------------------------------------------------------------------

%% Application des methodes classiques
clc; clear all; close all;

%SINR versus number of snasphots
%Number of elements in the array
N = 10;
%Inter-element spacing (in wavelength)
d = 0.5;
pos = d * (0:N-1)';

%Signal of interest
thetas = 0/180*pi;	%angle of arrival	
Ps = 1;			%power
as = exp(1i*2*pi*pos*sin(thetas));	%steering vector

%White noise
SNR = 0;		%SNR
sigma2 = Ps*10^(-SNR/10);	%white noise power

%Interference
thetaj = [-20;15]/180*pi;	%angles of arrival	
INR = [20;15];			%interference to noise ratio
Pj = sigma2 * 10.^(INR/10);		%interference power
J = length(thetaj);

%Interference + noise covariance matrix
Aj = exp(1i*2*pi*pos*sin(thetaj'));	%interference steering matrix N|J
Cth = Aj*diag(Pj)*Aj' + sigma2*eye(N);	%covariance matrix

%Total covariance matrix
Rth = Ps*(as*as') + Cth;

%Looked angle
theta0 = 0/180*pi;
a0 = exp(1i*2*pi*pos*sin(theta0));

%Blocking matrix
B = null(a0');

%Optimum SINR
sinr_opt = 10*log10( Ps * abs(as'*(Cth\as)) );

%CBF
w_cbf = a0/(a0'*a0);
sinr_cbf = 10*log10( Ps * abs(w_cbf'*as)^2 / abs(w_cbf'*Cth*w_cbf) );

%Number snapshots
K = 25;

%Number of Monte-Carlo simulations
niter = 50;

%Arrays to store results
sinr_mpdr_df = zeros(niter,1);
sinr_mvdr_df = zeros(niter,1);
sinr_mpdr_gsc = zeros(niter,1);
sinr_mvdr_gsc = zeros(niter,1);
Aw_mpdr = zeros(niter,1);
Aw_mvdr = zeros(niter,1);
   
for iter=1:niter
    
    S = sqrt(Ps/2) * as * (randn(1,K)+1i*randn(1,K));
    IN = Aj * diag(sqrt(Pj/2)) * (randn(J,K)+1i*randn(J,K)) + sqrt(sigma2/2)*(randn(N,K)+1i*randn(N,K));
    Y = S + IN;

    %Estimated covariance matrices
    Rsmi = (Y*Y')/K; 
    Csmi = (IN*IN')/K; 

    %Weight vector direct form
    %MPDR
    w_mpdr_df= Rsmi\a0; 
    w_mpdr_df = w_mpdr_df / (a0'*w_mpdr_df);
    sinr_mpdr_df(iter,1) = Ps * abs(w_mpdr_df'*as)^2 / abs(w_mpdr_df'*Cth*w_mpdr_df);
    Aw_mpdr(iter) = 1 / (norm(w_mpdr_df)^2);

    %MVDR
    w_mvdr_df= Csmi\a0; 
    w_mvdr_df = w_mvdr_df / (a0'*w_mvdr_df);
    sinr_mvdr_df(iter,1) = Ps * abs(w_mvdr_df'*as)^2 / abs(w_mvdr_df'*Cth*w_mvdr_df);	
    Aw_mvdr(iter) = 1 / (norm(w_mvdr_df)^2);

  
    %GSC-MPDR
    d = w_cbf' * Y; %signal in main channel 1|K
    Z = B' * Y;         %signal in auxilliary channels N-1|K
    Rz = (Z*Z')/K; rdz = Z*d'/K;
    wa= Rz\rdz;
    w_mpdr_gsc = w_cbf - B * wa;
    sinr_mpdr_gsc(iter,1) = Ps*abs(w_mpdr_gsc'*as)^2 / abs(w_mpdr_gsc'*Cth*w_mpdr_gsc);
  

    %GSC-MVDR
    d = w_cbf' * IN; %signal in main channel 1|K
    Z = B' * IN;         %signal in auxilliary channels N-1|K
    Rz = (Z*Z')/K; rdz = Z*d'/K;
    wa= Rz\rdz;
    w_mvdr_gsc = w_cbf - B * wa;
    sinr_mvdr_gsc(iter,1) = Ps*abs(w_mvdr_gsc'*as)^2 / abs(w_mvdr_gsc'*Cth*w_mvdr_gsc);  

end

%Plot SINR
plot((1:niter),sinr_opt*ones(niter,1),'k',(1:niter),10*log10(sinr_mpdr_df),'r',(1:niter),10*log10(sinr_mvdr_df),'b')
ylabel('dB','FontSize',12);
title('Comparison SINR MPDR-MVDR','FontSize',14);
legend('opt','MPDR','MVDR');

%White noise gain
figure
plot((1:niter),10*log10(N)*ones(niter,1),'k',(1:niter),10*log10(Aw_mpdr),'r',(1:niter),10*log10(Aw_mvdr),'b')
ylabel('dB','FontSize',12);
title('White noise array gain MPDR-MVDR','FontSize',12);
legend('CBF','MPDR','MVDR');


%Diagrams
tab_theta = (-90:0.5:90)/180*pi;        %Angles where to evaluate beampatterns
A = exp(1i*2*pi*pos*sin(tab_theta));    %Steering matrix
G_cbf = 20*log10(abs(w_cbf'*A));
G_mpdr = 20*log10(abs(w_mpdr_df'*A));
G_mvdr = 20*log10(abs(w_mvdr_df'*A));
figure
plot(tab_theta*180/pi,G_cbf,'g-.',tab_theta*180/pi,G_mpdr,'b--',tab_theta*180/pi,G_mvdr,'r:','linewidth',2);
title('Comparison CBF-MPDR-MVDR','fontsize',14);
ylabel('dB','FontSize',12);
xlabel('Angle of Arrival (degrees)','fontsize',12);
legend('CBF','MPDR','MVDR');
axis([-90 90 -70 10]);


%SINR
mean_sinr_mpdr_df=10*log10(mean(sinr_mpdr_df));
mean_sinr_mvdr_df=10*log10(mean(sinr_mvdr_df));
mean_sinr_mpdr_gsc=10*log10(mean(sinr_mpdr_gsc));
mean_sinr_mvdr_gsc=10*log10(mean(sinr_mvdr_gsc));

disp(['SINR_opt = ',num2str(sinr_opt),' dB']);
disp(['SINR_cbf = ',num2str(sinr_cbf),' dB']);
disp(['SINR_mpdr_df = ',num2str(mean_sinr_mpdr_df),' dB']);
disp(['SINR_mvdr_df = ',num2str(mean_sinr_mvdr_df),' dB']);
disp(['SINR_mpdr_gsc = ',num2str(mean_sinr_mpdr_gsc),' dB']);
disp(['SINR_mvdr_gsc = ',num2str(mean_sinr_mvdr_gsc),' dB']);

%% Application de la methode 1
r = 10;
w_methode1 = formateurav( Rth, a0, r );
for i = 1:r
    SINR_methode1(i) = 10*log10(Ps * abs(a0'*w_methode1(:,i).^2 / abs(w_methode1(:,i)'*Rth*w_methode1(:,i))));
end
figure
plot(SINR_methode1);
title('SINR methode1')

%% Application de la methode 1
w_methode2 = gradconj( Rth, a0, r );
for i = 1:r
    SINR_methode2(i) = 10*log10(Ps * abs(a0'*w_methode2(:,i).^2 / abs(w_methode2(:,i)'*Rth*w_methode2(:,i))));
end
figure
plot(SINR_methode2);
title('SINR methode1')
