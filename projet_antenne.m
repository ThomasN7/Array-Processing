%--------------------------------------------------------------------------
%   Projet Traitement d'antennes
%   Pierre Apap - Thomas Nicot - Stanislas Le Grelle
%   3 EN
%--------------------------------------------------------------------------

%% Application de la methode 1
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

%Looked angle
theta0 = 0/180*pi;
a0 = exp(1i*2*pi*pos*sin(theta0));

%Interferences
thetaj = [-20;15]/180*pi;	%angles of arrival	
INR = [20;15];			%interference to noise ratio
Pj = sigma2 * 10.^(INR/10);		%interference power
J = length(thetaj);

%Interference + noise covariance matrix
Aj = exp(1i*2*pi*pos*sin(thetaj'));	%interference steering matrix N|J
Cth = Aj*diag(Pj)*Aj' + sigma2*eye(N);	%covariance matrix

%Total covariance matrix
Rth = Ps*(as*as') + Cth;

%Plot methode 1
r = 10;
w_methode1 = formateurav( Rth, a0, r );
SINR_methode1 = zeros(1,r);
for i = 1:r
    SINR_methode1(i) = 10*log10(Ps * abs(w_methode1(:,i)'*as)^2 / abs(w_methode1(:,i)'*Cth*w_methode1(:,i)));
end
figure
plot(SINR_methode1);
title('SINR methode1')

% Gain AV
tab_theta = (-90:0.5:90)/180*pi;        %Angles where to evaluate beampatterns
A = exp(1i*2*pi*pos*sin(tab_theta));    %Steering matrix
G_methode1 = 20*log10(abs(w_methode1(:,9)'*A));
figure
plot(tab_theta*180/pi,G_methode1,'b:')
hold on

% Gain MPDR
w_mpdr_df= (Rth\a0); 
w_mpdr_df = w_mpdr_df / (a0'*w_mpdr_df);
sinr_mpdr_df = 10*log10( Ps * abs(w_mpdr_df'*as)^2 / abs(w_mpdr_df'*Cth*w_mpdr_df) );
G_mpdr = 20*log10(abs(w_mpdr_df'*A));
plot(tab_theta*180/pi,G_mpdr,'r-.')
legend('Formateur AV','MPDR')
text(-80,0,'interferences à -20° et 15°')

%% Application de la methode 2

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

%Looked angle
theta0 = 0/180*pi;
a0 = exp(1i*2*pi*pos*sin(theta0));

% Interference 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetaj = [-20;15]/180*pi;	%angles of arrival	
INR = [20;15];			%interference to noise ratio
Pj = sigma2 * 10.^(INR/10);		%interference power
J = length(thetaj);

%Interference + noise covariance matrix
Aj = exp(1i*2*pi*pos*sin(thetaj'));	%interference steering matrix N|J
Cth = Aj*diag(Pj)*Aj' + sigma2*eye(N);	%covariance matrix

%Total covariance matrix
Rth = Ps*(as*as') + Cth;

%Plot methode 2
%r = ceil(size(INR,1)/2);
r = 10;
w_methode2 = gradconj( Rth, a0, r );
for i = 1:r
    SINR_methode2(i) = 10*log10(Ps * abs(a0'*w_methode2(:,i))^2 / abs(w_methode2(:,i)'*Cth*w_methode2(:,i)));
end
figure
plot(SINR_methode2);

% Interference 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetaj = [-20;15;50]/180*pi;	%angles of arrival	
INR = [20;15;15];			%interference to noise ratio
Pj = sigma2 * 10.^(INR/10);		%interference power
J = length(thetaj);

%Interference + noise covariance matrix
Aj = exp(1i*2*pi*pos*sin(thetaj'));	%interference steering matrix N|J
Cth = Aj*diag(Pj)*Aj' + sigma2*eye(N);	%covariance matrix

%Total covariance matrix
Rth = Ps*(as*as') + Cth;

%Plot methode 2
%r = ceil(size(INR,1)/2);
r = 10;
w_methode2 = gradconj( Rth, a0, r );
for i = 1:r
    SINR_methode2(i) = 10*log10(Ps * abs(a0'*w_methode2(:,i))^2 / abs(w_methode2(:,i)'*Cth*w_methode2(:,i)));
end
hold on
plot(SINR_methode2,'y');

% Interference 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetaj = [-20;15;50;-45]/180*pi;	%angles of arrival	
INR = [20;15;15;10];			%interference to noise ratio
Pj = sigma2 * 10.^(INR/10);		%interference power
J = length(thetaj);

%Interference + noise covariance matrix
Aj = exp(1i*2*pi*pos*sin(thetaj'));	%interference steering matrix N|J
Cth = Aj*diag(Pj)*Aj' + sigma2*eye(N);	%covariance matrix

%Total covariance matrix
Rth = Ps*(as*as') + Cth;

%Plot methode 2
%r = ceil(size(INR,1)/2);
r = 10;
w_methode2 = gradconj( Rth, a0, r );
for i = 1:r
    SINR_methode2(i) = 10*log10(Ps * abs(a0'*w_methode2(:,i))^2 / abs(w_methode2(:,i)'*Cth*w_methode2(:,i)));
end
hold on
plot(SINR_methode2,'c');

% Interference 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetaj = [-20;15;50;-45;70]/180*pi;	%angles of arrival	
INR = [20;15;15;10;25];			%interference to noise ratio
Pj = sigma2 * 10.^(INR/10);		%interference power
J = length(thetaj);

%Interference + noise covariance matrix
Aj = exp(1i*2*pi*pos*sin(thetaj'));	%interference steering matrix N|J
Cth = Aj*diag(Pj)*Aj' + sigma2*eye(N);	%covariance matrix

%Total covariance matrix
Rth = Ps*(as*as') + Cth;

%Plot methode 2
%r = ceil(size(INR,1)/2);
r = 10;
w_methode2 = gradconj( Rth, a0, r );
for i = 1:r
    SINR_methode2(i) = 10*log10(Ps * abs(a0'*w_methode2(:,i))^2 / abs(w_methode2(:,i)'*Cth*w_methode2(:,i)));
end
hold on
plot(SINR_methode2,'r');
legend([int2str(size(INR,1)) ' interferences'])

% Interference 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetaj = [-20;30;15;50;-45;70]/180*pi;	%angles of arrival	
INR = [20;10;15;15;10;25];			%interference to noise ratio
Pj = sigma2 * 10.^(INR/10);		%interference power
J = length(thetaj);

%Interference + noise covariance matrix
Aj = exp(1i*2*pi*pos*sin(thetaj'));	%interference steering matrix N|J
Cth = Aj*diag(Pj)*Aj' + sigma2*eye(N);	%covariance matrix

%Total covariance matrix
Rth = Ps*(as*as') + Cth;

%Plot methode 2
%r = ceil(size(INR,1)/2);
r = 10;
w_methode2 = gradconj( Rth, a0, r );
for i = 1:r
    SINR_methode2(i) = 10*log10(Ps * abs(a0'*w_methode2(:,i))^2 / abs(w_methode2(:,i)'*Cth*w_methode2(:,i)));
end
hold on
plot(SINR_methode2,'g');
legend([int2str(size(INR,1)) ' interferences'])

% Interference 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetaj = [-20;-70;30;15;50;-45;70]/180*pi;	%angles of arrival	
INR = [20;15;10;15;15;10;25];			%interference to noise ratio
Pj = sigma2 * 10.^(INR/10);		%interference power
J = length(thetaj);

%Interference + noise covariance matrix
Aj = exp(1i*2*pi*pos*sin(thetaj'));	%interference steering matrix N|J
Cth = Aj*diag(Pj)*Aj' + sigma2*eye(N);	%covariance matrix

%Total covariance matrix
Rth = Ps*(as*as') + Cth;

%Plot methode 2
%r = ceil(size(INR,1)/2);
r = 10;
w_methode2 = gradconj( Rth, a0, r );
for i = 1:r
    SINR_methode2(i) = 10*log10(Ps * abs(a0'*w_methode2(:,i))^2 / abs(w_methode2(:,i)'*Cth*w_methode2(:,i)));
end
hold on
plot(SINR_methode2,'k');
legend([int2str(size(INR,1)) ' interferences'])

% Interference 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetaj = [-20;-70;30;15;50;-45;70;-60]/180*pi;	%angles of arrival	
INR = [20;15;10;15;15;10;25;10];			%interference to noise ratio
Pj = sigma2 * 10.^(INR/10);		%interference power
J = length(thetaj);

%Interference + noise covariance matrix
Aj = exp(1i*2*pi*pos*sin(thetaj'));	%interference steering matrix N|J
Cth = Aj*diag(Pj)*Aj' + sigma2*eye(N);	%covariance matrix

%Total covariance matrix
Rth = Ps*(as*as') + Cth;

%Plot methode 2
%r = ceil(size(INR,1)/2);
r = 10;
w_methode2 = gradconj( Rth, a0, r );
for i = 1:r
    SINR_methode2(i) = 10*log10(Ps * abs(w_methode2(:,i)'*as)^2 / abs(w_methode2(:,i)'*Cth*w_methode2(:,i)));
end

hold on
plot(SINR_methode2,'m');
legend('2 interférences','3 interférences','4 interférences','5 interférences','6 interférences','7 interférences','8 interférences')
title('Algorithme du gradient conjugué')
xlabel('Nombre d''itérations')
ylabel('SINR')


