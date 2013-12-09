function [ w ] = gradconj( R, a0, r )
%GRADCONJ   Resolution par la methode du gradient conjugue
%   Detailed explanation goes here
%   [ w ] = gradconj( R, a0, r )
%   R est la mtrice de covariance, a0 le vecteur directeur, r le nombre
%   d'iterations de l'agorithme

taille = length(a0);
% Omega
w = zeros(taille, r);
% Beta
beta = zeros(r,1); %beta(1) = 0;
% U
u = zeros(taille,r);
u(:,1) = a0;
% E
e = zeros(taille,r);
% Pour i = 1
    i = 1;
    z = R*u(:,i);
    c = norm(a0)^2 / (u(:,i)'*z); % car e0 = a0
    e(:,i) = a0 - c*z;
    w(:,i) = zeros(r,1) + c*u(:,i); % car w0 = 0
% Pour i = 2
    i = 2;
    beta(i) = norm(e(:,i-1))^2 / norm(a0)^2; % car e0 = a0
    u(:,i) = e(:,i-1)+ beta(i)*u(:,i-1);
    z = R*u(:,i);
    c = norm(e(:,i-1))^2 / (u(:,i)'*z);
    e(:,i) = e(:,i-1) - c*z;
    w(:,i) = w(i-1) + c*u(:,i);
for i = 3:r % Pour i allant de 3 a r
    %
    beta(i) = norm(e(:,i-1))^2 / norm(e(:,i-2))^2;
    u(:,i) = e(:,i-1)+ beta(i)*u(:,i-1);
    %
    z = R*u(:,i);
    c = norm(e(:,i-1))^2 / (u(:,i)'*z);
    e(:,i) = e(:,i-1) - c*z;
    w(:,i) = w(i-1) + c*u(:,i);
end     
% fonction non verifiee
end