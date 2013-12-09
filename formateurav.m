function [ w ] = formateurav( R, a0, r )
%FORMATEURAV    Methode iterative de formateur de faisceaux w
%   Formateurs de faisceaux w_n generee par recursion
%   [ w ] = formateurav( R, a0, r )
%   R est la mtrice de covariance, a0 le vecteur directeur, r le nombre
%   d'iterations de l'agorithme

taille = length(a0);
w = zeros(taille, r+1);
w(:,1) = a0/(a0'*a0);
Proj = eye(taille) - a0*a0' ./ (a0'*a0);
for i = 2:r+1
    g = Proj*R*w(:,i-1);
    if g == zeros(taille, 1)
        return %break
    else
        mu = (g'*R*w(:,i-1))/(g'*R*g);
        w(:,i) = w(:,i-1) - mu*g;
    end
end
% fonction non verifee
end