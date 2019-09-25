function [X,Y] = GetSample(N,D)

% function [X,Y] = GetSample(N,D)
%
% Get N paired samples (see details below)
%
% Inputs:
%
%  N - number of samples (Nx=Ny=N)
%  D - the "true" difference in the means (optional)
%
% Outputs:
%
%  X,Y - Nx1 vectors with random samples
%        drawn from the following distribution:
%
%       X(i) ~  N(mx,Sx^2),       mx=4, Sx=1
%       Y(i) ~  N(X(i)+D,Sd^2),   Sd=.5
%
% When no D input is given, D=1 and the random number
% generator is reset

if(nargin<2) D=.5; randn('seed',1); end

mx =  4;
Sx =  1;
Sd = .8;

X = mx +     randn(N,1)*Sx;
Y = X  + D + randn(N,1)*Sd;
