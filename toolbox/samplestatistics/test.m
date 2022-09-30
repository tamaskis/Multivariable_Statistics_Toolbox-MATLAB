clear;clc;close all;

n = 20;

mu = zeros(n,1);
Sigma = eye(n);

x = randmvn(mu,Sigma,1000);

tic
for i = 1:100
    mu = sample_mean_old(x);
end
toc

tic
for i = 1:100
    mu = sample_mean_new(x);
end
toc

mu = sample_mean_old(x)
mu = sample_mean_new(x)

%X = [1, 2, 3;
%     4, 5, 6;
%     7, 8, 9];
%X*ones(3,1)