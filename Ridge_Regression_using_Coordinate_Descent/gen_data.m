function [X,Y] = gen_data(mu, sigma_data,d, n, K, sigma_noise, B)
    X = normrnd(mu,sigma_data,d,n);
    W = zeros(d,1);
    W(1:K,1) = 10;
    E = normrnd(0,sigma_noise, n,1);
    Y = X'*W + E + B;
end