function [foury] = FouryCompSingle(s,t,u,mu,K)

foury = zeros(length(s),1);

for k = 1:K
        foury = foury + sqrt(2*pi)*mu*exp(-2*pi^2*mu^2*s.^2).*u(k).*exp(2*1i*pi*s*t(k));
end



