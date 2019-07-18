function [foury,Foury] = FouryComp(s,t,u,mu)

[L,K] = size(t);
Foury = zeros(length(s),L);
for l = 1:L
    for k = 1:K
        Foury(:,l) = Foury(:,l)+u(l,k)*sqrt(2*pi)*mu(l)*exp(-2*pi^2*mu(l)^2*s.^2).*exp(2*1i*pi*s*t(l,k));
    end
end

foury = sum(Foury,2);


