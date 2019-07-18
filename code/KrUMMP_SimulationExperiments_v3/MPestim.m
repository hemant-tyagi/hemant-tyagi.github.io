function [testim,uestim] = MPestim(fourytilde,s0,m,K)

%---------------------------------
% Create Toeplitz matrices H0, H1
%---------------------------------
r0 = [fourytilde(m+1:2*m)];
c0 = [fourytilde(m+1:-1:2)];
H0 = toeplitz(c0,r0);
r1 = [fourytilde(m:2*m-1)];
c1 = [fourytilde(m:-1:1)];
H1 = toeplitz(c1,r1);

%---------------------------------------------------------------------------
% Project H0,H1 onto top K left singular vectors of H0,H1 and 
% find generalized eigenvalues of H1, H0. Finally extract spike locations
% from gen. eig. values.
%---------------------------------------------------------------------------
[U,S,V]=svd(H0,0); % svd of H0
U = U(:,1:K); % m x K left singular matrix
Ahat = U'*H0*U; Bhat = U'*H1*U; % K x K matrices
gen_eigvals = eig(Bhat,Ahat); % find gen eig values of (Bhat,Ahat). Note: Bhat.v = lambda.Ahat.v, jere lambda is generalized eigenvalue.
gen_eigvals = gen_eigvals./abs(gen_eigvals); % project on to unit circle
[theta,~] = cart2pol(real(gen_eigvals), imag(gen_eigvals)); % Cartesian to Polar coordinates
theta1 = mod(-theta,2*pi); 
testim = theta1./(2*pi); % Extract the estimated spike locations

%---------------------------------
% Now estimate spike amplitudes
%---------------------------------
Vhat = zeros(m,K);
for i = 1:m
   Vhat(i,:) = (exp(-2*1i*pi*testim).^(i-1))'; % Form Vandermonde matrix using estimated {t_{j,k}}
end
uestim = pinv(Vhat)*(r0); % And now obtain the amplitude estimates as least squares solution
uestim = uestim.*exp(-2*1i*pi*s0*testim); % Take care of shift due to s0


