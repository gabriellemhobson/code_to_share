% code for HW 4 Presentation
% heat flow due to a plume

n = 2^12;
x = linspace(-5e6,5e6,n);
y = linspace(-5e6,5e6,n);
z = linspace(0,1e5,n);

kx = -n/2:n/2-1;
ky = -n/2:n/2-1;

vx = 1e-10; % m/s made up value
vy = 0; %assuming 1d plate motion
sigma = 50^3/(2*sqrt(2*log(2))); % m
A = 0.948; % W m^-2
k = 3.3; % W m^-2
kappa = 8e-7; %m^2 s^-1


% define heat source function
q = A*exp(-(x.*x + y.*y)/(2*sigma^2)); % no delta function?

% plot heat source function
figure
plot(x,q)



% take fourier transform of heat source fcn
Q = fftshift(fft2(q));
% define p
p = sqrt( (kx.*kx + ky.*ky) - 1i*(vx.*kx+vy.*ky)/(2*pi*kappa));

%remove singularity in p
sing = find(p==0);
p(sing) = p(sing+1);

% z0 = 1i*p; % so many questions about this???

% define temp function in wavenumber domain
T = Q.*exp(-2*pi.*p.*z)./(2*1i.*p);
% take inverse ft to get temp function in real domain
t = ifft2(fftshift(T));

plot(x,t)

