% code for HW 4 Presentation
% heat flow due to a plume

clear all;close all;

plt=1; % do you want to plot or not

n = 2^11;
Lx = 1e7; % length of x domain in m
Ly = 1e7; % length of x domain in m
Lz = 1e5; % length of z domain in m
x = linspace(-Lx/2,Lx/2,n);
y = linspace(-Ly/2,Ly/2,n);
z = linspace(0,Lz,n);

dx = Lx/n; %???

kx = (-n/2:n/2-1)/Lx;
ky = (-n/2:n/2-1)/Ly;

vx = 10/100/(365*24*60*60); % cm/yr to m/s 
vy = 0; %assuming 1d plate motion
sigma = 50^3/(2*sqrt(2*log(2))); % m
A = 0.948; % W m^-2
k = 3.3; % W m^-2
kappa = 8e-7; %m^2 s^-1


% define heat source function
q = A*exp(-(x.*x + y.*y)/(2*sigma^2)); % no delta function?

% plot heat source function
if plt == 1
    figure
    lw = 2;
    fs = 16;
    plot(x,q,'LineWidth',lw)
    xlim([-1e6 1e6])
    ylim([0 1])
    xlabel('x (m)','FontSize',fs)
    ylabel('q (W m^2)','FontSize',fs)
    title('Heat source','FontSize',fs)
end

% take fourier transform of heat source fcn
Q = fftshift(fft2(fftshift(q)));

% plot fft of heat source fcn
figure
plot(kx,Q,'LineWidth',lw)
%xlim([-1e6 1e6])
%ylim([0 1])
xlabel('kx','FontSize',fs)
ylabel('Q','FontSize',fs)
title('Heat source in wavenumber domain','FontSize',fs)


% define p
p = sqrt( (kx.*kx + ky.*ky) - 1i*(vx.*kx+vy.*ky)/(2*pi*kappa));

%remove singularity in p
sing = find(p==0);
p(sing) = p(sing+1); % I don't like this. 
                     % I'm going to think of something better


% z0 = 1i*p; % so many questions about this???

% define temp function in wavenumber domain
T = Q.*exp(-2*pi.*p.*z)./(2*1i.*p);
% take inverse ft to get temp function in real domain
t = ifft2(fftshift(T));

plot(x,t)

