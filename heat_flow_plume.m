% code for HW 4 Presentation
% heat flow due to a plume

clear all; close all;

nx = 4096; % grid refinement
ny = 1024*4;
Lx = 4e6; % length of x domain in m
Ly = 4e6; % length of y domain in m
x = linspace(-Lx/2,Lx/2,nx);
y = linspace(-Ly/2,Ly/2,ny);

kx = (-nx/2:nx/2-1)/Lx;
ky = (-ny/2:ny/2-1)/Ly;

figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
colorbar('eastoutside')
caxis([0 500]);
title('Temperature distribution over a plume')

vel =[0.1 1 5 10];
for i = 1:4
    vx = vel(i)/100/(365*24*60*60); % cm/yr to m/s 
    vy = 0; %
    sigma = 50^3/(2*sqrt(2*log(2))); % m
    A = 0.948; % W m^-2
    k = 3.3; % W m^-2
    kappa = 8e-7; %m^2 s^-1

    % make the domain and wavenumber domain mesh
    [X,Y] = meshgrid(x,y);
    [KX,KY] = meshgrid(kx,ky);

    % define heat source function
    q = A*exp(-(X.*X + Y.*Y)/(2*sigma^2)); % no delta function?

    % plot heat source function
    figure
    lw = 2;
    fs = 16;
    imagesc(x,y,abs(q))
    xlim([-1e6 1e6])
    ylim([-1e6 1e6])
    axis square 
    colormap hot
    xlabel('x (m)','Interpreter','latex','FontSize',fs)
    ylabel('y (m)','Interpreter','latex','FontSize',fs)
    title('Heat source viewed from above','Interpreter','latex','FontSize',fs)

    figure
    plot(x,q(ny/2,:),'LineWidth',lw)
    xlim([-1e6 1e6])
    ylim([0 1])
    xlabel('x (m)','Interpreter','latex','FontSize',fs)
    ylabel('q ($W m^2$)','Interpreter','latex','FontSize',fs)
    title('Gaussian heat source','FontSize',fs,'Interpreter','latex')



    %% take fourier transform of heat source fcn
    Q = fftshift(fft2(fftshift(q)));

    % plot fft of heat source fcn
    figure
    plot(kx,Q(ny/2,:),'LineWidth',lw)
    xlim([-1e-4 1e-4])
    ylim([0 18000])
    xlabel('kx','Interpreter','latex','FontSize',fs)
    ylabel('Q','Interpreter','latex','FontSize',fs)
    title('Heat source in wavenumber domain','FontSize',fs,'Interpreter','latex')


    figure
    imagesc(x,y,real(Q))
    axis square 
    colormap hot
    xlim([-1e4 1e4])
    ylim([-1e4 1e4])
    xlabel('x (m)','FontSize',fs,'Interpreter','latex')
    ylabel('y (m)','FontSize',fs,'Interpreter','latex')
    title('FT of heat source, viewed from above','FontSize',fs,'Interpreter','latex')

    %% define p
    p = sqrt( (KX.*KX + KY.*KY) + 1i*(vx.*KX+vy.*KY)/(2*pi*kappa));

    %% remove singularity in p
    sing = find(p==0);
    p(sing) = p(sing+1); % I don't like this. 

    % z0 = 0; % so many questions about this???
    % [FX,Z] = meshgrid(x,z);
    % z0 = zeros(n,n);

    z0 = 8e3;
    z = z0/2;

    % define temp function in wavenumber domain
    % T = Q.*exp(-2*pi*p*abs(z0-z))- exp(-2*pi*p*abs(z0+z))./(2*1i*p);
    T = Q.*(exp(-2*pi*p*abs(z0-z))- exp(-2*pi*p*abs(z0+z)))./(4*pi*k*p);

    figure
    pcolor(KX,KY,real(T))
    colorbar
    title('Temperature in wavenumber domain')
    shading flat

    % take inverse ft to get temp function in real domain
    t = ifftshift(ifft2(fftshift(T)));
    
    
    figure(1)
    subplot(2,2,i)
    pcolor(X,Y,real(t))
    shading flat
    caxis([0 500])
    title(['$v_x = $ ' num2str(vel(i)) ' cm/yr'],'Interpreter','latex','FontSize',fs)
    xlim([-1e6 1e6])
    ylim([-0.5e6 0.5e6])
    ax.TickLabelInterpreter='latex';
    xlabel('x (m)','Interpreter','latex','FontSize',fs)
    ylabel('y (m)','Interpreter','latex','FontSize',fs)

end

figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
hp4 = get(subplot(2,2,4),'Position');
colorbar('Position', [[0.92, 0.11, 0.02, 0.815]])
caxis([0 500])
