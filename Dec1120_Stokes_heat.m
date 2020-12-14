%advection 2D
% implement sticky air

clear all
% close all
L = 500*1e3; %length of x domain
D = 500*1e3; %length of z domain
Air = 100*1e3; %thickness of free surface 'sticky air'
Eta_amb = 1e21; %ambient viscosity of domain [Pas]
Rho_amb = 3300; %ambient density of domain [kg/m^3]
Eta_plume = 1e20; %viscosity of plume
Rho_plume = 3200; %density of plume
Eta_air = 1e17; %viscosity of sticky air
Rho_air = 1; %density of sticky air
Ta = 1300; %Adiabatic temperature [k]
Alpha = 3*1e-5;
% %====================% %
%  set numerical grids
% %====================% %
%set eulerian grid
nx = 51; %number of x nodes
nz = 61; %number of z nodes
nx1 = nx+1; %plus ghost nodes
nz1 = nz+1; %plus ghost nodes
dx = L/(nx-1); %spacing of the x doordinates
dz = (D+Air)/(nz-1); %spacing of the z coordinates
gz = 10; %gravity acceleration
gx = 0; %no horizontal gravity
% Rho_all = zeros(nz,nx);
% Eta_all = zeros(nz,nx);
[x,z,xvx,zvx,xvz,zvz,xp,zp] = vectorsetup(D+Air,L,dx,dz); %get vectors
[x2d,z2d] = meshgrid(x,z); %create grid
%setup initial visc and density profile
rplume = 100000;

[indP,indvx,indvz] = index_setup(numsetup(nz1,nx1));

%setup marker grid
nxm = 4; %number of x point clouds within node
nzm = 4; %number of z point clouds within node
dxm = dx/nxm; %spacing between each x marker
dzm = dz/nzm; %spacing between each z marker
nxm_all = nxm*(nx-1); %total number of x markers in one row
nzm_all = nzm*(nz-1); %total number of z markers in one collumn
ztemp = 0:dz:Air+D;
marknum = nxm_all*nzm_all; %total number of markers
% setup initial material grid onto markers
xm=zeros(1,marknum); % Horizontal coordinates, m
zm=zeros(1,marknum); % Vertical coordinates, m
matm = zeros(1,marknum); %material properties
Tm = zeros(1,marknum); % Marker temperature K
m=1; % Marker counter
% above surface

% Define properties of materials: 
% %       1mantle 2plume 3air
rhom   = [3300   3200   1     ]; % Density, kg/m^3
etam   = [1e+21  1e+20  1e+18 ]; % Viscosity, Pa s
cpm = [3.3e+6 3.2e+6 3.3e+6]; % Volumetric heat capacity, kg/m^3
% alpham = [3e-5   2e-5   0     ]; % Thermal expansion, 1/K
km     = [3      2      3000  ]; % Thermal conductivity, W/m/K
% hrm    = [2e-8   3e-8   0     ]; % Radiogenic heat production, W/m^3


for jm=1:1:nxm_all
    for im=1:1:nzm_all
        % Define marker coordinates
        xm(m)=dxm/2+(jm-1)*dxm+(rand-0.5)*dxm;
        zm(m)=dzm/2+(im-1)*dzm+(rand-0.5)*dzm;
        % Marker properties
        rmark=((xm(m)-L/2)^2+(zm(m)-(D+Air)/2)^2)^0.5;
        if(zm(m)>0.9*(Air+D))
%         if rmark<rplume
            matm(m) = 2; %plume material
            Tm(m) = 1800; % plume temperature
        else
            matm(m) = 1; %mantle material
            Tm(m) = 1500; %mantle temperature
        end
        % Sticky air (to have internal free surface)
        if(zm(m)<Air)
            matm(m) = 3; %air material
            Tm(m) = 298; %air temperature [k]
        end
        % Update marker counter
        m=m+1;
    end
end

%initialise arrays
Epsxz = zeros(nz,nx); %strain rate on the ordinary grid
Sigxz = zeros(nz,nx); %deviatoric stress on the ordinary grid
Epsxx = Epsxz; %strain rate in 
Sigxx = Sigxz; %deviatoric stress the middle of grid/pressure nodes
Hs = Sigxx; %shear heating, on the pressure nodes
Ha = Hs; %adiabatic heating, on pressure nodes

% Boundary conditions: free slip=-1; No Slip=1
bcleft=-1;
bcright=-1;
bctop=-1;
bcbottom=-1;

%initialise numerical loop
time = 0;
nt = 200;
vpratio=1/3; % Weight of averaged velocity for moving markers
dt = 1e10;
dxzmax = 0.5;
for ti = 1:nt
%     time = time+dt;
    
    %interpolate marker to real grid
    [Eta_out,Eta_mid,Rho_vz,Rho_vx,Rho_mid,k_vz,k_vx,CP_mid,T1]...
    = marker2grid(nx,nz,nx1,nz1,marknum,matm,xm,zm,dx,dz,x,z,xvx,zvx,...
    xvz,zvz,xp,zp,Tm,rhom,etam,km,cpm);

    Pscale = Eta_amb/(dx+dz)*2; %pressure scaling coefficient, minimum viscosity    
    
    %compute stokes continuity
    [P_out,vx_out,vz_out] = stokes_continuity(nx,nz,nx1,nz1,indvx,indvz...
    ,indP,Eta_out,Eta_mid,Pscale,gz,dx,dz,bctop,bcbottom,bcleft,bcright,Rho_vz,dt,gx,Rho_vx);

    % Calculate deviatoric stress and strain components
    for j = 1:1:nx
    for i = 1:1:nz
        Epsxz(i,j) = 1/2 * ((vx_out(i+1,j)-vx_out(i,j))/dz...
            + (vz_out(i,j+1)-vz_out(i,j))/dx); %shear strain on normal grid
        Sigxz(i,j) = 2*Eta_out(i,j)*Epsxz(i,j); %deviatoric stress on normal grid
    end 
    end
    for j = 2:1:nx %calculate normal stress/strain components on the middle of grid/pressure nodes
    for i = 2:1:nz
        Epsxx(i,j) = (vx_out(i,j)-vx_out(i,j-1))/dx; %Normal strain rate
        Sigxx(i,j) = 2*Eta_mid(i,j)*Epsxx(i,j); %deviatoric stress
    end
    end
    
    %compute shear heating and adiabatic heating
    for j = 2:1:nx 
    for i = 2:1:nz
        Hxz = (Epsxz(i-1,j-1)*Sigxz(i-1,j-1) + Epsxz(i-1,j)*Sigxz(i-1,j) +...
            Epsxz(i-1,j)*Sigxz(i-1,j) + Epsxz(i,j)*Sigxz(i,j))/4; %average of xz products
        Hs(i,j) = 2*Hxz + 2*Epsxx(i,j)*Sigxx(i,j); %shear heating
        % now compute adiabatic heating
        Ha(i,j) = (vz_out(i,j)+vz_out(i-1,j))/2 * (Rho_vz(i,j)+Rho_vz(i-1,j))/2 ...
            *Ta*gz*Alpha;
    end
    end
       


    % averaging velocities on centre nodes
    vx_mid = zeros(nz1,nx1);
    vz_mid = zeros(nz1,nx1);
    for j = 2:1:nx %solve for ordinary nodes only
    for i = 2:1:nz
        vx_mid(i,j) = 1/((1/vx_out(i,j)+1/vx_out(i,j-1))/2);% vx; (current+left)/2
        vz_mid(i,j) = 1/((1/vz_out(i,j)+1/vz_out(i-1,j))/2);% vz; (current+above)/2
    end
    end
    %applying free-slip boundary conditions
    %Top
    vx_mid(1,2:nx-1) = -bctop*vx_mid(2,2:nx-1);
    vz_mid(1,:)      = -vz_mid(2,:);
    %bottom
    vx_mid(nz1,2:nx-1) = -bcbottom*vx_mid(nz,2:nx-1);
    vz_mid(nz1,:)      = -vz_mid(nz,:);
    %left
    vx_mid(:,1)       = -vx_mid(:,2);
    vz_mid(2:nz-1,1)  = -bcleft*vz_mid(2:nz-1,2);
    %right
    vx_mid(:,nx1)=-vx_mid(:,nx);
    vz_mid(2:nz-1,nx1)=-bcright*vz_mid(2:nz-1,nx); % Free slip    

    % Define timestep
    dt=1e+30;
    maxvx=max(max(abs(vx_out)));
    maxvz=max(max(abs(vz_out)));
    if(dt*maxvx>dxzmax*dx)
        dt=dxzmax*dx/maxvx;
    end
    if(dt*maxvz>dxzmax*dz)
        dt=dxzmax*dz/maxvz;
    end
    vxm=zeros(4,1);
    vzm=zeros(4,1);
    %advect the marker coordinates. Eq. 8.19
    % interpolate grid velocities into markers
    for m=1:1:marknum
        %compute 4th order Runge-kutta velocities
        %create temporary variables for markers
        mat = matm(m);
        xmRK = xm(m);
        zmRK = zm(m);
        for rk = 1:1:4
        %interpolate vx_mid and vz_mid
        j=fix((xmRK-xp(1))/dx)+1;
        i=fix((zmRK-zp(1))/dz)+1;
        if(j<1)
            j=1;
        elseif(j>nx)
            j=nx;
        end
        if(i<1)
            i=1;
        elseif(i>nz)
            i=nz;
        end
        %compute distances
        dxm1 = abs(xmRK-xp(j));
        dzm1 = abs(zmRK-zp(i));
        %compute weights
        wtij = (1-dxm1/dx)*(1-dzm1/dz);
        wtij1 = (dxm1/dx)*(1-dzm1/dz);
        wti1j = (1-dxm1/dx)*(dzm1/dz);
        wti1j1 = dxm1*dzm1/dx/dz;
        %compute velocity on node i,j
        vxm(rk) = vx_mid(i,j)*wtij+vx_mid(i,j+1)*wtij1+vx_mid(i+1,j)*wti1j+...
            vx_mid(i+1,j+1)*wti1j1;
        vzm(rk) = vz_mid(i,j)*wtij+vz_mid(i,j+1)*wtij1+vz_mid(i+1,j)*wti1j+...
            vz_mid(i+1,j+1)*wti1j1;
        
        % Interpolate vx
        % Define i,j indexes for the upper left node
        j=fix((xmRK-xvx(1))/dx)+1;
        i=fix((zmRK-zvx(1))/dz)+1;
        if(j<1)
            j=1;
        elseif(j>nx-1)
            j=nx-1;
        end
        if(i<1)
            i=1;
        elseif(i>nz)
            i=nz;
        end
        %calculate distances
        dxm1 = xmRK-xvx(j);
        dzm1 = zmRK-zvx(i);
        %compute weights of distances from real nodes
        wtij = (1-dxm1/dx)*(1-dzm1/dz);
        wtij1 = (dxm1/dx)*(1-dzm1/dz);
        wti1j = (1-dxm1/dx)*(dzm1/dz);
        wti1j1 = dxm1*dzm1/dx/dz;
        %compute velocity on node i,j
        vxm(rk) =vxm(rk)*vpratio + (1-vpratio)*vx_out(i,j)*wtij...
            +vx_out(i,j+1)*wtij1+vx_out(i+1,j)*wti1j+...
            vx_out(i+1,j+1)*wti1j1;
        
        % Interpolate vz
        % Define i,j indexes for the upper left node
        j=fix((xmRK-xvz(1))/dx)+1;
        i=fix((zmRK-zvz(1))/dz)+1;
        if(j<1)
            j=1;
        elseif(j>nx)
            j=nx;
        end
        if(i<1)
            i=1;
        elseif(i>nz-1)
            i=nz-1;
        end
        %calculate distances
        dxm1 = xmRK-xvz(j);
        dzm1 = zmRK-zvz(i);
        %compute weights of distances from real nodes
        wtij = (1-dxm1/dx)*(1-dzm1/dz);
        wtij1 = (dxm1/dx)*(1-dzm1/dz);
        wti1j = (1-dxm1/dx)*(dzm1/dz);
        wti1j1 = dxm1*dzm1/dx/dz;
        %compute velocity on node i,j
        vzm(rk) =vzm(rk)*vpratio + (1-vpratio)*vz_out(i,j)*wtij+vz_out(i,j+1)*wtij1+vz_out(i+1,j)*wti1j+...
            vz_out(i+1,j+1)*wti1j1;
        
        %change coordinates of xm & zm for upper orders
        if (rk==1 || rk==2)
            xmRK = xm(m) + vxm(rk)*dt/2;
            zmRK = zm(m) + vzm(rk)*dt/2;
        elseif (rk==3)
            xmRK = xm(m) + vxm(rk)*dt;
            zmRK = zm(m) + vzm(rk)*dt;
        end
        end
        %compute 4th order RK effective velocity
        vxm_eff = (vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4))/6;
        vzm_eff = (vzm(1)+2*vzm(2)+2*vzm(3)+vzm(4))/6;
        %move markers
        xm(m) = xm(m)+ dt*vxm_eff;   
        zm(m) = zm(m)+ dt*vzm_eff;
    end
        vx_mid = zeros(nz1,nx1);
vz_mid = zeros(nz1,nx1);
% averaging velocities on centre nodes
for j = 2:1:nx %solve for ordinary nodes only
for i = 2:1:nz
    vx_mid(i,j) = 1/((1/vx_out(i,j)+1/vx_out(i,j-1))/2);% vx; (current+left)/2
    vz_mid(i,j) = 1/((1/vz_out(i,j)+1/vz_out(i-1,j))/2);% vz; (current+above)/2
end
end
time = time +dt; %time elapsed

figure(1);colormap('Jet');clf
subplot(1,3,1)
pcolor(x2d,z2d,log10(Eta_out));
shading flat;
axis ij image;
colorbar
title('colormap of log_1_0Viscosity')
hold on
quiver(xp(3:5:nx1),zp(3:5:nz1),vx_mid(3:5:nz,3:5:nx1),vz_mid(3:5:nz1,3:5:nx1),'k')

subplot(1,3,2)
pcolor(xp(2:end-1),zp(2:end-1),Hs(2:end,2:end));axis ij image;shading interp, colorbar
title('Shear heat distribution')

subplot(1,3,3)
pcolor(xp(2:end-1),zp(2:end-1),Ha(2:end,2:end));axis ij image;shading interp, colorbar
title('Adiabatic heat distribution')

% figure(1);colormap('Jet');clf
% title(['time elapsed ' num2str(time) 'seconds'])
% subplot(2,2,1)
% pcolor(x2d,z2d,log10(Eta_out));% caxis([17 21])
% shading flat;
% axis ij image;
% colorbar
% title('colormap of log_1_0Viscosity')
% hold on
% quiver(xp(3:5:nx1),zp(3:5:nz1),vx_mid(3:5:nz,3:5:nx1),vz_mid(3:5:nz1,3:5:nx1),'k')
% subplot(2,2,2)
% pcolor(xp(2:end-1),zp(2:end-1)...
%     ,P_out(2:end-1,2:end-1)./1e6);axis ij image;shading interp, colorbar
% 
% 
% subplot(2,2,3)
% 
% pcolor(xvx,zvx,vx_out(1:nz1,1:nx));axis ij image;
% shading interp, colorbar
% title(['colormap of vx, time elapsed ',num2str(time/3.154e+7/1e6),' Million years'])
% 
% subplot(2,2,4)
% pcolor(xvz,zvz,vz_out(1:nz,1:nx1));axis ij image;shading interp, colorbar
%         
% figure(2);colormap('Jet');clf
% subplot(2,1,1)
% pcolor(xp(2:end-1),zp(2:end-1),Hs(2:end,2:end));axis ij image;shading interp, colorbar
% title('Shear heat distribution')
% 
% subplot(2,1,2)
% pcolor(xp(2:end-1),zp(2:end-1),Ha(2:end,2:end));axis ij image;shading interp, colorbar
% title('Adiabatic heat distribution')
end

function [indP,indvx,indvz] = index_setup(Number_all)
indP = Number_all.*3; %A matrix indexing of P
indvx = indP-2;
indvz = indP-1;
end

function [x,z,xvx,zvx,xvz,zvz,xp,zp] = vectorsetup(D,L,dx,dz)
%setup x and y vectors for vx, vz & P
x=0:dx:L; % Horizontal coordinates of basic grid points, m
z=0:dz:D; % Vertical coordinates of basic grid points, m
xvx=0:dx:L; % Horizontal coordinates of vx grid points, m
zvx=-dz/2:dz:D+dz/2; % Vertical coordinates of vx grid points, m
xvz=-dx/2:dx:L+dx/2; % Horizontal coordinates of vy grid points, m
zvz=0:dz:D; % Vertical coordinates of vy grid points, m
xp=-dx/2:dx:L+dx/2; % Horizontal coordinates of P grid points, m
zp=-dz/2:dz:D+dz/2; % Vertical coordinates of P grid points, m
end


function [P_out,vx_out,vz_out] = stokes_continuity(nx,nz,nx1,nz1,indvx,indvz...
    ,indP,visc_all,visc_mid,k,gz,dx,dz,bctop,bcbottom,bcleft,bcright,Rho_vz,dt,gx,Rho_vx)
N = nx1*nz1*3; %total number of unknowns to solve (vx,vz,P for each node)
%We solve implicitly in the form A*c = RHS, where c is the solution
%implemented a sticky air solver on the vz grid
A = sparse(N,N);
RHS = zeros(N,1);

for j = 1:1:nx1
for i = 1:1:nz1
    %solve x equation
    %boundary conditions of vx
    if(i==1 || i==nz1 || j==1 || j==nx || j==nx1)
        A(indvx(i,j),indvx(i,j)) = 1; % A matrix coefficient
        RHS(indvx(i,j)) = 0; % RHS
        %Top boundary
        if (i==1 && j>1 && j<nx)
            A(indvx(i,j),indvx(i+1,j)) = bctop; %only solve for the bottom of the top boundary
        end
        %Bottom boundary    
        if (i==nz1 && j>1 && j<nx)
            A(indvx(i,j),indvx(i-1,j)) = bcbottom; % above the bottom boundary
        end
    % now solve internal points on the real grid
    else
    % A matrix coefficients
    Eta1 = visc_all(i-1,j);
    Eta2 = visc_all(i,j);
    EtaP1 = visc_mid(i,j);
    EtaP2 = visc_mid(i,j+1);
    drhodx = (Rho_vx(i,j+1)-Rho_vx(i,j-1))/2/dx;
    drhodz = (Rho_vx(i+1,j)-Rho_vx(i-1,j))/2/dz;
    
    
    A(indvx(i,j),indvx(i,j-1)) = 2*EtaP1/dx^2;%vx left of current node
    A(indvx(i,j),indvx(i-1,j))     = Eta1/dz^2;%vx  above current node
    A(indvx(i,j),indvx(i,j))       = -2*(EtaP1+EtaP2)/dx^2-...
                                    (Eta1+Eta2)/dz^2 -...
                                    drhodx*gx*dt;%vx current node
    A(indvx(i,j),indvx(i+1,j))     = Eta2/dz^2;%vx below current node
    A(indvx(i,j),indvx(i,j+1)) = 2*EtaP2/dx^2;%vx right of current node
    
    A(indvx(i,j),indvz(i,j)) = -Eta2/dx/dz -drhodz*gx*dt; %vz bottomleft
    A(indvx(i,j),indvz(i,j+1)) = Eta2/dx/dz -drhodz*gx*dt; %vz bottomright      
    A(indvx(i,j),indvz(i-1,j)) = Eta1/dx/dz -drhodz*gx*dt; %vz topleft
    A(indvx(i,j),indvz(i-1,j+1)) = -Eta1/dx/dz -drhodz*gx*dt; %vz topright     
    A(indvx(i,j),indP(i,j))        = k/dx; %P1; current node
    A(indvx(i,j),indP(i,j+1))  = -k/dx; %P2; right of current node
    % RHS
    RHS(indvx(i,j)) = Rho_vx(i,j)*gx;
    end
    
    %solve z equation
    %boundary conditions of vz
    if(j==1 || j==nx1 || i==1 || i==nz || i==nz1)
        A(indvz(i,j),indvz(i,j)) = 1; % A matrix coefficient
        RHS(indvz(i,j)) = 0; % RHS
        %left boundary
        if (j==1 && i>1 && i<nz)
            A(indvz(i,j),indvz(i,j+1)) = bcleft; %solve for right of the leftmost bodes
        end
        %right boundary    
        if (j==nx1 && i>1 && i<nz)
            A(indvz(i,j),indvz(i,j-1)) = bcright; % above the bottom boundary
        end
    %solve internal points    
    else
        
        Eta1 = visc_all(i,j-1);
        Eta2 = visc_all(i,j);
        EtaP1 = visc_mid(i,j);
        EtaP2 = visc_mid(i+1,j);
        drhodx = (Rho_vz(i,j+1)-Rho_vz(i,j-1))/2/dx;
        drhodz = (Rho_vz(i+1,j)-Rho_vz(i-1,j))/2/dz;
    % A matrix coefficients
    A(indvz(i,j),indvz(i,j-1)) = Eta1/dx^2;%vx1 left of current node
    A(indvz(i,j),indvz(i-1,j))     = 2*EtaP1/dz^2;%vx2 above current node
    A(indvz(i,j),indvz(i,j))       = -2*(EtaP1+EtaP2)/dz^2-...
                      (Eta1+Eta2)/dx^2-...
                      drhodz*gz*dt;%vx3 current node
    A(indvz(i,j),indvz(i+1,j))     = 2*EtaP2/(dz^2);%vx4 below current node
    A(indvz(i,j),indvz(i,j+1)) = Eta2/(dx^2);%vx5 right of current node
    A(indvz(i,j),indvx(i,j)) = -Eta2/dx/dz-drhodx*gz*dt/4;% topright
    A(indvz(i,j),indvx(i+1,j)) = Eta2/dx/dz-drhodx*gz*dt/4; %bottomright
    A(indvz(i,j),indvx(i,j-1)) = Eta1/dx/dz-drhodx*gz*dt/4; %topleft
    A(indvz(i,j),indvx(i+1,j-1)) = -Eta1/dx/dz -drhodx*gz*dt/4;
    
    A(indvz(i,j),indP(i,j))        = k/dz; %P1; current node
    A(indvz(i,j),indP(i+1,j))  = -k/dz; %P2; bottom of current node
    % RHS
    RHS(indvz(i,j)) = -gz*Rho_vz(i,j);
    end
    
    %solve P continuity equation
    %boundary points
    % P equation External points
    if(i==1 || j==1 || i==nz1 || j==nx1 ||...
            (i==2 && j==2))
        % Boundary Condition
        % 1*P=0
        A(indP(i,j),indP(i,j))=1; % Left part
        RHS(indP(i,j))=0; % Right part
         % Real BC
        if(i==2 && j==2)
            A(indP(i,j),indP(i,j))=1*k; %Left part
            RHS(indP(i,j))=1e+9; % Right part   
        end
    %now solve internal points    
    else
    %A matrix coefficients    
    A(indP(i,j),indvx(i,j-1))=-1/dx; % left of current node
        A(indP(i,j),indvx(i,j))=1/dx; % current node
        A(indP(i,j),indvz(i-1,j))=-1/dz; % above current node
        A(indP(i,j),indvz(i,j))=1/dz; % below current node
    
    %RHS
    RHS(indP(i,j)) = 0;
    end       
end 
end
c = A\RHS; %get solution vector
%extrapolate into individual matrices
for j = 1:1:nx1
for i = 1:1:nz1
    
    P_out(i,j) = c(indP(i,j))*k; %output pressure
    vx_out(i,j) = c(indvx(i,j)); %output vx
    vz_out(i,j) = c(indvz(i,j)); %output vz
    
end
end
end

function [Eta_out,Eta_mid,Rho_vz,Rho_vx,Rho_mid,k_vz,k_vx,CP_mid,T1]...
    = marker2grid(nx,nz,nx1,nz1,marknum,matm,xm,zm,dx,dz,x,z,xvx,zvx,...
    xvz,zvz,xp,zp,Tm,rhom,etam,km,cpm)
% Interpolate RHO, ETA from markers
%ordinary nodes
ETASUM=zeros(nz,nx);
WTSUM=zeros(nz,nx);
% vz grid
RHOSUM2 = zeros(nz,nx1);
KSUM2 = zeros(nz,nx1);
WTSUM2 = zeros(nz,nx1);
%vx grid
RHOSUM3 = zeros(nz1,nx);
KSUM3 = zeros(nz1,nx);
WTSUM3 = zeros(nz1,nx);
%Pgrid
TSUM = zeros(nz1,nx1);
RHOPSUM = zeros(nz1,nx1);
CPSUM = zeros(nz1,nx1);
ETAPSUM = zeros(nz1,nx1);
WTSUMP = zeros(nz1,nx1);
WTSUMT = zeros(nz1,nx1);


for m=1:1:marknum
    mat = matm(m); % identify material of the marker
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((zm(m)-z(1))/dz)+1;
    if(j<1)
        j=1;
    elseif(j>nx-1)
        j=nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>nz-1)
        i=nz-1;
    end
    %distance between nodes and markers
    dxm1 = abs(xm(m)-x(j));
    dzm1 = abs(zm(m)-z(i));
    %compute weights of distances from real nodes
    wtij = (1-dxm1/dx)*(1-dzm1/dz);
    wtij1 = (dxm1/dx)*(1-dzm1/dz);
    wti1j = (1-dxm1/dx)*(dzm1/dz);
    wti1j1 = dxm1*dzm1/dx/dz;
    % update properties
    %ij
    ETASUM(i,j) = ETASUM(i,j) +etam(mat)*wtij;
    WTSUM(i,j) = WTSUM(i,j)+wtij;
    %i,j+1
    ETASUM(i,j+1) = ETASUM(i,j+1) +etam(mat)*wtij1;
    WTSUM(i,j+1) = WTSUM(i,j+1)+wtij1;
    %i+1,j
    ETASUM(i+1,j) = ETASUM(i+1,j) +etam(mat)*wti1j;
    WTSUM(i+1,j) = WTSUM(i+1,j)+wti1j;
    %i+1,j
    ETASUM(i+1,j+1) = ETASUM(i+1,j+1) +etam(mat)*wti1j1;
    WTSUM(i+1,j+1) = WTSUM(i+1,j+1)+wti1j1;
    
    %now interpolate density onto the vz grod (for the z stokes equation)
    j=fix((xm(m)-xvz(1))/dx)+1;
    i=fix((zm(m)-zvz(1))/dz)+1;
        if(j<1)
            j=1;
        elseif(j>nx)
            j=nx;
        end
        if(i<1)
            i=1;
        elseif(i>nz-1)
            i=nz-1;
        end
        
        %compute distances
        dxm1 = abs(xm(m)-xvz(j));
        dzm1 = abs(zm(m)-zvz(i));
        %compute weights
        wtij = (1-dxm1/dx)*(1-dzm1/dz);
        wtij1 = (dxm1/dx)*(1-dzm1/dz);
        wti1j = (1-dxm1/dx)*(dzm1/dz);
        wti1j1 = dxm1*dzm1/dx/dz;
        %update properties for the vz grid
        %ij
        RHOSUM2(i,j) = RHOSUM2(i,j) +rhom(mat)*wtij;
        KSUM2(i,j) = KSUM2(i,j) + wtij*km(mat);
        WTSUM2(i,j) = WTSUM2(i,j)+wtij;
        %i,j+1
        RHOSUM2(i,j+1) = RHOSUM2(i,j+1) +rhom(mat)*wtij1;
        KSUM2(i,j+1) = KSUM2(i,j+1) + wtij1*km(mat);
        WTSUM2(i,j+1) = WTSUM2(i,j+1)+wtij1;
        %i+1,j
        RHOSUM2(i+1,j) = RHOSUM2(i+1,j) +rhom(mat)*wti1j;
        KSUM2(i+1,j) = KSUM2(i+1,j) + wti1j*km(mat);
        WTSUM2(i+1,j) = WTSUM2(i+1,j)+wti1j;
        %i+1,j
        RHOSUM2(i+1,j+1) = RHOSUM2(i+1,j+1) +rhom(mat)*wti1j1;
        KSUM2(i+1,j+1) = KSUM2(i+1,j+1) + wti1j1*km(mat);
        WTSUM2(i+1,j+1) = WTSUM2(i+1,j+1)+wti1j1;
        
        % interpolate density onto the vx grid (for the x stokes equation)
        j=fix((xm(m)-xvx(1))/dx)+1;
        i=fix((zm(m)-zvx(1))/dz)+1;
        if(j<1)
            j=1;
        elseif(j>nx-1)
            j=nx-1;
        end
        if(i<1)
            i=1;
        elseif(i>nz)
            i=nz;
        end
        
        %compute distances
        dxm1 = abs(xm(m)-xvx(j));
        dzm1 = abs(zm(m)-zvx(i));
        %compute weights
        wtij = (1-dxm1/dx)*(1-dzm1/dz);
        wtij1 = (dxm1/dx)*(1-dzm1/dz);
        wti1j = (1-dxm1/dx)*(dzm1/dz);
        wti1j1 = dxm1*dzm1/dx/dz;
        %update properties for the vz grid
        %ij
        RHOSUM3(i,j) = RHOSUM3(i,j) +rhom(mat)*wtij;    
        KSUM3(i,j) = KSUM3(i,j) + wtij*km(mat);
        WTSUM3(i,j) = WTSUM3(i,j)+wtij;
        %i,j+1
        RHOSUM3(i,j+1) = RHOSUM3(i,j+1) +rhom(mat)*wtij1;
        KSUM3(i,j+1) = KSUM3(i,j+1) + wtij1*km(mat);
        WTSUM3(i,j+1) = WTSUM3(i,j+1)+wtij1;
        %i+1,j
        RHOSUM3(i+1,j) = RHOSUM3(i+1,j) +rhom(mat)*wti1j;
        KSUM3(i+1,j) = KSUM3(i+1,j) + wti1j*km(mat);
        WTSUM3(i+1,j) = WTSUM3(i+1,j)+wti1j;
        %i+1,j
        RHOSUM3(i+1,j+1) = RHOSUM3(i+1,j+1) +rhom(mat)*wti1j1;
        KSUM3(i+1,j+1) = KSUM3(i+1,j+1) + wti1j1*km(mat);
        WTSUM3(i+1,j+1) = WTSUM3(i+1,j+1)+wti1j1;

        %Compute on the P grid
     % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((zm(m)-zp(1))/dz)+1;
    if(j<1)
        j=1;
    elseif(j>nx)
        j=nx;
    end
    if(i<1)
        i=1;
    elseif(i>nz)
        i=nz;
    end
    % Compute distances
    dxm1=xm(m)-xp(j);
    dzm1=zm(m)-zp(i);
    %compute weights
    wtij = (1-dxm1/dx)*(1-dzm1/dz);
    wtij1 = (dxm1/dx)*(1-dzm1/dz);
    wti1j = (1-dxm1/dx)*(dzm1/dz);
    wti1j1 = dxm1*dzm1/dx/dz;
    %(i,j)
    WTSUMP(i,j) = WTSUMP(i,j) + wtij;
    WTSUMT(i,j) = WTSUMT(i,j) + wtij*rhom(mat)*cpm(mat);
    TSUM(i,j) =TSUM(i,j) + Tm(m)*wtij*rhom(mat)*cpm(mat);
    RHOPSUM(i,j) = RHOPSUM(i,j) + wtij*rhom(mat);
    CPSUM(i,j) = CPSUM(i,j) + wtij*cpm(mat);
    ETAPSUM(i,j) = ETAPSUM(i,j) +wtij*etam(mat);
    %i,j+1
    WTSUMP(i,j+1) = WTSUMP(i,j+1)+wtij1;
    WTSUMT(i,j+1) = WTSUMT(i,j+1) + wtij1*rhom(mat)*cpm(mat);
    TSUM(i,j+1) =TSUM(i,j+1) + Tm(m)*wtij1*rhom(mat)*cpm(mat);
    RHOPSUM(i,j+1) = RHOPSUM(i,j+1) + wtij1*rhom(mat);
    CPSUM(i,j+1) = CPSUM(i,j+1) + wtij1*cpm(mat);
    ETAPSUM(i,j+1) = ETAPSUM(i,j+1) +wtij1*etam(mat);
    %i+1,j
    WTSUMP(i+1,j) = WTSUMP(i+1,j)+wti1j;
    WTSUMT(i+1,j) = WTSUMT(i+1,j) + wti1j*rhom(mat)*cpm(mat);
    TSUM(i+1,j) =TSUM(i+1,j) + Tm(m)*wti1j*rhom(mat)*cpm(mat);
    RHOPSUM(i+1,j) = RHOPSUM(i+1,j) + wti1j*rhom(mat);
    CPSUM(i+1,j) = CPSUM(i+1,j) + wti1j*cpm(mat);
    ETAPSUM(i+1,j) = ETAPSUM(i+1,j) +wti1j*etam(mat);
    %i+1,j
    WTSUMP(i+1,j+1) = WTSUMP(i+1,j+1)+wti1j1;
    WTSUMT(i+1,j+1) = WTSUMT(i+1,j+1) + wti1j1*rhom(mat)*cpm(mat);
    TSUM(i+1,j+1) = TSUM(i+1,j+1) + Tm(m)*wti1j1*rhom(mat)*cpm(mat);
    RHOPSUM(i+1,j+1) = RHOPSUM(i+1,j+1) + wti1j1*rhom(mat);
    CPSUM(i+1,j+1) = CPSUM(i+1,j+1) + wti1j1*cpm(mat);
    ETAPSUM(i+1,j+1) = ETAPSUM(i+1,j+1) +wti1j1*etam(mat);
end


%compute Eta and Rho on real grid
Eta_out = ETASUM./WTSUM;
% staggered grids
Rho_vz = RHOSUM2./WTSUM2;
Rho_vx = RHOSUM3./WTSUM3;
k_vz = KSUM2./WTSUM2;
k_vx = KSUM3./WTSUM3;
T1 = TSUM./WTSUMT;
Rho_mid = RHOPSUM./WTSUMP;
CP_mid = CPSUM./WTSUMP;
Eta_mid = ETAPSUM./WTSUMP;
end

function Number = numsetup(nz,nx)
%setup numbering system
Number = zeros(nz,nx);
num = 1;
for j=1:nx    
for i=1:nz
Number(i,j) = num;
num = num+1;
end
end
end

