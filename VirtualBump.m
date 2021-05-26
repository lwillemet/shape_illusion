% This script enables to calculate the displacements and interfacial pressure generated when the
% finger slide on a surface with a virtual bump (created with a varying friction
% coefficient)
%%
clear all
close all

global kext kt Kong kn B Ndx posi nu

disp('Initializing')
%% --------------------Parametrization--------------------
% % Skin parameters
E_ext = 1.54e6;         %% young's modulus (Pa) of the epidermis
d = 5e-3;               %% diameter of the contact area
e_ext = 20e-6;          %% thickness of the stratum corneum
e_epi = 0.15e-3;        %% thickness of the epidermis
E_tis = 0.025e6;        %% young's modulus (Pa) of the skin tissue
nu = 0.48;              %% poisson's ratio (Fung)

% % Skin structure
r = 8e-3;               %% finger radius in m
Ndx = 251;              %% Number of elements along the membrane

% % External load
P = 0.5;                %% initial pression in N on the bone
Q = 1;                  %% initial traction in N on the bone

%% ----------------Shape function of the finger----------------
d0 = 4e-4;                  %% initial distance from the surface
S = @(x) r-sqrt(r^2-x.^2);  %% shape function of the finger surface
dx = 1e-6;                  %% spatial discretization
x = -r+1e-6:dx:r-1e-6;
xtest = -r+1e-6:2*(r-1e-6)/(Ndx-1):r-1e-6;
xtest = xtest';
p = [x',S(x)'];
[q,dL] = curvspace(p,Ndx);  %% allows to have an equally spaced discretization along the curve
x_eq = q(:,1);
x_eq(floor(Ndx/2)+1) = 0;
xi = x_eq;                  %% abscisses of the equally spaced points

bone_i = [r;0];             %% initial position of the bone

%% ------------------Matrix of dependency------------------
    %%Stiffness matrix (spring)
kext = E_ext*d*e_ext/dL;        %% stiffness of the spring = membrane
kt = E_tis*pi*(r)^2/2/(r);      %% stiffness of the spring = tissue
kn = kt;                        %% stiffness of the spring = nail

% % Matrix of the spring dependencies relevant to the nail
idx_i =  [1, 2,Ndx+1,1,2, 1,    Ndx+1];
idx_j =  [1, 1,1,    2,2, Ndx+1,Ndx+1];
values = [-2,1,1,    1,-1,1,    -1];
Kong = sparse(idx_i,idx_j,values,2*(Ndx+1),2*(Ndx+1));

%% -------------------------Contact-------------------------
kc = 100*kt;                                %% stiffness of the contact spring

period = 4e-3;
posi = [r;S(xi);-2*period;xi-2*period];     %% initial position of the elements

% % Virtual bump
AVir = 0.4;                                 %% friction coefficient difference
Amean = 0.6;
Amin = Amean-AVir/2;
fclbVir = @(x) (AVir/2*(sin(2*pi/period*x+pi)+1)+Amin).*(x>=-period/2 & x<=period/2)+...
    Amean.*(x<-period/2) + Amean.*(x>period/2); %% variation of the static friction coefficient

    % % Dahl friction model
sigma0 = 1e3;                 %% rest stiffness of the bristle (N/m)
n = 0.7;                      %% exponent: material dependent parameter
                               	% 0<=n<=1 for brittle materials
                                
    % % Damping matrix (dashpot)
zeta = 0.1;            
B = -zeta*speye(2*(Ndx+1));

% % Temporal parameters
Fsmin = 1/(2*zeta/kc);    %% b>KT/2 (Colgate)
Fs = 2e5;               %% sampling frequency
dt = 1/Fs;
time = 0:dt:0.5;          %% time vector
Ndt = length(time);


%% -------------------------RK4 function:-------------------------
    % Uc is the current displacement vector
    % Fc the current force vector
    % tc the current instant
    % s is the sign of the current x
% % see at the bottom of the script for the function f_U

%% -----------------------Resolution-----------------------
% % Initial conditions
UVir = zeros(2*(Ndx+1),length(time));
posVir(1:2*(Ndx+1),1) = posi;
UptVir = zeros(Ndx+1,length(time));
UpnVir = zeros(Ndx+1,length(time));
FVir = zeros(2*(Ndx+1),length(time));
FVir(:,1:Ndt) = repmat([-P;zeros(Ndx,1);0;zeros(Ndx,1)],1,Ndt);
contact_areaVir = NaN(Ndx,Ndt);
muVir = NaN(Ndx,Ndt);

t=2;
% % Impedance loading
disp('Computing')
while(t<=2 || abs(sum(FVir(1:Ndx+1,t-1),1))>1/10)      %% Stop condition: global net force <10%
    % % Force resolution
    % % Contact modeling with a high stiffness spring -> force
    for i=1:Ndx+1
        % % Computation of the speeds
        if(t>2)
            UpnVir(i,t-1) = (UVir(i,t-1)-UVir(i,t-2))/dt;
            UptVir(i,t-1) = (UVir(i+Ndx+1,t-1)-UVir(i+Ndx+1,t-2))/dt;
        end
        % % if element i enters in contact
        if(posVir(i,t-1) <0)
            contact_areaVir(i-1,t) = posVir(i+Ndx+1,t-1);
            Fr = kc*(abs(posVir(i,t-1)));
            FVir(i,t) = FVir(i,t) + Fr;
            FVir(i+Ndx+1,t) = (FVir(i+Ndx+1,t) + FVir(i+Ndx+1,t-1) - dt*...
                    sigma0*UptVir(i,t-1)*(abs(1+FVir(i+Ndx+1,t-1)/(fclbVir(posVir(i+Ndx+1,t-1))*Fr)*sign(UptVir(i,t-1)))).^n*...
                    sign(1+FVir(i+Ndx+1,t-1)/(fclbVir(posVir(i+Ndx+1,t-1))*Fr)*sign(UptVir(i,t-1))));
            if(i == floor(Ndx/2)+2)
                FVir(i+Ndx+1,t) = 0;
            end
            muVir(i-1,t) = abs(FVir(i+Ndx+1,t)./FVir(i,t));
        else
            FVir(:,t) = FVir(:,t);
        end
    end
    
    % % Runge-Kutta4 : displacement resolution
    k1 = f_U(t*dt,UVir(:,t-1),FVir(:,t));
    k2 = f_U((t+0.5)*dt,UVir(:,t-1)+0.5*dt*k1,FVir(:,t));
    k3 = f_U((t+0.5)*dt,UVir(:,t-1)+0.5*dt*k2,FVir(:,t));
    k4 = f_U((t+1)*dt,UVir(:,t-1)+dt*k3,FVir(:,t));
    UVir(:,t) = UVir(:,t-1)+ dt/6*(k1+2*k2+2*k3+k4);
    
    % % Adding shape function -> final position of each element
    posVir(:,t) = UVir(:,t) + posi;
    
    t = t+1;
end
ti=t;
%%
% % We add traction Q on the bone
disp('Adding tangential stress...')
FVir(:,ti:Ndt) = repmat([-P;zeros(Ndx,1);Q;zeros(Ndx,1)],1,Ndt-ti+1);
UVir(:,ti:Ndt) = zeros(2*(Ndx+1),Ndt-ti+1);
    % % IDEM
for t=ti:Ndt
    for i=1:Ndx+1
        if(t>2)
            UpnVir(i,t-1) = (UVir(i,t-1)-UVir(i,t-2))/dt;
            UptVir(i,t-1) = (UVir(i+Ndx+1,t-1)-UVir(i+Ndx+1,t-2))/dt;
        end
        if(posVir(i,t-1) <0)
            contact_areaVir(i-1,t) = posVir(i+Ndx+1,t-1);
            Fr = kc*(abs(posVir(i,t-1)));
            FVir(i,t) = FVir(i,t) + Fr;
            FVir(i+Ndx+1,t) = (FVir(i+Ndx+1,t) + FVir(i+Ndx+1,t-1) - dt*...
                    sigma0*UptVir(i,t-1)*(abs(1+FVir(i+Ndx+1,t-1)/(fclbVir(posVir(i+Ndx+1,t-1))*Fr)*sign(UptVir(i,t-1)))).^n*...
                    sign(1+FVir(i+Ndx+1,t-1)/(fclbVir(posVir(i+Ndx+1,t-1))*Fr)*sign(UptVir(i,t-1))));
            muVir(i-1,t) = abs(FVir(i+Ndx+1,t)./FVir(i,t));
        else
            FVir(:,t) = FVir(:,t);
        end
    end
    
    % % Runge-Kutta4 : displacement resolution
    k1 = f_U(t*dt,UVir(:,t-1),FVir(:,t));
    k2 = f_U((t+0.5)*dt,UVir(:,t-1)+0.5*dt*k1,FVir(:,t));
    k3 = f_U((t+0.5)*dt,UVir(:,t-1)+0.5*dt*k2,FVir(:,t));
    k4 = f_U((t+1)*dt,UVir(:,t-1)+dt*k3,FVir(:,t));
    UVir(:,t) = UVir(:,t-1)+ dt/6*(k1+2*k2+2*k3+k4);
    
    % % Adding shape function -> final position of each element
    posVir(:,t) = UVir(:,t) + posi;
    
end
% delete(f)
%% -----------------------Plotting-----------------------
h = figure('units','normalized','outerposition',[0 0 1 1],'Color','w');
hold on
axis off
% filename = 'virtual_bump.gif';
% frame = getframe(h);
% im = frame2im(frame);
% [imind,cm] = rgb2ind(im,256);
% imwrite(imind,cm,filename,'gif','LoopCount',Inf,'DelayTime',5e-4);
ndraw = 2000;
vec_plot = floor((Ndx+1)/2)+1-50:2:floor((Ndx+1)/2)+1+50;
sz = 20;
for k=1:Ndt
    if (~mod(k-1,ndraw))
        clf;
        subplot 211
        axis off
        hold on
        area([-3/2*period 3/2*period],[0 0],-1e-4,'Facecolor','k')
        cd = [NaN;muVir(vec_plot(1:end-1),k)];
        cd(isnan(cd)) = zeros(length(cd(isnan(cd))),1);
        scatter(posVir(Ndx+1+vec_plot,k),posVir(vec_plot,k),sz,cd,'filled');
        caxis([0 1])
        
        quiver(posVir(Ndx+1+vec_plot,k),posVir(vec_plot,k),...
            FVir(Ndx+1+vec_plot,k),FVir(vec_plot,k),0.5,'Color','b')
        
        xlim([-3/2*period 3/2*period])
        xlabel('x (m)')
        axis equal
        ylabel('y (m)')
        
        subplot 212
        plot(posVir(Ndx+3:end,k),abs(FVir(Ndx+3:end,k)));
        xlabel('Position')
        ylabel('Q/fP')
        xlim([-3/2*period 3/2*period])
        ylim([0 0.05])
        
        drawnow
        % Capture the plot as an image
%         frame = getframe(h);
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
%         if k==1
%             imwrite(imind,cm,filename,'gif','LoopCount',Inf,'DelayTime',5e-3);
%         else
%             imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',5e-3);
%         end
        pause(1e-5);
    end
end

%% ----------------------Strains----------------------
vec_plot = floor((Ndx+1)/2)+1-50:2:floor((Ndx+1)/2)+1+50;
sz = 10;
position = [-4,-2,0,2,4].*1e-3;
for jj=1:length(position)
    [~,instants(jj)] = min(abs(posVir(Ndx+2,:) - position(jj)));
end
% % Plot
figure; hold on
for kk=1:length(instants)
    k = instants(kk);
    
    subplot(2,length(instants),kk)
    axis off
    hold on
    area([-2*period 2*period],[0 0],-1e-4,'Facecolor','k')
    plot(posVir(Ndx+1+vec_plot,k),posVir(vec_plot,k),'-b');
    cd = [NaN;muVir(vec_plot(1:end-1),k)];
    cd(isnan(cd)) = zeros(length(cd(isnan(cd))),1);
    scatter(posVir(Ndx+1+vec_plot,k),posVir(vec_plot,k),sz,cd,'filled');
    
    caxis([0 1])
    quiver(posVir(Ndx+1+vec_plot,k),posVir(vec_plot,k),...
        FVir(Ndx+1+vec_plot,k),FVir(vec_plot,k),0.5,'Color','b')
    xlim([-5e-3 5e-3])
    ylim([-1e-4 4e-3])
    
    subplot(2,length(instants),kk+length(instants))
    hold on
    % muP = area(pos(Ndx+1+vec_plot,k),-fclb.*F(vec_plot,k),-1e-4,'Facecolor',[0 0.5 0.5]);
    % alpha(muP,0.3);
    % plot(pos(Ndx+1+vec_plot(1:end),k),F(Ndx+1+vec_plot,k));
    plot(posVir(Ndx+1+vec_plot(1:end-1),k).*1e3,diff(UVir(Ndx+1+vec_plot,k))./diff(posVir(Ndx+1+vec_plot,k)));
    
    xlabel('Position (mm)')
    ylabel('Tangential strain')
    xlim([-5 5])
    ylim([-0.5 0.5])

end