%This is how we comment 

function hello_world()

% Define Pauli matrices
sigma_x= [0 1/2; 1/2 0];
sigma_y=[0 -1i/2; 1i/2 0];
sigma_z=[1/2 0; 0 -1/2];
unit=[1 0;0 1];

%Calculate two spin operators
% 
% Lx=kron(sigma_x,unit);
% Ly=kron(sigma_y,unit);
% Lz=kron(sigma_z,unit);
% 
% Sx=kron(unit,sigma_x)
% Sy=kron(unit,sigma_y)
% Sz=kron(unit,sigma_z)

%Sx*Sy-Sy*Sx-1i*Sz
%Lx*Sz-Sz*Lx
%Lz*Lx-Lx*Lz-1i*Ly

%(sigma_x*sigma_y-sigma_y*sigma_x)-1i*sigma_z

%LxSx=kron(sigma_x,sigma_x)
%Lx*Sx

%Build Hamiltonian

omega_ONE=2*pi*1000; %Zeeman Frequency spin 1 
omega_TWO=2*pi*1000; %Zemman Frequency spin 2

H_Z=omega_ONE*(kron(sigma_z,unit)) + omega_TWO*(kron(unit,sigma_z))
%{


H_J
%initial state
rho=[0 1;0 0];

%Build Propogators
time_step=0.125/norm(H);
P_left=expm(+1i*H*time_step);
P_right=expm(-1i*H*time_step);

%Simulation
nsteps=1024;
for n=1:nsteps
    
    mu_x(n)=real(trace(rho*sigma_x));
    mu_y(n)=real(trace(rho*sigma_y));
    mu_z(n)=real(trace(rho*sigma_z));
    rho=P_right*rho*P_left;
    
end

plot3(mu_x,mu_y,mu_z); grid on; box on;
% xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
%}

end
