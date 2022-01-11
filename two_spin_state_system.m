

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

omega_ONE=2*pi*10000; %Zeeman Frequency spin 1 
omega_TWO=2*pi*10000; %Zemman Frequency spin 2

H_Z=omega_ONE*(kron(sigma_z,unit)) + omega_TWO*(kron(unit,sigma_z)) %Zeeman coupling two spin states

J=15

H_J=2*pi*J*(kron(sigma_x,sigma_x)+kron(sigma_y,sigma_y) + kron(sigma_z,sigma_z)) %J coupling two spin states


H = H_Z + H_J %Total Hamiltonian

%initial state
rho=[0 1 0 0;0 0 0 0;0 0 0 1; 0 0 0 0;]


%Build Propogators
time_step=0.125/norm(H);
P_left=expm(+1i*H*time_step);
P_right=expm(-1i*H*time_step);

%Total operators of two spin state system
sigma_x_total=kron(sigma_x,unit) + kron(unit,sigma_x)
sigma_y_total=kron(sigma_y,unit) + kron(unit,sigma_y)
sigma_z_total=kron(sigma_z,unit) + kron(unit,sigma_z)


%Simulation
nsteps=1024;
for n=1:nsteps
    
    mu_x(n)=real(trace(rho*sigma_x_total));
    mu_y(n)=real(trace(rho*sigma_y_total));
    mu_z(n)=real(trace(rho*sigma_z_total));
    rho=P_right*rho*P_left;
   
    
end

time=[1:nsteps]*time_step
 plot(time,mu_x); grid on; box on;
% xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);


end
