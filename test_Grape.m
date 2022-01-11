function test_Grape()
 % Magnetic field
    sys.magnet=9.4;

    % Spin system
    sys.isotopes={'15N','1H','13C','13C','13C','15N'};

    % Chemical shifts, ppm
    inter.zeeman.scalar={0.0, 0.0, 55.0, 30.0, 175.0, 0.0};

    % Scalar couplings, Hz (literature values)
    inter.coupling.scalar=cell(6);
    inter.coupling.scalar{1,3}=-11; 
    inter.coupling.scalar{2,3}=140;
    inter.coupling.scalar{3,4}=35;
    inter.coupling.scalar{3,5}=55;
    inter.coupling.scalar{3,6}=7;
    inter.coupling.scalar{5,6}=-15;

    % Basis set
    bas.formalism='sphten-liouv';
    bas.approximation='IK-0';
    bas.level=4;

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    
    
      % Set up and normalise the initial state
    rho_init=state(spin_system,{'Lz'},{2});
    rho_init=rho_init/norm(full(rho_init),2);
    % Set up and normalise the target state
    rho_targ=state(spin_system,{'Lz'},{5});
    rho_targ=rho_targ/norm(full(rho_targ),2);
    
    
    
    
    % Drift Hamiltonian
    H=hamiltonian(assume(spin_system,'nmr'));
    
    % Control operators
    LpH=operator(spin_system,'L+','1H');
    LpC=operator(spin_system,'L+','13C');
    LpF=operator(spin_system,'L+','15N');
    LxH=(LpH+LpH')/2; LyH=(LpH-LpH')/2i;
    LxC=(LpC+LpC')/2; LyC=(LpC-LpC')/2i;
    LxF=(LpF+LpF')/2; LyF=(LpF-LpF')/2i;
    
    
    % Define control parameters
    control.drifts={{H}};                           % Drift
    control.operators={LxH,LyH,LxC,LyC,LxF,LyF};    % Controls
    control.rho_init={rho_init};                    % Starting state
    control.rho_targ={rho_targ};                    % Destination state
    control.pwr_levels=2*pi*1e3;                    % Pulse power
    control.pulse_dur=10e-3;                        % Pulse duration
    control.pulse_nsteps=50;                        % Time points
    control.penalties={'SNS'};                      % Penalty
    control.p_weights=100;                          % Penalty weight
    control.method='lbfgs';                         % Optimisation method
    control.freeze=zeros(6,50);                     % Freeze mask
    control.phase_cycle=[0  0  0  0*pi/2  -0*pi/2;
                         0  0  0  1*pi/2  -1*pi/2;
                         0  0  0  2*pi/2  -2*pi/2;
                         0  0  0  3*pi/2  -3*pi/2]; % Phase cycle
    
    % Plots during optimisation
    control.plotting={'correlation_order','local_each_spin','xy_controls'};
    
    
    
      % Spinach housekeeping
    spin_system=optimcon(spin_system,control);
    
    
    % Take a random guess for the waveform
    guess=randn(numel(control.operators),control.pulse_nsteps)/3;
    
    
        % Run the optimization
    fminnewton(spin_system,@grape_xy,guess);
end