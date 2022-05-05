function [MPO_Odd, MPO_Even] = CreateTransverseFieldIsingModelMPO(tau,L,J,B)

% Philip Daniel Blocher
% Center for Quantum Information and Control
% University of New Mexico, New Mexico, USA.
% blocher@unm.edu

% Purpose: Creates an odd-bond / even-bond MPO for the time evolution
% operator of the Transverse Field Ising Model with L sites. The parameter
% tau is the step size, while J (B) is the NN interaction (transverse
% field) strength, and L is the chain length.

% Last updated: March 2021.
Version = '20210407';
% Update log:
% 20210310: Creation data.
% 20210407: Fixed error in matrix representations of Hamiltonians.
% CreateHamiltonianMPO requires the 'usual' basis order.

%% Define the bulk and edge bond Hamiltonians
% hB denotes the bulk Hamiltonian, hLE the left edge Hamiltonian, and hRE
% the right edge Hamiltonian. In the TFIM model, hLE = hB. Note that we are
% using the basis {|uu>, |ud>, |du>, |dd>}, where |u> is spin up and |d> is
% spin down.

hLE = [-J,  0, -B,  0;
        0,  J,  0, -B;
       -B,  0,  J,  0;
        0, -B,  0, -J];

hB = hLE;

hRE = [-J, -B, -B,  0;
       -B,  J,  0, -B;
       -B,  0,  J, -B;
        0, -B, -B, -J];
    

if L == 2
    hLE = hRE;
end


%% Creation of the Odd and Even bond MPOs
% We use the CreateHamiltonianMPO.m algorithm for creating odd and even 
% bond MPOs. The degrees of freedom of each spin site is:
d = 2;

% Return the MPO_Odd and MPO_Even:
[MPO_Odd, MPO_Even] = CreateHamiltonianMPO(hB,hLE,hRE,tau,d,L);