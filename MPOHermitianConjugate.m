function [outputMPO] = MPOHermitianConjugate(inputMPO,dVector)

% Philip Daniel Blocher
% Center for Quantum Information and Control
% Department of Physics and Astronomy, University of New Mexico
% Email: blocher@unm.edu (current mail on Google Scholar)

% Description: This algorithm returns the Hermitian conjugate of the input
% matrix product operator (MPO) inputMPO.
% Last updated: August 2021.


%% Initialization of environment
% Read out the number of sites:
N = size(dVector,1);

% Initialize the outputMPO to have the same overall structure as the
% inputMPO:
outputMPO = inputMPO;


%% Hermitian conjugate
for cIndex = 1:1:N
    di = dVector(cIndex);
    
    for sigma = 1:1:di
        for sigmaprime = 1:1:di
            outputMPO(sigma,cIndex,sigmaprime) = {conj(double(inputMPO{sigmaprime,cIndex,sigma}))};
        end
    end
    
end

end