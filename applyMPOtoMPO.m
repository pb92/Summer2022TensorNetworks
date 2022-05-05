function [outputMPO] = applyMPOtoMPO(inputMPO1,inputMPO2,dVector)

% Philip Daniel Blocher
% Center for Quantum Information and Control
% Department of Physics and Astronomy, University of New Mexico
% Email: blocher@unm.edu (current mail on Google Scholar)

% Description: This algorithm applies a matrix product operator (MPO) to
% another MPO (MPO1 is applied to MPO2 from the left). This yields a new
% MPO that is then returned.
% Last updated: August 2021.

%% Description of inputs and output:
% Input:
% inputMPO1, inputMPO2 are (dmax times N times dmax)-arrays, each entry
% being a matrix.
% N is the number of sites in the 1D chain, while dmax is the largest
% number of degrees of freedom for a site in the 1D chain.
% dVector is a (N times 1)-vector containing the degrees of freedom for
% each site.

% Output:
% outputMPO outputs the resulting MPO from applying MPO1 to MPO2 from the
% left. The resulting array is again an MPO.

%% Initialization of environment:

% Define variables:
dmax = max(dVector);
N = size(dVector, 1);

% Initialize output storage:
outputMPO = cell(dmax,N,dmax);

%% Algorithm:

% For each site in the chainÆ
for k = 1:1:N

    % Read out physical degrees of freedomÆ
    dk = dVector(k);

    % For-loop over sigma:
    for sigma = 1:1:dk
        
        % For-loop over sigma'
        for sigmaprime = 1:1:dk
        
            Mtilde = 0;

            % Calculate the new MPO matrix:
            for sigmaprimeprime = 1:1:dk
                W = cell2mat(inputMPO1(sigma,k,sigmaprimeprime));
                M = cell2mat(inputMPO2(sigmaprimeprime,k,sigmaprime));
                Mtilde = Mtilde + kron(W,M);
            end
            
            % Store the new MPO matrix:
            outputMPO(sigma,k,sigmaprime) = {Mtilde};
            
        end
    end
end

end