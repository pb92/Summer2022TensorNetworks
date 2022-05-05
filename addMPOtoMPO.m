function [outputMPO] = addMPOtoMPOOBC(inputMPO1,inputMPO2,dVector)

% Philip Daniel Blocher
% Center for Quantum Information and Control
% Department of Physics and Astronomy, University of New Mexico
% Email: blocher@unm.edu (current mail on Google Scholar)

% Description: This algorithm adds a matrix product operator (MPO) to
% another MPO (MPO1 is added to MPO2 from the left). This yields a new
% MPO that is then returned.
% Last updated: September 2021.

%% Description of inputs and output:
% Input:
% inputMPO1, inputMPO2 are (dmax times N times dmax)-arrays, each entry
% being a matrix.
% N is the number of sites in the 1D chain, while dmax is the largest
% number of degrees of freedom for a site in the 1D chain.
% dVector is a (N times 1)-vector containing the degrees of freedom for
% each site.

% Output:
% outputMPO outputs the resulting MPO from adding MPO1 to MPO2 from the
% left. The resulting array is again an MPO.

%% Initialization of environment:

% Define variables:
dmax = max(dVector);
L = size(dVector, 1);

% Initialize output storage:
outputMPO = cell(dmax,L,dmax);

%% Algorithm:

for n = 1:1:L
    dn = dVector(n);
    
    for sigma = 1:1:dn
        
        for sigmaprime = 1:1:dn
            
            % Read out the two matrices from MPO1 and MPO2 to be added
            % together:            
            M = cell2mat(inputMPO1(sigma,n,sigmaprime));
            W = cell2mat(inputMPO2(sigma,n,sigmaprime));
            
            % We consider the matrices M and W as diagonal blocks of the
            % resulting matrix outputMatrix = diag(M,N) (for the first site
            % we instead have outputMatrix = [M, N], and for the final site
            % we have outputMatrix = [M, N]^T).
            
            if n == 1
                % At the first site we need to form the matrix 
                % outputMatrix = [M, N]:
                
                outputMatrix = [M, W];
                
                outputMPO(sigma,n,sigmaprime) = {outputMatrix};
                
                
            elseif n == L
                % At the final site we need to create the matrix
                % outputMatrix = [M, N]^T = [M; N]:
                
                outputMatrix = [M; W];
                
                outputMPO(sigma,n,sigmaprime) = {outputMatrix};
                
                
            else
                % We are looking at bulk sites and thus need to form
                % outputMatrix = diag(M,N). First we read out the
                % dimensions of matrices M and W:
                
                [m1,m2] = size(M);
                [w1,w2] = size(W);
            
                % Then we create a larger matrix in which M and W will fit
                % on the diagonal:
                outputMatrix = zeros(m1+w1,m2+w2);

                % Finally we add the two matrices to the diagonal:
                outputMatrix(1:m1,1:m2) = M;
                outputMatrix(m1+1:m1+w1,m2+1:m2+w2) = W;
                
                % The matrix is stored in the outputArray:
                outputMPO(sigma,n,sigmaprime) = {outputMatrix};
                
            end
        end 
    end 
end

end