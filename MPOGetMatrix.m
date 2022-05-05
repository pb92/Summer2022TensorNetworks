function [returnMatrix] = MPOGetMatrix(inputMPO,dVector)

% Philip Daniel Blocher
% Center for Quantum Information and Control
% Department of Physics and Astronomy, University of New Mexico
% Email: blocher@unm.edu (current mail on Google Scholar)

% Description: This algorithm converts a matrix product operator (MPO) 
% back into a matrix containing the entries of the MPO.
% Last updated: November 2021.

%% Description of input and output:
% Input:
% inputMPO is the matrix product operator array.
% dVector is an (N times 1)-vector whose entries are the number of physical
% degrees of freedom for that particular site.

% Output:
% StateVector is a (n times 1)-vector, where each entry is a particular
% coefficient of the density matrix, and where n is the total number of
% entries in the density matrix.


%% Initialize working environment.
% Length of chain:
N = size(dVector,1);

% Count the Hilbert space size and initialize the density matrix storage:
dtot = 1;
for i = 1:1:N
    dtot = dtot * dVector(i);
end

returnMatrix = NaN(dtot,dtot);


%% Reconstruct density matrix from MPO decomposition.

% Initialize list of physical (unprimed) indices:
index = ones(size(dVector));

for i = 1:1:dtot

    % Initialize list of physical (primed) indices:
    index_prime = ones(size(dVector));
    
    for k = 1:1:dtot
    
        % Coefficient will ultimately be the density matrix coefficient. 
        % We multiply matrices according to index.
        coefficient = inputMPO{index(1),1,index_prime(1)};
        for j = 2:1:N
            coefficient = coefficient * inputMPO{index(j),j,index_prime(j)};
        end

        % Now take the trace of coefficient (which is only relevant when we
        % have periodic boundary conditions) and store the value.
        returnMatrix(i,k) = trace(coefficient);

        % Now update the index vector such that we calculate the next
        % coefficient.
        index_prime(1) = index_prime(1)+1;
        for j = 1:1:N
           if index_prime(j) > dVector(j) && k ~= dtot
              index_prime(j)=1;
              index_prime(j+1)=index_prime(j+1)+1;
           end
        end 
    end
    
    % Now update the index vector such that we calculate the next
    % coefficient.
    index(1) = index(1)+1;
    for j = 1:1:N
       if index(j) > dVector(j) && i ~= dtot
          index(j)=1;
          index(j+1)=index(j+1)+1;
       end
    end     
   
end

end