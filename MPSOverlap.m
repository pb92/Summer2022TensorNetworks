function [overlapValue] = MPSOverlap(MPS1, MPS2, dVector)

% Philip Daniel Blocher
% Quantum Optics Group
% Department of Physics and Astronomy, Aarhus University
% pblocher@phys.au.dk
% Algorithm for calculating the overlap between two states in their MPS
% representation, <MPS1|MPS2>.
% Last updated: November 2017.

%% Description of input and output:
% Input:
% MPS1 is a matrix product state and will be used as a bra <MPS1|.
% MPS2 is likewise an MPS and will be used as a ket |MPS2>.
% dVector is a vector in which each entry is the number of degrees of
% freedom on each site.

% Output:
% overlapValue is the value of the inner product <MPS1|MPS2>, i.e. the
% overlap of the two MPS being considered.


%% Initialization of environment:

% Define variables:
L = size(dVector, 1);

% Define a storage for the calculation:
MatrixStorage = cell(1, L+1);
MatrixStorage(1) = {1};

% For each site:
for i = 1:1:L
    
    % Define storage variable for site:
    D = 0;
    di = dVector(i);
    
    for j = 1:1:di
        % Read out matrices:
        A = cell2mat(MPS1(j,i));
        B = cell2mat(MPS2(j,i));
        C = cell2mat(MatrixStorage(i));
        
        D = D + A'*C*B;
    end
    MatrixStorage(i+1) = {D};
end

overlapValue = cell2mat(MatrixStorage(end));
end