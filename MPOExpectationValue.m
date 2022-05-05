function [expectationValue] = MPOExpectationValue(inputMPO, inputMPS, dVector)

% Philip Daniel Blocher
% Quantum Optics Group
% Department of Physics and Astronomy, Aarhus University
% pblocher@phys.au.dk
% Algorithm for calculating the expectation value of a matrix product
% operator (MPO) in a matrix product state (MPS) |MPS>.
% Last updated: November 2017.

%% Description of input and output:
% Input:
% MPS is a matrix product state.
% dVector is a vector in which each entry is the number of degrees of
% freedom on each site.

% Output:
% expectationValue is the result of <MPS| MPO |MPS>.

% This function utilizes the functions 'applyMPOtoMPS.m' and 'MPSOverlap.m'
% in order to do the needed calculations.


%% Algorithm:

[resultingMPS] = applyMPOtoMPS(inputMPO, inputMPS, dVector);

[expectationValue] = MPSOverlap(inputMPS, resultingMPS, dVector);


end