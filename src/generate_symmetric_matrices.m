function A = generate_symmetric_matrices(n)
% GENERATE_SYMMETRIC_MATRICES Generate a list of diagonal matrices with random entries.
%
%   f_vec = GENERATE_SYMMETRIC_MATRICES(n) creates a cell array of n 
%   symmetric dense matrices. The i-th matrix in the list will be of size i x i 
%   with random entries uniformly disitributed between 0 and 1.
%
%   Input:
%       n       - Number of matrices to generate. 
%                 An integer where 1 <= n.
%
%   Output:
%       f_vec   - A 1 x n cell array where the i-th cell contains a symmetric 
%                 matrix of size i x i with random entries.
%
%   Example:
%       matrices = GENERATE_SYMMETRIC_MATRICES(3);
%       This will generate a cell array containing 3 matrices:
%       1x1, 2x2, and 3x3, each with random entries.
%   Author: Nathan Rousselot

    A = randn(n,n);
    A = A + A';
end