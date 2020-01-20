%
%  Fast and Accurate Pseudoinverse for Real-world Sparse Matrices
%
%  This software may be used only for research evaluation purposes.
%  For other purposes (e.g., commercial), please contact the authors.
%

function [U, S, VT] = lra(A, r)
% Low-rank approximation
%
% Input
%   - A: input matrix
%   - r: target rank
%
% Ouput
%   - [U, S, VT]: low-rank results based on SVD such that A \simeq U * S * VT 

[m, n] = size(A);

if r > n
    r = n;
end

[U, S, V] = svd(full(A), 'econ');
VT = V';

% low-rank approximation
U = U(:, 1:r);
VT = VT(1:r, :);
S = S(1:r, 1:r);
S = sparse(S);

end
