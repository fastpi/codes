%
%  Fast and Accurate Pseudoinverse for Real-world Sparse Matrices
%
%  This software may be used only for research evaluation purposes.
%  For other purposes (e.g., commercial), please contact the authors.
%

function [V, pinvS, UT, rank] = FastPI(A, alpha)
% Compute pseudoinverse of given matrix A with target rank ratio alpha
% This implementation is for overdetermined cases
% For underdetermined cases, make A transposed before using it
%
% Input
%   - A: (m x n) design matrix
%   - alpha: target rank ratio
%
% Output:
%   - SVD results for the pseudoinverse such that pinvA = (V * pinvS) * UT
%   - rank: target rank computed by ceil(n * alpha)

AOrig = A;
k = 0.01;              % hub selection ratio
n = size(A, 2);        % # of features
r = ceil(n * alpha);   % target rank

% step 1. reorder and partition
[A, row, col] = MatrixReordering(AOrig, k);
m1 = row.n1;
n1 = col.n1;

A11 = A(1:m1, 1:n1);
A12 = A(1:m1, n1+1:end);
A21 = A(m1+1:end, 1:n1);
A22 = A(m1+1:end, n1+1:end);

% step 2. compute the SVD of A11
[U, S, VT, s] = reordered_svd(A11, alpha);

% step 3. update the SVD for A21
P = [S * VT; A21];
[Ut, St, VtT] = lra(P, s);
Ut1 = Ut(1:s, :);
Ut2 = Ut(s+1:end, :);
U = [U * Ut1; Ut2];
S = St;
VT = VtT;

% step 4. update the SVD for [A12; A22]
T = [A12; A22];
P = [U * S, T];
[Ut, St, VtT] = lra(P, r);
VtT1 = VtT(:, 1:s);
VtT2 = VtT(:, s+1:end);
VT = [VtT1 * VT, VtT2];
U = Ut;
S = St;

U = U(row.invind, :);

r = size(S, 1);
vec = diag(S);
I = find(vec ~= 0);
vec(I) = 1 ./ vec(I);
pinvS = spdiags(vec(:), 0, r, r);

VT = VT(:, col.invind);
UT = U';
V = VT';
rank = r;

end


function [U, S, VT, s] = reordered_svd(A11, alpha)

rowsum = full(sum(A11, 2));
row_zeroind = find(rowsum == 0);
row_non_zeroind = find(rowsum > 0);
row_newind = [row_non_zeroind; row_zeroind];

colsum = full(sum(A11, 1));
col_zeroind = find(colsum == 0);
col_non_zeroind = find(colsum > 0);
col_newind = [col_non_zeroind(:); col_zeroind(:)];

row_invind(row_newind) = 1:length(row_newind);
col_invind(col_newind) = 1:length(col_newind);

row_nzn = length(row_non_zeroind);
col_nzn = length(col_non_zeroind);
row_zn  = length(row_zeroind);
col_zn  = length(col_zeroind);

rA11 = A11(row_newind, col_newind);

T = rA11(1:row_nzn, 1:col_nzn);

s = ceil(size(A11, 2) * alpha);

if s >= col_nzn
    st = size(T, 2);
    [U, S, VT] = lra(T, st);
    
    augU = [U, sparse(row_nzn, col_zn); sparse(row_zn, st), speye(row_zn, col_zn)];
    augS = [S, sparse(st, row_zn); sparse(col_zn, st), sparse(col_zn, row_zn)];
    augVT = [VT, sparse(st, col_zn); sparse(row_zn, col_nzn), speye(row_zn, col_zn)];
    U = augU(:, 1:s);
    S = augS(1:s, 1:s);
    VT = augVT(1:s, :);
else
    s = ceil(size(T, 2) * alpha);
    [U, S, VT] = lra(T, s);
    
    ts = size(S, 1);
    U = [U; sparse(row_zn, ts)];
    VT = [VT, sparse(ts, col_zn)];
end

U = U(row_invind, :);
VT = VT(:, col_invind);
S = sparse(S);


end