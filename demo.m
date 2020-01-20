%
%  Fast and Accurate Pseudoinverse for Real-world Sparse Matrices
%
%  This software may be used only for research evaluation purposes.
%  For other purposes (e.g., commercial), please contact the authors.
%

load('data.mat');
alpha = 0.1;

fprintf('Start psuedoinverse computation...\n');

tic;
[V, pinvS, UT, rank] = FastPI(A, alpha);
% pinvA = (V * pinvS) * UT
fastpi_time = toc;

fprintf('Psuedoinverse computation completed in %.4f sec...\n', fastpi_time);