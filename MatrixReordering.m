%
%  Fast and Accurate Pseudoinverse for Real-world Sparse Matrices
%
%  This software may be used only for research evaluation purposes.
%  For other purposes (e.g., commercial), please contact the authors.
%

function [Ak, row, col] = MatrixReordering(AOrig, k)

[row.n, col.n] = size(AOrig);
row.nk = ceil(row.n * k);
col.nk = ceil(col.n * k);

nIters = 0;

row.gccsize = row.n;
col.gccsize = col.n;

row.gccind = 1:row.n;
col.gccind = 1:col.n;

row.totalind = 1:row.n;
col.totalind = 1:col.n;

row.cur_lpos = 1;
row.cur_rpos = row.n;

col.cur_lpos = 1;
col.cur_rpos = col.n;

tmpA = AOrig | AOrig;

i = 1;

while nIters == 0 || row.gccsize > row.nk || col.gccsize > col.nk
	nIters = nIters + 1;
	%fprintf('%d-th iteration...\n', nIters);
	A = tmpA(row.gccind, col.gccind);

	% compute row/col spokes/gcc/hubs
	[row, col] = RemoveHighDegree(A, row, col, k);

	% update ordering
    if row.gccsize > row.nk
        [row] = UpdateOrdering(row);
    end
    
    if col.gccsize > col.nk
        [col] = UpdateOrdering(col);
    end

	row.gccsize = length(row.gccind);
	col.gccsize = length(col.gccind);
    
%     clf;
%     hold on;
%     box on;
%     ax = gca;
%     ax.LineWidth = 1.3;
%     pbaspect([1 1 1]);
%     daspect([0.3 1.3 1]);
%     set(gcf,'color','w');
%     %spy(AOrig(fliplr(row.newind), fliplr(col.newind)));
%     set(gca,'FontSize', 30);
%     ylabel('Instance Nodes');
%     xlabel('Feature Nodes');
%     hold off;
% %     figure_path = sprintf('exp/reorder/FIG/REORDERING_ITER_%d.pdf', i);
% %     export_fig(figure_path, '-png', '-transparent', '-m2.5');
% %     i = i + 1;
%     if sf == true
%         figure_path = sprintf('exp/reorder/FIG/REORDERING_ITER_%d.png', i);
%         export_fig(figure_path, '-png', '-transparent', '-m2.5');
%         i = i + 1;
%     end
    
end

row.newind = row.newind(:);
col.newind = col.newind(:);
row.newind = flipud(row.newind);
col.newind = flipud(col.newind);

row.n2 = row.cur_lpos - 1;
col.n2 = col.cur_lpos - 1;
row.n1 = row.n - row.n2;
col.n1 = col.n - col.n2;

Ak = AOrig(row.newind, col.newind);

row.invind = zeros(row.n, 1);
row.invind(row.newind) = (1:row.n)';

col.invind = zeros(col.n, 1);
col.invind(col.newind) = (1:col.n)';

end

function [row] = UpdateOrdering(row)

row.totalind(row.cur_lpos:row.cur_lpos + row.topind_size - 1) = row.gccind(row.topind);
row.cur_lpos = row.cur_lpos + row.topind_size;
row.totalind(row.cur_rpos - row.disind_size + 1:row.cur_rpos) = row.gccind(row.disind);
row.cur_rpos = row.cur_rpos - row.disind_size;
row.newind = row.totalind;
row.gccind = row.gccind(row.newgccind);
row.newind(row.cur_lpos:row.cur_rpos) = row.gccind;

end

function [row, col] = RemoveHighDegree(A, row, col, k)

[nRows, nCols] = size(A);

% select hub nodes
row.topk = GetTopkNodes(A, row, 'row');
col.topk = GetTopkNodes(A, col, 'col');

[cclabels] = ComputeConnComp(A, A', row.topk, col.topk);

row.cclabels = cclabels(1:nRows);
col.cclabels = cclabels(nRows+1:end);

[row.C, row.ia, row.ic] = unique(row.cclabels);
row.H = accumarray(row.ic, 1); % # of row nodes for each CC
row.counts = [row.C, row.H];

[col.C, col.ia, col.ic] = unique(col.cclabels);
col.H = accumarray(col.ic, 1); % # of col nodes for each CC
col.counts = [col.C, col.H];

[row] = ComputeIndices(row);
[col] = ComputeIndices(col);

end

function [row] = ComputeIndices(row)

[~, row_gcc_idx] = max(row.H);
row.gcc = row.C(row_gcc_idx);
row.topind = row.topk;
row.newgccind = find(row.cclabels==row.gcc);
row.disind = find(row.cclabels~=row.gcc); 
row.disind = setdiff(row.disind, row.topind, 'stable');
[~, I] = sort(row.cclabels(row.disind));
row.disind = row.disind(I);
row.topind_size = length(row.topind);
row.disind_size = length(row.disind);

end

function [topks] = GetTopkNodes(A, data, type)
direction = 1;
if strcmp(type, 'row')
	direction = 2;
end
[~, I] = sort(sum(A, direction), 'descend');

if length(I) >= data.nk
    topks = I(1:data.nk);
else
    topks = I;
end

topks = topks(:);
end
