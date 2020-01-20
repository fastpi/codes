%
%  Fast and Accurate Pseudoinverse for Real-world Sparse Matrices
%
%  This software may be used only for research evaluation purposes.
%  For other purposes (e.g., commercial), please contact the authors.
%

function mexCompile()

compile_str = {
   gen_command('.', 'ComputeConnComp');
};

for i=1:length(compile_str)
	fprintf('Evaluate:\n\t%s\n', compile_str{i});
	eval(compile_str{i});
	fprintf('Done!\n');
end

end

function [cmd] = gen_command(SOURCE_HOME_PATH, FILE_NAME)
    cmd = sprintf('mex(''CXXFLAGS=-ansi -std=c++0x -D_SNU_SOURCE -fPIC -pthread -O2 -DNDEBUG -g'', ''%s/%s.cpp'', ''-largeArrayDims'', ''-outdir'', ''%s'')', SOURCE_HOME_PATH, FILE_NAME, SOURCE_HOME_PATH);
end
