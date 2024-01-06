clc;
clearvars;
syms x;
numIterations = 123; % Set the total number of iterations
for i = 1:numIterations
    % Your processing or computations here
    fprintf('Completed %d/%d', i, numIterations);
    pause(0.1); % Just for demonstration, replace with your actual code
    fprintf(repmat('\b', 1, numel(['Completed %d/%d', i, numIterations])));
end
fprintf('Processing complete!\n');
functions = {2^sym(x)};
precision = 200;
start = 0;
rec_degree = 3;
poly_degree = 3;
solution = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, true); 
disp(solution)

function ret = substitute(n, func_val)
    ret = func_val(:,1) - 4/sym(n)*sym(n+1) * func_val(:,2) + 4/sym(n-1)*sym(n+1) * func_val(:,3);
end
