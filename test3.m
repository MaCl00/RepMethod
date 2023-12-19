clc;
clearvars;
% A small testcase for a recursive equation.
% The recursive equation is a_n = a_(n-1) * H_n / H_(n-1).
% H_n=psi(n+1)-psi(1) is the n-th term of the Harmonic sequence.
% The solution is obivously H_n itself.
% The guesses are x!, H_x, 1, x, ..., x^5 and pairwise products, such as 
% x^2*H_x.
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
solution = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true); 
disp(solution)

function ret = substitute(n, func_val)
    ret = func_val(:,1) - 4/sym(n)*sym(n+1) * func_val(:,2) + 4/sym(n-1)*sym(n+1) * func_val(:,3);
end
