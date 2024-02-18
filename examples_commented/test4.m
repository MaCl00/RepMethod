clc;
clearvars;
% A testcase demonstrating the zero-check.
% 0 is always a solution to a homogeneous recursive relation, thus, we try
% to avoid this solution. However, it can happen, that a linear combination
% in the guesses is zero. If a solution is found, a check is run, to
% confirm, that the solution is not zero. 
% In this case, we guess 2^x, 3^x and 6^x. Since the programm also guesses
% pairwise products, 2^x*3^x is guessed as well. Thus, we have a trivial 
% solution in our guesses.
% Having a zero in the guesses slows the calculation significantly.
% The rest of the testcase is identical to test3
syms x;
functions = {gamma(sym(x)+1), 2^sym(x), 3^sym(x), 6^sym(x)};
precision = 300;
start = 0;
rec_degree = 3;
poly_degree = 5;
[func_name, solutions] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, "verbose", "quitSearch", "addMul"); 
% solutions has for each solution a vector with the coefficients
% func_name has the corresponding functions
disp(solutions)
disp(func_name)

function ret = substitute(n, func_val)
    ret = func_val(:,1) - 4/sym(n)*sym(n+1) * func_val(:,2) + 4/sym(n-1)*sym(n+1) * func_val(:,3);
end

