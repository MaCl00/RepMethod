clc;
clearvars;

syms x;
functions = {factorial(sym(x)), 3^sym(x), (-1)^sym(x), (psi(sym(x)+1)- psi(sym(1)))};
precision = 30; 
start = 0;
rec_degree = 3;
poly_degree = 7;
solution = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true);

function ret = substitute(~, func_val)
    ret = func_val(:,1) - 4 * func_val(:,2) + 3 * func_val(:, 3);
end
