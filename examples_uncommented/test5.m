clc;
clearvars;
syms x;
functions = {factorial(sym(x)), (2)^sym(x), (-1)^sym(x), psi(sym(x)+1)- psi(sym(1))};
precision = 100; 
start = 0;
rec_degree = -1;
poly_degree = 1;
[function_name, solution] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, true, true);
disp(solution)
function ret = substitute(~, func_val)
    ret = func_val(:,1) - sum(func_val(:,2:end), 2);
end
