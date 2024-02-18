clc;
clearvars;
syms x;
functions = {psi(sym(x)+1)- psi(sym(1))};
precision = 200; 
start = 2;
rec_degree = -1;
poly_degree = 30;
[function_name, solution] = repertoire(functions, @substitute, sym(x)+1, precision, start, rec_degree, poly_degree, "verbose", "quitSearch");
disp(function_name)
disp(solution)
function ret = substitute(n, func_val)
    ret = func_val(:,1) - 2/sym(n)*sum(func_val(:,2:end), 2);
end
