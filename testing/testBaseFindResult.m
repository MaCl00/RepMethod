clc;
clearvars;
syms x;
functions = {};
start = 0;

for h=1:6
    max_poly = 50;
    if h==2
        functions = {factorial(sym(x)+1)};
        max_poly = 25;
    elseif h==3
        functions = {log(sym(x)+1)};
        max_poly = 25;
    elseif h==4
        functions = {log(sym(x)+1), factorial(sym(x)+1)};
        max_poly = 20;
    elseif h==5
        functions = {2^sym(x),3^sym(x)};
        max_poly = 20;
    elseif h==6
        functions = {2^sym(x),3^sym(x),5^sym(x)};
        max_poly = 10;
    end
    for j=1:3
        homogeneous = 0;
        num_s = 1;
        rec_degree = 2;
        if j==1
            rec_eq = @substitute;
            num_s = 2;
            rec_degree = 3;
        elseif j==2
            rec_eq = @substitute2;
        else
            rec_eq = @substitute3;
            homogeneous = -(sym(x)-1)^4+4*(sym(x)-1)^2+5*sym(x)-3;
        end
        precision = 9;
        timestamp = datestr(now, 'yyyy_mm_dd-HH_MM_SS');
        fid = fopen(['timings_', timestamp, '.txt'], 'w');
        for k=4:max_poly
            poly_degree = k;
            solved = false;
            function_name = 0;
            while ~solved
                precision = precision + 1;
                disp([h, j,k,precision]);
                [function_name, solution] = repertoire(functions, rec_eq, homogeneous, precision, start, rec_degree, poly_degree, false, false);
                solved = size(solution, 2) == num_s;
                if size(solution, 1) == 0
                    solved = false;
                end
            end
            exp = 0;
            for i=1:size(function_name, 2)
                exp = exp + function_name{i};
            end
            fprintf(fid, '%d\t%d\t'+string(exp)+'\n', poly_degree, precision);
            precision = precision - 1;
        end
        fclose(fid);
    end
end

function ret = substitute(n, func_val)
    ret = sym(n)*(sym(n)-1)*func_val(:,1) - 2*(sym(n)+1)*(sym(n)-1)* func_val(:,2)+(sym(n)+1)*sym(n)*func_val(:, 3);
end
function ret = substitute2(n, func_val)
    ret = (sym(n)-2)*func_val(:,1) - (sym(n)-1)*func_val(:,2);
end
function ret = substitute3(n, func_val)
    ret = func_val(:,1) - (sym(n)-1)*func_val(:,2);
end
