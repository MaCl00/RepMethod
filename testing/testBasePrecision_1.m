clc;
clearvars;
syms x;
functions = {};
starting = 758;
ending = 820;
start = 0;
rec_degree = 3;
poly_degree = 15;
precision = 3;
timestamp = datestr(now, 'yyyy_mm_dd-HH_MM_SS');
[~, solution] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true);
fid = fopen(['timings_', timestamp, '.txt'], 'w');
for i=starting:ending
    precision = i; 
    disp(size(solution));
    disp(i);
    tic
    [function_name, solution] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true);
    elapsedTime = toc;
    fprintf(fid, '%d\t%f\n', precision, elapsedTime);
end
fclose(fid); % Close the file
if 1==2
poly_degree = 10;
timestamp = datestr(now, 'yyyy_mm_dd-HH_MM_SS');
[~, ~] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true);
fid = fopen(['timings_', timestamp, '.txt'], 'w');
for i=starting:ending
    precision = i; 
    disp(size(solution));
    disp(i);
    tic
    [function_name, solution] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true);
    elapsedTime = toc;
    fprintf(fid, '%d\t%f\n', precision, elapsedTime);
end
fclose(fid); % Close the file
precision = 10;
poly_degree = 15;
timestamp = datestr(now, 'yyyy_mm_dd-HH_MM_SS');
[a, b] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true);
disp(a);
disp(b);
fid = fopen(['timings_', timestamp, '.txt'], 'w');
for i=starting:ending
    precision = i; 
    disp(size(solution));
    disp(i);
    tic
    [function_name, solution] = repertoire(functions, @substitute, 0, precision, start, rec_degree, poly_degree, false, true);
    elapsedTime = toc;
    fprintf(fid, '%d\t%f\n', precision, elapsedTime);
end
fclose(fid); % Close the file
end
function ret = substitute(n, func_val)
    ret = sym(n)*(sym(n)-1)*func_val(:,1) - 2*(sym(n)+1)*(sym(n)-1)* func_val(:,2)+(sym(n)+1)*sym(n)*func_val(:, 3);
end
