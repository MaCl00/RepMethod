function solution = repertoire(func_guess, recursive_relation, nonHomogeneous, precision, start_point, rec_degree, num_poly, verbose)
    digits(precision);
    syms x;
    
    num_func = length(func_guess);
    % The number of required data points
    % 1 + 2p + p*f where p is the number of polynomials and f is the number
    % of functions
    num_substitution = 1 + num_poly * ((num_func + num_func^2) / 2 + 1);
    disp(num_substitution);
    num_function_values = num_substitution + rec_degree;
    function_values = zeros(num_substitution-1,num_function_values);
    % ---------------------Calculating function-values---------------------
    range = num2cell(start_point:start_point+num_function_values-1);
    % A matrix where enties are the input functions evaluated at the points
    fprintf('Evaluating function-guesses');
    sequences = cellfun(@(j) cellfun(@(i) subs(func_guess{j}, x, sym(i)), range), num2cell(1:num_func), 'UniformOutput', false);
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    disp('Finished evaluating function-guesses');
    fprintf('Calculating product-functions');
    function_values(1:num_func, :) = vertcat(sequences{:});
    current_row = num_func + 1;
    N = length(sequences);
    % Adding f*g if f and g are in the input
    counter = 1;
    prod_func = cell(1, (N*N-N)/2);
    digits(precision);
    c1 = zeros(1, num_func);
    for i=1:num_func
        c1(i)= subs(func_guess{i}, x, num_function_values);
    end
    disp(c1(5));
    digits(precision * 100);
    disp(c1(5));
    c2 = zeros(1, num_func);
   
    for i=1:num_func
        c2(i)= subs(func_guess{i}, x, num_function_values);
    end
    disp(c2(5));
    c1 = c2 - c1;
    for i = 1:N
        for j = i+1:N
            function_values(current_row, :) = function_values(i, :) .* function_values(j, :);
            current_row = current_row + 1;
            prod_func{counter} = func_guess{i} * func_guess{j};
            counter = counter + 1;
        end
    end
    func_guess = [func_guess, prod_func];
    % fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    for i = 1:num_func
        fprintf('%.2e', c1(i));
        disp("\n");
    end
    disp('Finished calculating product-functions');
    fprintf("Calculating polynomes");
    % Calculating polynomes up to the degree "rec_degree"
    polynomes = cell(1, num_poly);
    poly_name = cell(1, num_poly);
    if num_poly >= 1
        polynomes{1} = cellfun(@(i) 1, range);
        poly_name{1} = prod(sym(x)+1:sym(x));
    end
    for j=1:num_poly-1
        fprintf('%5d of %5d)\n', j+1, num_poly);
        polynomes{j+1} = cellfun(@(i) polynomes{j}(i- start_point + 1) * (sym(i) - j + 1), range);
        poly_name{j+1} = prod(sym(x)-j+1:sym(x));
        fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    end
    polyMatrix = vertcat(polynomes{:});
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    disp('Finished calculating polynomes');
    fprintf("Calculating function-polynom products");
    % Adding f*p for f is in the input and p is a polynom previous
    % calculated
    counter = 1;
    function_name = cell(1, (current_row-1)*(num_poly-1));
    for i = 1:(current_row-1)
        for j= 2:num_poly
            function_values(current_row, :) = function_values(i, :) .* polyMatrix(j, :);
            current_row = current_row + 1;
            function_name{counter} = func_guess{i} *poly_name{j};
            counter = counter + 1;
        end
    end
    function_values(current_row:current_row + num_poly - 1, :) = polyMatrix;
    current_row = current_row + num_poly;
    function_name = [func_guess, function_name, poly_name];
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    disp('Finished calculating polynomes');
    disp("Finished calulating function-values");
    %--------------------Checking best trivial solution--------------------
    % Will probably change, is not relavant right now
    fprintf("Calculating best null in input");
    
    A = function_values';
    [~, ~, V2] = svd(A);
    x2 = V2(:, end);
    y2 = A * x2;
    n_best = norm(y2, 'fro');
    % solution = 0;
    % disp(n_best)
    % counter = 1;
    % while counter <= length(function_name)
    %     if abs(norm(x2(counter)*sequences(counter), 'fro')) <= 6
    %         A(:, counter) = [];            % Delete i-th column of A
    %         x2(counter) = [];              % Delete i-th entry in x2
    %         function_name(counter) = [];   % Delete i-th entry in function_name
    %     else
    %         counter = counter + 1;
    %     end
    % end
    % disp(length(function_name))
    % [~, ~, V2] = svd(A);
    % x2_original = V2(:, end);
    % n_best_original = norm(A * x2_original, 'fro');
    % 
    % for i = 1:size(A, 2)
    %     A_removed_i = A;
    %     disp(size(A));
    %     A_removed_i(:, i) = [];  % Remove the i-th column
    %     disp(size(A_removed_i));
    %     [~, ~, V2_i] = svd(A_removed_i);
    %     x2_i = V2_i(:, end);
    %     y2_i = A_removed_i * x2_i;
    %     n_best_i = norm(y2_i, 'fro');
    %     disp(string(i) + ": " + string(n_best_i/n_best_original));
    % end
    % 
    % solution = 0;
    % disp(n_best)
    % for i=1:length(function_name)
    %     if abs(norm(x2(i)*sequences(i), 'fro')) > 0.1
    %         formatted_x = sprintf('%.2f', x2(i));
    %         solution = solution + function_name{i} * str2double(formatted_x);
    %     end
    % 
    % end
    % 
    % disp(solution);
    % A(:,abs(x2) == min(abs(x2))) = [];
    % while ~isempty(A)
    %     [~, ~, V2] = svd(A);
    %     x2 = V2(:, end);
    %     y2 = A * x2;
    %     n = norm(y2, 'fro');
    %     if n < n_best
    %         n_best = n;
    %     end
    %     A(:,abs(x2) == min(abs(x2))) = [];
    % end
    fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b")
    fprintf("Finished calculating best null in input\n")
    %-------Substituting function values into the recursive equation-------
    [N, M] = size(function_values);
    M = M-rec_degree;
    substituted = sym(zeros(N, M));
    fprintf("Calculating substitution(")
    idx = 1:rec_degree;
    for i = 1:M
        fprintf('%5d of %5d)\n', i, M);
        n = (sym(i)-1)+rec_degree;
        s = recursive_relation(n, function_values(:,i+(rec_degree-idx)));
        substituted(:, i) = s;
        fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    end
    fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b")
    disp("Finished calulating substitution in recursive equation");
    fprintf("Calculating results");
    %--------------------------------Result--------------------------------
    % Finding a vector that minimizes the matrix of the substituted values
    matrix = substituted';
    % Adding a column to the matrix in case of a nonhomogeneous equation
    if nonHomogeneous ~= 0
        nonHColumn = zeros(M, 1);
        for i = 1:M
            n = sym(i)-1+rec_degree;
            nonHColumn(i) = subs(nonHomogeneous, x, n);
        end
        matrix = [matrix, nonHColumn];
    end
    [~, ~, V] = svd(matrix);
    x = V(:, end);
    y = matrix * x;
    if nonHomogeneous ~= 0
        x = x(1:end-1);
    end
    y3 = function_values' * x;
    if verbose
        fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        disp("---------------------Results---------------------")
        disp("How good is the best null? (larger values are better) " + string(n_best));
        disp("How good is the solution? (smaller values are better) " + string(norm(y, 'fro')));
        disp("How close is the solution to null? (larger values are better) " + string(norm(y3, 'fro')));
    end
    solution = x;
    example = 0;
    first = true;
    lead = 0;
    for i=1:length(function_name)
        if abs(norm(x(i)*function_values(i, :)', 'fro')) > 0.1
            if first
                lead = x(i);
                first = false;
            end
            x(i)=x(i)/lead;
            formatted_x = sprintf('%.2e', x(i));
            example = example + function_name{i} * str2double(formatted_x);
        end
    end
    disp(example);
end