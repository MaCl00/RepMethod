function [function_name, solution] = repertoire(func_guess, recursive_relation, nonHomogeneous, precision, start_point, rec_degree, num_poly, varargin)
    digits(precision);
    % Flags
    verbose = any(cellfun( @(x) isequal( x, 'verbose' ), varargin ));
    continue_search = ~any(cellfun( @(x) isequal( x, 'quitSearch' ), varargin ));
    multiply = any(cellfun( @(x) isequal( x, 'addMul' ), varargin ));
    mono = any(cellfun( @(x) isequal( x, 'mono' ), varargin ));
    num_func = length(func_guess);
    syms x;
    % The number of required data points
    polyMul = num_poly;
    if num_poly == 0
        polyMul = 1;
    end
    num_substitution = 1 + polyMul * ((num_func + num_func^2) / 2 + 1);
    if ~multiply
        num_substitution = 1 + polyMul * (num_func + 1);
    end
    if num_poly == 0
        num_substitution = num_substitution - 1;
    end
    if nonHomogeneous ~= 0
        num_substitution = num_substitution +1;
    end
    summing = false;
    if rec_degree == -1
        rec_degree = 2;
        summing = true;
    end
    num_function_values = num_substitution + rec_degree;
    % ---------------------Calculating function-values---------------------
    range = num2cell(start_point:start_point+num_function_values-1);
    % A matrix where enties are the input functions evaluated at the points
    if mono
        [function_values, function_name] = generateGuessesMono(func_guess, range, num_poly, verbose, multiply);
    else
        [function_values, function_name] = generateGuesses(func_guess, range, num_poly, verbose, multiply);
    end
    %-------Substituting function values into the recursive equation-------
    [N, M] = size(function_values);
    M = M-rec_degree;
    substituted = sym(zeros(N, M));
    if verbose
        fprintf("Calculating substitution")
    end
    idx = 1:rec_degree;
    for i = 1:M
        n = (sym(i)-2)+rec_degree + start_point;
        if summing
            s = recursive_relation(n, function_values(:,(i+2)-(1:(i+1))));
        else
            s = recursive_relation(n, function_values(:,i+(rec_degree-idx)));
        end
        substituted(:, i) = s;
    end
    if verbose
        fprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b")
        disp("Finished calulating substitution in recursive equation");
        disp("Extracting Solutions");
    end
    %--------------------------------Result--------------------------------
    % Finding a vector that minimizes the matrix of the substituted values
    matrix = substituted';
    % Adding a column to the matrix in case of a nonhomogeneous equation
    hasNonH = false;
    solution = [];
    function_values = function_values';
    
    if nonHomogeneous ~= 0
        hasNonH = true;
        nonHColumn = zeros(M, 1);
        for i = 1:M
            n = sym(i)-2+rec_degree+start_point;
            nonHColumn(i) = subs(nonHomogeneous, x, n);
        end
        solution = [solution, zeros(size(matrix, 2), 1)];
        matrix = [matrix, nonHColumn];
    end
    is_valid_solution = true;
    selected_columns = 1:size(matrix, 2);
    while is_valid_solution
        reduced_matrix = matrix(:,selected_columns);
        [~, ~, V] = svd(reduced_matrix);
        x = V(:, end);
        y = eval(log10(norm(reduced_matrix * x, "fro")));
        if verbose
            if y < -precision/2
                disp("Confirmed a solution, search starts.");
            elseif continue_search
                if y < 0
                    disp("Low likelihood of further solutions within the matrix. Search will continue. Expect longer run-time.");
                else
                    disp("No further solution expected within the matrix. Search will continue. Expect long run-time.");
                end
            else
                disp("Exiting search. No further solutions expected within the matrix. If you want to continue the search anyway, remove the flag 'quitSearch'.")
            end
            
        end
        if ~continue_search && y >= -precision/2
            break;
        end
        [result, is_valid_solution] = SolutionExtraction(reduced_matrix, x, precision, 0, verbose);
        solution_vec = zeros(size(matrix, 2), 1);
        solution_vec(selected_columns) = result;
        if is_valid_solution
            [max_value, index] = max(abs(result));
            if hasNonH && result(end) ~= 0
                hasNonH = false;
                index = length(result);
                solution_vec = solution_vec./-solution_vec(end);
                solution_vec(end) = [];
                result(end) = [];
                solution(:,1) = solution_vec;
            else
                if nonHomogeneous ~= 0
                    solution_vec(end) = [];
                end
                % null check
                close_to_null = log10(norm(function_values*solution_vec, 'fro'));
                if close_to_null < -precision/2
                    if verbose
                        disp("Found a linear combination that is zero in the guesses:")
                        example = 0;
                        solution_vec = solution_vec ./ max_value;
                        for i=1:length(solution_vec)
                            if solution_vec(i) ~= 0
                                formatted_x = sprintf('%.2e', solution_vec(i));
                                example = example + function_name{i} * str2double(formatted_x);
                            end
                        end
                        disp(example);
                        disp("Try to avoid this, run time can increase drastically.")
                    end
                    for i=1:size(solution, 2)
                        coef = solution(index,i);
                        if coef ~= 0
                            solution(:,i) = solution(:,i) - coef/max_value * solution_vec;
                        end
                    end
                    selected_columns(index) = [];
                    continue;
                end
                solution = [solution, solution_vec];
            end
            selected_columns(index) = [];
            if verbose
                disp("Calculated a solution: ");
                example = 0;
                solution_vec = solution_vec ./ max_value;
                for i=1:length(solution_vec)
                    if solution_vec(i) ~= 0
                        formatted_x = sprintf('%.2e', solution_vec(i));
                        example = example + function_name{i} * str2double(formatted_x);
                    end
                end
                disp(example);
            end
        elseif verbose
            disp("Found no more solution, end of search.");
        end
    end
    rows_to_keep = true(1, size(solution, 1));
    for i=1:length(function_name)
        if size(solution, 2)>0 && all(solution(i,:) == 0)
            rows_to_keep(i) = false;
        end
    end
    solution = solution(rows_to_keep, :);
    function_name = function_name(rows_to_keep);
    % if 1 == 2
    %     y = matrix * solution;
    %     if nonHomogeneous ~= 0
    %         solution = solution(1:end-1);
    %     end
    %     y3 = function_values' * solution;
    %     example = 0;
    %     solution = solution ./ max(abs(solution));
    %     for i=1:length(solution)
    %         if solution(i) ~= 0
    %             formatted_x = sprintf('%.2e', solution(i));
    %             example = example + function_name{i} * str2double(formatted_x);
    %         end
    %     end
    %     if verbose
    %         disp("---------------------Results---------------------")
    %         if checkNull
    %             disp("How good is the best null? (larger values are better) " + sprintf('%.2e', n_best));
    %         end
    %         disp("How good is the solution? (smaller values are better) " + sprintf('%.2e', norm(y, 'fro')));
    %         disp("How close is the solution to null? (larger values are better) " +  sprintf('%.2e', norm(y3, 'fro')));
    %         disp(example);
    %     end
    % end
end