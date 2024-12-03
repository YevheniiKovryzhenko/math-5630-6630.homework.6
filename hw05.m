% Author: Yevhenii Kovryzhenko / yzk0058@auburn.edu
% Date: 2024-09-01
% Assignment Name: hw05

classdef hw05
    methods (Static)

        function ret = p1(data, powers)
            m = length(data);     % Number of data points
            num_powers = length(powers); % Number of powers given (should be m-1)
        
            % Check that the number of powers is correct
            if num_powers ~= m - 1
                error('The number of powers should be one less than the number of data points.');
            end
        
            % Set up the system of equations
            A = zeros(m, m);
            b = zeros(m, 1);
            b(1) = 1; % First element of b is 1 for the extrapolation condition
            
            % Construct the matrix A based on the powers
            for i = 1:m
                A(i, 1) = 1; % First column corresponds to the constant term f(0)                
                for j = 1:num_powers
                    if (i == 1)
                        A(i, j+1) = 1;
                    else
                        A(i, j+1) = 2^(-(j) * powers(i-1)); % Powers of 2 according to the error terms
                    end
                end
            end
        
            % Display matrix A and vector b for debugging
            % disp('Matrix A:');
            % disp(A);
            % disp('Vector b:');
            % disp(b);
            
            % Solve for the coefficients
            coeffs = A \ b;
        
            % Display coefficients for debugging
            % disp('Coefficients:');
            % disp(coeffs);
        
            % Calculate the extrapolated value for f(0) using the weighted sum
            ret = sum(coeffs .* data(:));
        end

        function R = p2(beta)
            % Compute the value of the series 
            %   sum_{k=0}^(\infty) ((-1)^k /(2k + 1)^{beta})
            % using Richardson extrapolation for faster convergence.
            %
            %:param: beta: a real number on (0, 1].
            %:return: R: the value of the series
        
            % Set parameters for extrapolation
            maxN = 3000;         % Base number of terms to sum initially
            numExtrap = 13;       % Number of levels of extrapolation (useful levels)
            alpha = 0.2;           % Approximate error rate (assume 1 for the series)
        
            % Compute the partial sums S_N, S_{2N}, S_{4N}, etc.
            partial_sums = zeros(1, numExtrap);
            for i = 1:numExtrap
                N = maxN * 2^(i-1);
                partial_sums(i) = compute_partial_sum(N, beta);
            end
        
            % Apply Richardson extrapolation
            for k = 1:(numExtrap - 1)
                for j = 1:(numExtrap - k)
                    partial_sums(j) = (2^(alpha * k) * partial_sums(j+1) - partial_sums(j)) / (2^(alpha * k) - 1);
                end
            end
        
            % The extrapolated value is the first element
            R = partial_sums(1);


            function sum_val = compute_partial_sum(N, beta)
                % Compute the partial sum of the series for a given N and beta
                sum_val = 0;
                for k_ = 0:(N-1)
                    sum_val = sum_val + ((-1)^k_ / (2 * k_ + 1)^beta);
                end
            end
        end




        function coefs = p3(shifts)
            % Compute the coefficients of the finite difference scheme for f'(x)
            % using Lagrange interpolation.
            %
            % f'(x) \approx \frac{1}{h} (c_0 f(x_0) + c_1 f(x_1) + ... + c_n f(x_n)) + O(h^n)
            %
            %:param shifts: a vector of shifts (a_0, a_1, ..., a_n), 
            %               the nodes are x_i = x + a_i h
            %:return coefs: a vector of coefficients (c_0, c_1, ..., c_n)
        
            m = length(shifts);
            coefs = zeros(m, 1);  % Initialize coefficient vector
            
            % Loop over each shift a_i to determine its coefficient
            for i = 1:m
                % Initialize Lagrange basis polynomial derivative at x
                L_prime_at_x = 0;
                
                % Compute the Lagrange basis polynomial L_i(x)
                for j = 1:m
                    if j ~= i
                        % Compute the product terms for the Lagrange polynomial derivative
                        prod = 1 / (shifts(i) - shifts(j));  % Start of the product for d/dx L_i(x)
                        
                        for k = 1:m
                            if k ~= i && k ~= j
                                prod = prod * (0 - shifts(k)) / (shifts(i) - shifts(k));
                            end
                        end
                        
                        % Accumulate the result for the derivative of L_i(x) at x
                        L_prime_at_x = L_prime_at_x + prod;
                    end
                end
                
                % Store the coefficient for f(x_i)
                coefs(i) = L_prime_at_x;
            end
        end

        function coefs = p4(shifts, l)
            % Compute the coefficients of the finite difference scheme for f^{(l)}(x)
            % using the formula
            % f^{(l)}(x) \approx \frac{1}{h^l} (c_0 f(x_0) + c_1 f(x_1) + ... + c_n f(x_n)) + O(h^{n + 1 - l})
        
            % :param: shifts: a vector of shifts (a_0, a_1, ..., a_n), the nodes are x_i = x + a_i * h
            % :param: l: the order of the derivative
            % :return: coefs: a vector of coefficients (c_0, c_1, ..., c_n)
        
            m = length(shifts);
            coefs = zeros(m, 1);
        
            % Set up the system of equations
            % A matrix of powers of shifts
            A = zeros(m, m);
            for i = 1:m
                for j = 1:m
                    A(i, j) = shifts(j)^(i - 1);  % Powers of shifts
                end
            end
        
            % Right-hand side vector, representing the factorials for l-th derivative
            b = zeros(m, 1);
            b(l + 1) = factorial(l);
        
            % Solve the linear system to get the coefficients
            coefs = A \ b;
        end

    end
end