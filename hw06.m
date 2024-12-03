% Author: Yevhenii Kovryzhenko / yzk0058@auburn.edu
% Date: 2024-09-01
% Assignment Name: hw06


classdef hw06
    methods (Static)
        
        % Problem 1

        function ret = p1(func, a, b, n, option)
            % Implement composite quadrature rules for numerical integration
            % of a function over the interval [a, b] with n subintervals.
            % The option parameter determines the type of quadrature rule to use: 
            % 1 for the midpoint rule, 2 for the trapezoidal rule, and 3 for Simpson's rule.

            %:param func: The function to integrate, provided as a function handle.
            %:param a: The lower bound of the integration interval.
            %:param b: The upper bound of the integration interval.
            %:param n: The number of subintervals to use for the integration.
            %:param option: The type of quadrature rule to use (1, 2, or 3).
            %:return: The approximate integral of the function over the interval [a, b].

            % Calculate step size
            h = (b - a) / n;
            x = linspace(a, b, n + 1); % Partition the interval [a, b]
        
            % Composite Midpoint Rule
            if option == 1
                midpoints = (x(1:end-1) + x(2:end)) / 2; % Midpoints of each subinterval
                ret = h * sum(arrayfun(func,midpoints));
        
            % Composite Trapezoidal Rule
            elseif option == 2
                ret = (h / 2) * (func(a) + 2 * sum(arrayfun(func,x(2:end-1))) + func(b));
        
            % Composite Simpson's Rule
            elseif option == 3
                if mod(n, 2) ~= 0
                    error('n must be even for Simpson''s rule');
                end
                odd_sum = sum(arrayfun(func,x(2:2:end-1)));  % Odd indices
                even_sum = sum(arrayfun(func,x(3:2:end-2))); % Even indices
                ret = (h / 3) * (func(a) + 4 * odd_sum + 2 * even_sum + func(b));
        
            % Invalid Option
            else
                error('Invalid option: %d', option);
            end
        end


        % Problem 2

        function p2()
            % run with the following command: hw06.p2(). Do not edit this function.
            %
            % It checks the convergence of the composite quadrature rules implemented in p1.
            %
            % Here we use some examples, 
            % f_1(x) = exp(x) 
            % f_2(x) = (1 - x^2)^3, this function's first 2 derivatives at the endpoints are zero.
            % f_3(x) = (1 - x^2)^5, this function's first 4 derivatives at the endpoints are zero.
            % f_4(x) = (1 - x^2)^7, this function's first 6 derivatives at the endpoints are zero.

            % Run this function will plot the figures for the convergence of the composite quadrature rules.
            % Make comments about the results obtained from the plots. 
            %
            % > For instance, you can comment on the convergence rates of the quadrature rules, and how they compare to the theoretical rates.
            % > Here are a few example questions you can answer in your comments:
            % > Does the Runge phenomenon of f1 (Runge's function) lower the convergence rate?
            % > Does Simpson's rule have a faster convergence rate than the other two for these examples?
            % > Based on your observations, what kind of functions can have a faster convergence rate with the given composite quadrature rules?

            % Write your comments here.
            %
            % As the order in increased, the error is decreased and fewer
            % subintervals are required to mathc the same error.
            % Overall, second order convergence is the worst, however,
            % Simpson and trapesoidal rules may be performing even worse at low number of subintervals
            % All methods have linear trends for exponental function
            % (f1),with Simpson's rule being better than midpoint and
            % trapezoidal rules.
            %

            f = {  @(x)exp(x),  @(x) (1 - x.^2 ).^3, @(x)(1 - x.^2).^5,  @(x) (1 - x.^2).^7} ;  % Define the integrand
            exact = [exp(1) - exp(-1), 32/35, 512/693 , 4096/6435];  % Define the exact integral
            n = 2.^(1:8);  % Define the number of subintervals
            for k = 1 : length(f)

                error = zeros(3, length(n));  % Initialize the error matrix with zeros

                % Calculate the approximate integral and the error for each quadrature rule and number of subintervals
                for i = 1 : length(n)
                    error(1, i) = abs(hw06.p1(f{k},-1, 1, n(i), 1) - exact(k));
                    error(2, i) = abs(hw06.p1(f{k},-1, 1, n(i), 2) - exact(k));
                    error(3, i) = abs(hw06.p1(f{k},-1, 1, n(i), 3) - exact(k));
                end

                % Plot the error against the number of subintervals using a log-log scale
                figure(k);
    
                loglog(n, error(1, :), 'r-+', 'LineWidth', 2);
                hold on;
                loglog(n, error(2, :), 'g-d', 'LineWidth', 2);
                loglog(n, error(3, :), 'b-x', 'LineWidth', 2);

                loglog(n, 1./ n.^2, 'm--', 'LineWidth', 1);
                loglog(n, 1./ n.^4, 'k-.', 'LineWidth', 1);
                loglog(n, 1./ n.^6, 'm--d', 'LineWidth', 1);
                loglog(n, 1./ n.^8, 'k--o', 'LineWidth', 1);

                xlabel('Number of subintervals');
                ylabel('Absolute error');
                title(sprintf('Convergence of composite quadrature rules for %s', functions(f{k}).function));
                legend('Midpoint rule', 'Trapezoidal rule', 'Simpson''s rule', '2nd order convergence', '4th order convergence', '6th order convergence', '8th order convergence', 'Location', 'best');
                grid on;
                hold off;
            end

        end

        
        % Problem 3

        function ret = p3(func, a, b, N, option)
            % Use Richardson extrapolation to implement the Romberg integration method.
            %
            % :param func: The function to integrate, provided as a function handle.
            % :param a: The lower bound of the integration interval.
            % :param b: The upper bound of the integration interval.
            % :param N: Maximum number of iterations for Romberg integration (2^N subintervals).
            % :param option: The type of quadrature rule to use (1, 2, or 3). See p1.
            % :return: The approximate integral of the function over the interval [a, b].
        
            % Initialize a list to hold integral approximations
            approx_values = zeros(N, 1);
            
            % Compute integral approximations for 2^k subintervals (k = 1, ..., N)
            for k = 1:N
                num_subintervals = 2^k;
                approx_values(k) = hw06.compute_integral(func, a, b, num_subintervals, option);
            end
        
            % Determine the powers for Richardson extrapolation
            if option == 1 || option == 2
                powers = 2:2:(2 * (N - 1));
            elseif option == 3
                powers = 4:2:(2 * (N - 1) + 2);
            else
                error('Invalid option: %d', option);
            end
        
            % Apply Richardson extrapolation using the provided p1 function
            ret = hw05.p1(approx_values, powers);
        end

        function integral = compute_integral(func, a, b, n, option)
            % Helper function to compute the integral approximation
            % using the specified quadrature rule and number of subintervals.
            h = (b - a) / n; % Width of each subinterval
            x = linspace(a, b, n + 1); % Subinterval boundaries
        
            switch option
                case 1 % Midpoint rule
                    midpoints = (x(1:end-1) + x(2:end)) / 2;
                    integral = sum(func(midpoints)) * h;
                case 2 % Trapezoidal rule
                    integral = (h / 2) * (func(a) + 2 * sum(func(x(2:end-1))) + func(b));
                case 3 % Simpson's rule
                    if mod(n, 2) ~= 0
                        error('n must be even for Simpson''s rule');
                    end
                    integral = (h / 3) * (func(a) + 4 * sum(func(x(2:2:end-1))) + 2 * sum(func(x(3:2:end-2))) + func(b));
                otherwise
                    error('Invalid quadrature rule option: %d', option);
            end
        end

        % % Problem 4
        function ret = p4()
            % Construct the Gauss quadrature rule using the roots of the Legendre polynomial of degree 6.
            %
            % :return: A 6x2 matrix containing the roots and weights of the Gauss quadrature rule. 
            %          The first column contains the roots and the second column contains the corresponding weights.
        
            % Define the Legendre polynomial of degree 6
            legendre_6 = @(x) hw06.legendre_poly_6(x);
        
            % Define the derivative of the Legendre polynomial (for Newton's method)
            legendre_6_derivative = @(x) (1386 * x^5 - 1260 * x^3 + 210 * x) / 16;
        
            % Use symmetry of roots (roots are symmetric about 0)
            roots = zeros(6, 1);
            roots_guess = [-0.9, -0.5, 0]; % Initial guesses for the roots in the interval [-1, 0]
        
            % Find the roots in [0, 1] using fzero
            for i = 1:3
                roots(i) = fzero(legendre_6, roots_guess(i)); % Find roots in [-1, 0]
            end
        
            % Mirror the roots for [0, 1] to [0, -1]
            roots(4:6) = -roots(1:3);
        
            % Sort roots for convenience
            roots = sort(roots);
        
            % Compute the weights for each root
            weights = zeros(6, 1);
            for i = 1:6
                weights(i) = 2 / ((1 - roots(i)^2) * (legendre_6_derivative(roots(i)))^2);
            end
        
            % Combine roots and weights into a matrix
            ret = [roots, weights]; % Return the roots and weights of the Gauss quadrature rule
        end

        % Helper function for Legendre polynomial of degree 6
        function val = legendre_poly_6(x)
            % Compute the Legendre polynomial of degree 6 at the point x.
            %
            % :param x: The point at which to evaluate the Legendre polynomial.
            % :return: The value of the Legendre polynomial of degree 6 at the point x.
        
            val = (231 * x^6 - 315 * x^4 + 105 * x^2 - 5) / 16;
        end

        function ret = p5(n)
            % Construct the Gauss quadrature rule using the roots of the Legendre polynomial of degree n.
            %
            % :param n: The degree of the Legendre polynomial for the nodes of the Gauss quadrature rule.
            % :return: An nx2 matrix containing the roots and weights of the Gauss quadrature rule.
            %
            % To evaluate the Legendre polynomial or its derivative of a specific degree n, the handles are:
            % @(x) hw06.legendre_poly(n, x) and @(x) hw06.deriv_lengendre_poly(n, x).
        
            % Preallocate arrays for roots and weights
            roots = zeros(n, 1);
            weights = zeros(n, 1);
        
            % Define the Legendre polynomial and its derivative
            legendre_n = @(x) hw06.legendre_poly(n, x);
            legendre_n_derivative = @(x) hw06.deriv_lengendre_poly(n, x);
        
            % Initial guesses for roots in the interval [-1, 1] (e.g., Chebyshev nodes)
            guesses = cos(pi * (4 * (1:n) - 1) / (4 * n + 2));
        
            % Use Newton's method to find roots of the Legendre polynomial
            for i = 1:n
                x = guesses(i);
                for j = 1:100 % Limit iterations to ensure convergence
                    dx = -legendre_n(x) / legendre_n_derivative(x);
                    x = x + dx;
                    if abs(dx) < 1e-12 % Convergence tolerance
                        break;
                    end
                end
                roots(i) = x;
            end
        
            % Compute the weights for each root
            for i = 1:n
                weights(i) = 2 / ((1 - roots(i)^2) * (legendre_n_derivative(roots(i)))^2);
            end
        
            % Sort the roots and weights for consistency
            [roots, idx] = sort(roots);
            weights = weights(idx);
        
            % Combine roots and weights into a matrix
            ret = [roots, weights]; % Return the roots and weights of the Gauss quadrature rule
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                                                             %
        % Helper functions below. Do not modify. You can create your own helper functions if needed.                  %
        %                                                                                                             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function test_p5()
            % Unit tests for the p5 function
            fprintf('Running tests for p5...\n');
            
            % Tolerance for floating-point comparisons
            tol = 1e-12;
        
            %% Test 1: Basic functionality for n = 2
            result = hw06.p5(2);
            assert(size(result, 1) == 2 && size(result, 2) == 2, 'Test 1 Failed: Incorrect size for n=2');
            disp('Test 1 Passed');
        
            %% Test 2: Symmetry of roots for n = 4
            result = hw06.p5(4);
            roots = result(:, 1);
            assert(all(abs(sort(roots) + flip(sort(roots))) < tol), 'Test 2 Failed: Roots are not symmetric for n=4');
            disp('Test 2 Passed');
            
            %% Test 3: Sum of weights for n = 3 (Should equal 2)
            result = hw06.p5(3);
            weights = result(:, 2);
            assert(abs(sum(weights) - 2) < tol, 'Test 3 Failed: Sum of weights is not 2 for n=3');
            disp('Test 3 Passed');
        
            %% Test 4: Integration of f(x) = 1 over [-1, 1] for n = 5
            result = hw06.p5(5);
            roots = result(:, 1);
            weights = result(:, 2);
            integral = sum(weights .* ones(size(roots)));
            assert(abs(integral - 2) < tol, 'Test 4 Failed: Integration of f(x)=1 is incorrect for n=5');
            disp('Test 4 Passed');
        
            %% Test 5: Integration of f(x) = x^2 over [-1, 1] for n = 6
            result = hw06.p5(6);
            roots = result(:, 1);
            weights = result(:, 2);
            integral = sum(weights .* (roots.^2));
            true_value = 2/3; % Exact integral of x^2 over [-1, 1]
            assert(abs(integral - true_value) < tol, 'Test 5 Failed: Integration of f(x)=x^2 is incorrect for n=6');
            disp('Test 5 Passed');
        
            %% Test 6: Integration of f(x) = x^4 over [-1, 1] for n = 6
            result = hw06.p5(6);
            roots = result(:, 1);
            weights = result(:, 2);
            integral = sum(weights .* (roots.^4));
            true_value = 2/5; % Exact integral of x^4 over [-1, 1]
            assert(abs(integral - true_value) < tol, 'Test 6 Failed: Integration of f(x)=x^4 is incorrect for n=6');
            disp('Test 6 Passed');
        
            %% Test 7: Symmetry of weights for n = 8
            result = hw06.p5(8);
            weights = result(:, 2);
            assert(all(abs(weights - flip(weights)) < tol), 'Test 7 Failed: Weights are not symmetric for n=8');
            disp('Test 7 Passed');
        
            %% Test 8: Roots within the interval [-1, 1] for n = 10
            result = hw06.p5(10);
            roots = result(:, 1);
            assert(all(roots >= -1 & roots <= 1), 'Test 8 Failed: Roots are not within [-1, 1] for n=10');
            disp('Test 8 Passed');
        
            %% Test 9: Integration of f(x) = x^3 (odd function, should be 0) for n = 7
            result = hw06.p5(7);
            roots = result(:, 1);
            weights = result(:, 2);
            integral = sum(weights .* (roots.^3));
            assert(abs(integral) < tol, 'Test 9 Failed: Integration of f(x)=x^3 is incorrect for n=7');
            disp('Test 9 Passed');
        
            %% Test 10: Larger n = 20 and sum of weights
            result = hw06.p5(20);
            weights = result(:, 2);
            assert(abs(sum(weights) - 2) < tol, 'Test 10 Failed: Sum of weights is not 2 for n=20');
            disp('Test 10 Passed');
        
            fprintf('All tests passed successfully!\n');
        end


        % Helper functions for p5. The following function is used to evaluate the Legendre polynomial of degree n.
        function val = legendre_poly(n, x)
            % Compute the nth Legendre polynomial P_n at the point x.
            %
            % :param n: The degree of the Legendre polynomial.
            % :param x: The point at which to evaluate the Legendre polynomial.
            % :return: The value of the nth Legendre polynomial at the point x.

            if (n == 0)
                val = 1;
            elseif (n == 1)
                val = x;
            else
                val = hw06.legendre_poly(n-1, x) * x * (2 * n - 1)/n - (n - 1) * hw06.legendre_poly(n - 2, x) / n;
            end
        end

        function val = deriv_lengendre_poly(n, x)
            % Compute the derivative of the nth Legendre polynomial P_n at the point x.
            %   
            % :param n: The degree of the Legendre polynomial.
            % :param x: The point at which to evaluate the derivative of the Legendre polynomial.
            % :return: The value of the derivative of the nth Legendre polynomial at the point x.
            val = n / (x^2 - 1) * (x * hw06.legendre_poly(n, x) - hw06.legendre_poly(n - 1, x));
        end
    end
end

