function BlasiusSolver()
    % Define the domain
    eta_max = 10;
    num_points = 1000;
    eta = linspace(0, eta_max, num_points);
    
    % Solve using shooting method for different schemes
    [f_4, f_prime_4, f_double_prime_4] = blasius_shooting(eta, 'pade4');
    [f_6, f_prime_6, f_double_prime_6] = blasius_shooting(eta, 'pade6');
    [f_44, f_prime_44, f_double_prime_44] = blasius_shooting(eta, 'pade44');
    
    % Debug: Check values at specific points
    disp('Sample values for f_prime_4:');
    disp(f_prime_4(1:10)); % Display the first 10 values
    disp('Sample values for f_prime_6:');
    disp(f_prime_6(1:10)); % Display the first 10 values
    disp('Sample values for f_prime_44:');
    disp(f_prime_44(1:10)); % Display the first 10 values
    
    % Plot the Blasius boundary layer profiles
    figure;
    
    % Plot f'(η)
    subplot(3, 1, 1);
    plot(eta, f_prime_4, 'b', 'DisplayName', 'Padé 4');
    hold on;
    plot(eta, f_prime_6, 'r', 'DisplayName', 'Padé 6');
    plot(eta, f_prime_44, 'g', 'DisplayName', 'Padé [4/4]');
    title('Velocity Profile: f''(\eta)');
    xlabel('\eta');
    ylabel('f''(\eta)');
    legend show;
    grid on;
    
    % Plot f(η)
    subplot(3, 1, 2);
    plot(eta, f_4, 'b', 'DisplayName', 'Padé 4');
    hold on;
    plot(eta, f_6, 'r', 'DisplayName', 'Padé 6');
    plot(eta, f_44, 'g', 'DisplayName', 'Padé [4/4]');
    title('Stream Function: f(\eta)');
    xlabel('\eta');
    ylabel('f(\eta)');
    legend show;
    grid on;
    
    % Plot f''(η)
    subplot(3, 1, 3);
    plot(eta, f_double_prime_4, 'b', 'DisplayName', 'Padé 4');
    hold on;
    plot(eta, f_double_prime_6, 'r', 'DisplayName', 'Padé 6');
    plot(eta, f_double_prime_44, 'g', 'DisplayName', 'Padé [4/4]');
    title('Shear Stress: f''''(\eta)');
    xlabel('\eta');
    ylabel('f''''(\eta)');
    legend show;
    grid on;
end

function [f, f_prime, f_double_prime] = blasius_shooting(eta, method)
    % Initial conditions
    f0 = 0;
    f_prime0 = 0;
    f_double_prime0 = 0.3320573362151963; % The Blasius value at eta = 0

    % Set the differential equation
    ode = @(eta, F) [F(2); F(3); -0.5 * F(1) * F(3)];

    % Solve using ode45
    [eta, F] = ode45(ode, eta, [f0, f_prime0, f_double_prime0]);

    % Extract results
    f = F(:, 1);
    f_prime = F(:, 2);
    f_double_prime = F(:, 3);

    % Apply the specific Padé approximation method
    switch method
        case 'pade4'
            f = pade_approx(f, 4, 4);
            f_prime = pade_approx(f_prime, 4, 4);
            f_double_prime = pade_approx(f_double_prime, 4, 4);
        case 'pade6'
            f = pade_approx(f, 6, 6);
            f_prime = pade_approx(f_prime, 6, 6);
            f_double_prime = pade_approx(f_double_prime, 6, 6);
        case 'pade44'
            f = pade_approx(f, 4, 4);
            f_prime = pade_approx(f_prime, 4, 4);
            f_double_prime = pade_approx(f_double_prime, 4, 4);
        otherwise
            error('Unknown method: %s', method);
    end
end

function approx = pade_approx(data, p, q)
    % Padé approximation function
    % Calculate Padé approximation coefficients
    c = zeros(p + q + 1, 1);
    c(1) = 1;
    for k = 1:p
        c(k + 1) = c(k) * (p - k + 1) / (k * (q + k));
    end
    
    % Calculate Padé approximation
    n = length(data);
    approx = zeros(n, 1);
    for i = 1:n
        approx(i) = 0;
        for j = 1:min(i, p + 1)
            approx(i) = approx(i) + c(j) * data(i - j + 1);
        end
    end
end
