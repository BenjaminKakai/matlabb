function BlasiusSolver()
    % Define the domain
    eta_max = 10;
    num_points = 1000;
    eta = linspace(0, eta_max, num_points);
    
    % Solve using shooting method
    [f_4, f_prime_4, f_double_prime_4] = blasius_shooting(eta, 'pade4');
    [f_6, f_prime_6, f_double_prime_6] = blasius_shooting(eta, 'pade6');
    [f_44, f_prime_44, f_double_prime_44] = blasius_shooting(eta, 'pade44');
    
    % Debug: Print values to check for non-positive values
    disp('f_prime_4 values:');
    disp(f_prime_4);
    disp('f_prime_6 values:');
    disp(f_prime_6);
    disp('f_prime_44 values:');
    disp(f_prime_44);
    
    % Sanity checks
    if any(f_prime_4 <= 0)
        warning('f_prime_4 contains non-positive values.');
    end
    if any(f_prime_6 <= 0)
        warning('f_prime_6 contains non-positive values.');
    end
    if any(f_prime_44 <= 0)
        warning('f_prime_44 contains non-positive values.');
    end
    
    % Plot the Blasius boundary layer profiles
    figure;
    
    % Plot f'(η)
    subplot(3, 1, 1);
    plot(eta, f_prime_4, 'b', 'DisplayName', 'Optimized Padé 4');
    hold on;
    plot(eta, f_prime_6, 'r', 'DisplayName', 'Optimized Padé 6');
    plot(eta, f_prime_44, 'g', 'DisplayName', 'Optimized Padé [4/4]');
    title('Velocity Profile: f''(\eta)');
    xlabel('\eta');
    ylabel('f''(\eta)');
    legend show;
    grid on;
    
    % Plot f(η)
    subplot(3, 1, 2);
    plot(eta, f_4, 'b', 'DisplayName', 'Optimized Padé 4');
    hold on;
    plot(eta, f_6, 'r', 'DisplayName', 'Optimized Padé 6');
    plot(eta, f_44, 'g', 'DisplayName', 'Optimized Padé [4/4]');
    title('Stream Function: f(\eta)');
    xlabel('\eta');
    ylabel('f(\eta)');
    legend show;
    grid on;
    
    % Plot f''(η)
    subplot(3, 1, 3);
    plot(eta, f_double_prime_4, 'b', 'DisplayName', 'Optimized Padé 4');
    hold on;
    plot(eta, f_double_prime_6, 'r', 'DisplayName', 'Optimized Padé 6');
    plot(eta, f_double_prime_44, 'g', 'DisplayName', 'Optimized Padé [4/4]');
    title('Shear Stress: f''''(\eta)');
    xlabel('\eta');
    ylabel('f''''(\eta)');
    legend show;
    grid on;
end
