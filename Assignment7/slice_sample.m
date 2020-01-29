function samples = slice_sample(N_samples,burnIn,log_posterior,bita_ini,widths,step_out,varargin)

D = numel(bita_ini);  % the dimension of the distribution 
samples = zeros(D, N_samples); % samples (NxD)

if numel(widths) == 1  % width (1xD)
    widths = repmat(widths, D, 1);
end
log_Px = feval(log_posterior, bita_ini, varargin{:}); %evaluates a function using its name or its handle, and using the input arguments 

% Main loop
for ii = 1:(N_samples+burnIn)
   % fprintf('Iteration %d                 \r', ii - burnIn);
    log_uprime = log(rand) + log_Px;

    % Sweep through axes (simplest thing)
    for dd = 1:D
        x_l = bita_ini;
        x_r = bita_ini;
        xprime = bita_ini;

        % Create a horizontal interval (x_l, x_r) enclosing xx
        rr = rand;
        x_l(dd) = bita_ini(dd) - rr*widths(dd);
        x_r(dd) = bita_ini(dd) + (1-rr)*widths(dd);
        if step_out
            % Typo in early editions of book. Book said compare to u, but it should say u'
            while (feval(log_posterior, x_l, varargin{:}) > log_uprime)
                x_l(dd) = x_l(dd) - widths(dd);
            end
            while (feval(log_posterior, x_r, varargin{:}) > log_uprime)
                x_r(dd) = x_r(dd) + widths(dd);
            end
        end

        % Inner loop:
        % Propose xprimes and shrink interval until good one found
        zz = 0;
        while 1
            zz = zz + 1;
            %fprintf('Iteration %d   Step %d       \r', ii - burnIn, zz);
            xprime(dd) = rand()*(x_r(dd) - x_l(dd)) + x_l(dd);
            log_Px = feval(log_posterior, xprime, varargin{:});
            if log_Px > log_uprime
                break % this is the only way to leave the while loop
            else
                % Shrink in
                if xprime(dd) > bita_ini(dd)
                    x_r(dd) = xprime(dd);
                elseif xprime(dd) < bita_ini(dd)
                    x_l(dd) = xprime(dd);
                else
                    error('BUG DETECTED: Shrunk to current position and still not acceptable.');
                end
            end
        end
        bita_ini(dd) = xprime(dd);
    end

    % Record samples
    if ii > burnIn
        samples(:, ii - burnIn) = bita_ini(:);
    end
end
%fprintf('\n');
