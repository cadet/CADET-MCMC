
classdef OptAlgo < handle
% =============================================================================
% This is the class of the functions of OPTimization ALGOrithms (OptAlgo).
%
% =============================================================================

    properties (Constant)
        swarm       = 20;              % population size
        sample      = 10000;           % loop count for evolution (at least 100)
        dataPoint   = 10000;           % amount of data observation
        prior       = [];              % prior information (default empty)
        logScale    = true;            % log-scale for parameter domain (default true)
%        prior = load('prior.dat');     % load the prior distribution if possible
    end


    methods (Static = true, Access = 'public')

        function [xValue, yValue] = Markov_Chain_Monte_Carlo(FUNC, opt)
%------------------------------------------------------------------------------
% Markov Chain Monte Carlo (MCMC) simulation
%
% The central problem is that of determining the posterior probability for
% parameters given the data, p(thelta|y). With uninformative prior distribution
% and normalization constant p(y), the task is reduced to that of maximizing the
% likelihood of data, p(y|thelta).
%
% "Delayed rejection" was implemented in order to increase the acceptance ratio.
%
% If the Jacobian matrix can be obtained, it will be used to generate R
%
% Parameter:
%       - FUNC. The objective function from your field
%
% Returns:
%       - xValue. The optimal parameter set
%       - yValue. The value of objective function with optimal parameter set
%------------------------------------------------------------------------------


            startTime = clock;

            if nargin < 2
                error('OptAlgo.MCMC: There are no enough inputs \n');
            end

            % Get the MCMC options
            opt = OptAlgo.getOptionsMCMC(opt);

            % Preallocation
            accepted = 0; n0 = 1;
            lasti = 0; chaincov = []; chainmean = []; wsum = [];
            chain = zeros(opt.nsamples, opt.Nparams+1);

            % Inilization of a starting point
            [oldpar, SS, Jac] = OptAlgo.initPointMCMC(FUNC, opt);

            % Construct the error standard deviation: sigma square
            if SS < 0, sumSquare = exp(SS); else, sumSquare = SS; end
            sigmaSqu = sumSquare / (opt.nDataPoint - opt.Nparams);
            sigmaSqu_0 = sigmaSqu;

            if ~isempty(Jac)

                % Construct the R = chol(covMatrix) for candidate generatio
                [~, S, V] = svd(Jac);

                % Truncated SVD
                if opt.Nparams == 1, S = S(1); else, S = diag(S); end
                S(S < 1e-3) = 1e-3;

                % Fisher information matrix
                % covariance matrix = sigmaSqu * inv(Fisher information matrix) = v'*s^{-2}*v
                covMatrix = V * diag(1./S.^2)* V' * sigmaSqu;

                % Construct the R = chol(covMatrix) for candidate generation
                R = chol(covMatrix .* 2.4^2 ./ opt.Nparams);

            else

                % If the Jacobian matrix cannot be obtained,
                %   a set of samples is used to generate R
                [R, oldpar, SS] = OptAlgo.burnInSamples(FUNC, opt);

            end

%------------------------------------------------------------------------------

% Main loop
%------------------------------------------------------------------------------
            for j = 1:(opt.nsamples+opt.burn_in)

                if j > opt.burn_in
                    fprintf('Iter: %4d -------- Accept_ratio: %3d%% ---------- Minimum: %g ---------- \n',...
                        j, fix(accepted / j * 100), SS);
                    fprintf('%10.3g | ', OptAlgo.pTransfer('exp', oldpar)); fprintf('\n');
                end

                accept = false;

                % Generate the new proposal point with R
                newpar = oldpar + randn(1, opt.Nparams) * R;

                % Check the boundary limiation
                newpar( newpar < opt.bounds(1, :) ) = opt.bounds(1, newpar < opt.bounds(1, :));
                newpar( newpar > opt.bounds(2, :) ) = opt.bounds(2, newpar > opt.bounds(2, :));

                % Calculate the objective value of the new proposal
                newSS = feval( FUNC, OptAlgo.pTransfer('exp', newpar) );

                % The Metropolis probability
                if OptAlgo.logScale
                    rho12 = exp( -0.5 *((newSS - SS) / sigmaSqu) + sum(newpar) - sum(oldpar) ) * ...
                        OptAlgo.priorPDF(newpar) / OptAlgo.priorPDF(oldpar);
                else
                    rho12 = exp( -0.5 * (newSS - SS) / sigmaSqu ) * ...
                        OptAlgo.priorPDF(newpar) / OptAlgo.priorPDF(oldpar);
                end

                % The new proposal is accepted with Metropolis probability
                if rand <= min(1, rho12)
                    accept   = true;
                    oldpar   = newpar;
                    SS       = newSS;
                    accepted = accepted + 1;
                end

                % If the poposal is denied, a Delayed Rejection procedure is adopted
                %   in order to increase the acceptance ratio
                if ~accept && opt.delayReject

                    % Shrink the searching domain by the factor 1/10
                    newpar2 = oldpar + randn(1, opt.Nparams) * (R ./ 10);

                    % Check the boundary limitation of the new generated point
                    newpar2(newpar2 < opt.bounds(1, :)) = opt.bounds(1, newpar2 < opt.bounds(1, :));
                    newpar2(newpar2 > opt.bounds(2, :)) = opt.bounds(2, newpar2 > opt.bounds(2, :));

                    % Calculate the objective value of the new proposal
                    newSS2 = feval( FUNC, OptAlgo.pTransfer('exp', newpar2) );

                    if OptAlgo.logScale

                        rho32 = exp( -0.5 *((newSS - newSS2) / sigmaSqu) + sum(newpar) - sum(newpar2) ) * ...
                            OptAlgo.priorPDF(newpar) / OptAlgo.priorPDF(newpar2);

                        % The conventional version of calculation
%                        q2 = exp( -0.5 *((newSS2 - SS) / sigmaSqu) + sum(newpar2) - sum(oldpar) ) * ...
%                            OptAlgo.priorPDF(newpar2) / OptAlgo.priorPDF(oldpar);
%                        q1 = exp( -0.5 * (norm((newpar2 - newpar) * inv(R))^2 - norm((oldpar - newpar) * inv(R))^2) );

                        % The speed-up version of above calculation
                        q1q2 = exp( -0.5 *( (newSS2 - SS) / sigmaSqu + ...
                            (newpar2 - newpar) * (R \ (R' \ (newpar2' - newpar'))) - ...
                            (oldpar - newpar) * (R \ (R' \ (oldpar' - newpar'))) ) + ...
                            sum(newpar2) - sum(oldpar) ) * ...
                            OptAlgo.priorPDF(newpar2) / OptAlgo.priorPDF(oldpar);

                    else

                        rho32 = exp( -0.5 * (newSS - newSS2) / sigmaSqu ) * ...
                            OptAlgo.priorPDF(newpar) / OptAlgo.priorPDF(newpar2);

                        % The conventional version of calculation
%                        q2 = exp( -0.5 * (newSS2 - SS) / sigmaSqu ) * ...
%                            OptAlgo.priorPDF(newpar2) / OptAlgo.priorPDF(oldpar);
%                        q1 = exp( -0.5 * (norm((newpar2 - newpar) * inv(R))^2 - norm((oldpar - newpar) * inv(R))^2) );

                        % The speed-up version of above calculation
                        q1q2 = exp( -0.5 *( (newSS2 - SS) / sigmaSqu + ...
                            (newpar2 - newpar) * (R \ (R' \ (newpar2' - newpar'))) - ...
                            (oldpar - newpar) * (R \ (R' \ (oldpar' - newpar'))) ) ) * ...
                            OptAlgo.priorPDF(newpar2) / OptAlgo.priorPDF(oldpar);

                    end

                    rho13 = q1q2 * (1 - rho32) / (1 - rho12);

                    if rand <= min(1, rho13)
                        oldpar   = newpar2;
                        SS       = newSS2;
                        accepted = accepted + 1;
                    end

                end % if ~accept && delayReject

                % During the burn-in period, if the acceptance rate is extremly high or low,
                %   the R matrix is manually adjusted
                if j <= opt.burn_in
                    if mod(j, 50) == 0
                        if accepted/j < 0.05
                            fprintf('Acceptance ratio %3.2f smaller than 5 %%, scaled \n', accepted/j*100);
                            R = R ./ 5;
                        elseif accepted/j > 0.95
                            fprintf('Acceptance ratio %3.2f largeer than 95 %%, scaled \n', accepted/j*100);
                            R = R .* 5;
                        end
                    end
                end

                % After the burn-in period, the chain is stored
                if j > opt.burn_in
                    chain(j-opt.burn_in, 1:opt.Nparams) = oldpar;
                    chain(j-opt.burn_in, opt.Nparams+1) = SS;

                    temp = chain(j-opt.burn_in, :);
                    save('chainData.dat', 'temp', '-ascii', '-append');
                end

                if mod(j-opt.burn_in, opt.convergInt) == 0 && (j-opt.burn_in) > 0
                    % Updata the R according to previous chain
                    [chaincov, chainmean, wsum] = OptAlgo.covUpdate( chain(lasti+1:j-opt.burn_in, 1:opt.Nparams), ...
                        1, chaincov, chainmean, wsum );

                    lasti = j;
                    R = chol(chaincov + eye(opt.Nparams)*1e-7);

                    % Check the convergence condition
                    criterion = OptAlgo.Geweke( chain(1:j-opt.burn_in, 1:opt.Nparams) );
                    if all( abs(criterion) < opt.criterion ) || j == opt.nsamples+opt.burn_in
                        maxIter = j;
                        break
                    end
                end

                % Updata the sigma^2 according to the current objective value
                if SS < 0, sumSquare = exp(SS); else, sumSquare = SS; end
                sigmaSqu  = 1 / OptAlgo.GammarDistribution( 1, 1, (n0 + opt.nDataPoint)/2,...
                    2 / (n0 * sigmaSqu_0 + sumSquare) );

                save('sigmaSqu.dat', 'sigmaSqu', '-ascii', '-append');

            end % for j = 1:opt.nsamples

%------------------------------------------------------------------------------

% Post-process
%------------------------------------------------------------------------------
            clear chain;

            % Generate the population for figure plot
            Population = OptAlgo.conversionDataMCMC(maxIter, opt);

            OptAlgo.FigurePlot(Population, opt);

            [yValue, row] = min(Population(:, opt.Nparams+1));
            xValue = OptAlgo.pTransfer('exp', Population(row, 1:opt.Nparams));

            % Gather some useful information and store them
            result.optTime         = etime(clock, startTime) / 3600;
            result.covMatrix       = R' * R;
            result.Iteration       = maxIter;
            result.criterion       = criterion;
            result.accepted        = fix(accepted/maxIter * 100);
            result.xValue          = xValue;
            result.yValue          = yValue;
            result.sigma           = sqrt(sigmaSqu);

            fprintf('\n****************************************************** \n');
            save(sprintf('result_%2d.mat', fix(rand*100)), 'result');
            fprintf('Time %5.3g hours elapsed after %5d iterations \n', result.optTime, result.Iteration);
            fprintf('The minimal objective value found during sampling is: %10.3g \n', yValue);
            fprintf('The correspondingly optimal set is: ');
            fprintf(' %g |', xValue);
            fprintf('\nThe statistical information of MCMC is stored as result.mat \n');
            fprintf('The historical chain of MCMC is stored as Population.dat \n');
            fprintf('****************************************************** \n');

        end % Markov_Chain_Monte_Carlo

        function opt = getOptionsMCMC(obj)
%------------------------------------------------------------------------------
% The parameters for the Optimizer
%
% Return:
%       - opt.
%           + opt.Nparams. The number of optimized parameters
%           + opt.bounds. 2*nCol matrix of the parameter limitation
%           + opt.nsamples. The pre-defined maximal iteration
%           + opt.criterion. The error tolerance to stop the algorithm
%           + opt.burn_in. The burn-in period before adaptation begin
%           + opt.convergInt. The integer for checking of the convergence
%           + opt.rejectValue. This is specific used in the generation of R matrix
%               when Jacobian information is not available. In the generated sample,
%               the proposals whose objective value is larger than this will be rejected
%           + opt.nDataPoint. The number of observations in the measured data
%           + opt.deylayReject. By default DR= 1
%           + opt.Jacobian. If Jacobian matrix is available, set it to true
%------------------------------------------------------------------------------


            opt = [];

            opt.Nparams       = length(obj.params);
            opt.bounds        = OptAlgo.pTransfer('log', obj.paramBound)';
            opt.nsamples      = OptAlgo.sample;
            opt.criterion     = 0.0001;
            opt.burn_in       = 0;
            opt.convergInt    = 100;
            opt.rejectValue   = 1e5;
            opt.nDataPoint    = OptAlgo.dataPoint;
            opt.Jacobian      = false; % set it to true only when Jacobian matrix is available
            opt.delayReject   = true;

            [row, col] = size(opt.Nparams);
            if row > 1 || col > 1
                error('OptAlgo.getOptionsMCMC: The initialized dimension of the parameter set might be wrong \n');
            end

            [row, col] = size(opt.bounds);
            if col ~= opt.Nparams || row ~= 2
                error('OptAlgo.getOptionsMCMC: Please check your setup of the range of parameters \n');
            end

        end % getOptionsMCMC

        function [oldpar, SS, Jac]= initPointMCMC(FUNC, opt)
%------------------------------------------------------------------------------
% Generate the initially guessing point for the MCMC algorithm
%------------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgo.initPointMCMC: There are no enough input arguments \n');
            end

            oldpar = rand(1,opt.Nparams);

            oldpar(:, 1:opt.Nparams) = opt.bounds(1,:) + ...
                oldpar(:, 1:opt.Nparams) .* (opt.bounds(2,:) - opt.bounds(1,:));

            % Check the boundary limiation
%            oldpar(oldpar < opt.bounds(1, :)) = opt.bounds(1, oldpar < opt.bounds(1, :));
%            oldpar(oldpar > opt.bounds(2, :)) = opt.bounds(1, oldpar > opt.bounds(2, :));

            % Get the resudual value and Jacobian matrix of the guessing point
            Jac = [];
            if opt.Jacobian
                [SS, ~, Jac] = feval( FUNC, OptAlgo.pTransfer('exp', oldpar) );
            else
                SS = feval( FUNC, OptAlgo.pTransfer('exp', oldpar) );
            end

        end % initPointMCMC

        function [R, oldpar, SS] = burnInSamples(FUNC, opt)
%------------------------------------------------------------------------------
% It is used for generating samples when Jocabian matrix is not available
%------------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgo.burnInSamples: There are no enough input arguments \n');
            end

            Swarm = OptAlgo.swarm;

            % Generate random sample with swarm scale
            ParSwarm = rand(Swarm, opt.Nparams+1);

            % Keep all the swarm in the searching domain
            ParSwarm(:, 1:opt.Nparams) = repmat(opt.bounds(1,:), Swarm, 1) + ...
                ParSwarm(:, 1:opt.Nparams) .* repmat( (opt.bounds(2,:) - opt.bounds(1,:)), Swarm, 1 );

            % All the proposals are simulated
            ParSwarm(:, opt.Nparams+1) = arrayfun( @(idx) feval( FUNC, ...
                OptAlgo.pTransfer('exp', ParSwarm(idx, 1:opt.Nparams)) ), 1:Swarm);

            % If the objective value is too big in the case, the proposals are deleted
            ParSwarm(ParSwarm(:, opt.Nparams+1) > opt.rejectValue, :) = [];
            [SS, minRow] = min(ParSwarm(:, opt.Nparams+1));
            oldpar = ParSwarm(minRow, 1:opt.Nparams);

            % Calculate the covariance matrix
            [chaincov, ~, ~] = OptAlgo.covUpdate(ParSwarm(:, 1:opt.Nparams), 1, [], [], []);

            if isempty(chaincov)
                error('OptAlgo.burnInSamples:\n %g randomly generated samples are all rejected with opt.rejectValue = %g \n Please change your rejection value in OptAlgo.getOptionsMCMC function', Swarm, opt.rejectValue);
            end

            % Cholesky decomposition
            R = chol( chaincov + eye(opt.Nparams) * 1e-7 );

        end % burnInSamples

        function [xcov, xmean, wsum] = covUpdate(x, w, oldcov, oldmean, oldwsum)
%------------------------------------------------------------------------------
% Recursive update the covariance matrix
%------------------------------------------------------------------------------


            [n, p] = size(x);

            if n == 0
                xcov = oldcov;
                xmean = oldmean;
                wsum = oldwsum;
                return
            end

            if nargin < 2 || isempty(w)
                w = 1;
            end

            if length(w) == 1
                w = ones(n,1) * w;
            end

            if nargin > 2 && ~isempty(oldcov)

                for i = 1:n
                    xi     = x(i,:);
                    wsum   = w(i);
                    xmeann = xi;
                    xmean  = oldmean + wsum / (wsum + oldwsum) * (xmeann - oldmean);

                    xcov =  oldcov + wsum ./ (wsum + oldwsum - 1) .* (oldwsum / (wsum + oldwsum) ...
                        .* ((xi - oldmean)' * (xi - oldmean)) - oldcov);
                    wsum    = wsum + oldwsum;
                    oldcov  = xcov;
                    oldmean = xmean;
                    oldwsum = wsum;
                end

            else

                wsum  = sum(w);
                xmean = zeros(1,p);
                xcov  = zeros(p,p);

                for i = 1:p
                    xmean(i) = sum(x(:,i) .* w) ./ wsum;
                end

                if wsum > 1
                    for i = 1:p
                        for j = 1:i
                            xcov(i,j) = (x(:,i) - xmean(i))' * ((x(:,j) - xmean(j)) .* w) ./ (wsum - 1);
                            if (i ~= j)
                                xcov(j,i) = xcov(i,j);
                            end
                        end
                    end
                end

            end

        end % covUpdate

        function Population = conversionDataMCMC(maxIter, opt)
%------------------------------------------------------------------------------
% Load and convert the data in ascii format to the mat format
% then generate the population for statistics
%------------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgo.conversionDataMCMC: There are no enough input arguments \n');
            end

            load('chainData.dat');

            % Discard the former 50% chain, retain only the last 50%
            if maxIter < opt.nsamples + opt.burn_in
                idx = floor(0.5 * maxIter);
            else
                idx = floor(0.5 * (opt.nsamples + opt.burn_in));
            end

            if idx < length(chainData), idx = 0; end

            eval(sprintf('chainData(1:idx, :) = [];'));
            Population = chainData;

            save('population.dat', 'Population', '-ascii');

        end % conversionDataMCMC

        function z = Geweke(chain, a, b)
%------------------------------------------------------------------------------
% Geweke's MCMC convergence diagnostic
%------------------------------------------------------------------------------


            if nargin < 3
                a = 0.5;
                if nargin < 2
                    b = 0.8;
                end
            end

            n = length(chain);
            na = floor(a * n);
            nb = floor(b * n);

            if (n-na) / n >= 1
                error('OptAlgo.Geweke: Error with na and nb \n');
            end

            m1 = mean(chain(na:(nb-1),:));
            m2 = mean(chain(nb:end,:));

            % Spectral estimates for variance
            sa = OptAlgo.spectrum0(chain(na:(nb-1),:));
            sb = OptAlgo.spectrum0(chain(nb:end,:));

            z = (m1 - m2) ./ (sqrt( sa / (nb-na) + sb / (n-nb+1)));

        end % Geweke

        function s = spectrum0(x)
%------------------------------------------------------------------------------
% Spectral density at frequency zero
%------------------------------------------------------------------------------


            [m, n] = size(x);
            s = zeros(1,n);

            for i = 1:n
                spec = OptAlgo.spectrum(x(:,i),m);
                s(i) = spec(1);
            end

        end % spectrum0

        function [y, f] = spectrum(x, nfft, nw)
%------------------------------------------------------------------------------
% Power spectral density using Hanning window
%------------------------------------------------------------------------------


            if nargin < 3 || isempty(nw)
                nw = fix(nfft/4);
                if nargin < 2 || isempty(nfft)
                    nfft = min(length(x),256);
                end
            end

            noverlap = fix(nw/2);

            % Hanning window
            w = .5*(1 - cos(2*pi*(1:nw)' / (nw+1)));
            % Daniel
%            w = [0.5;ones(nw-2,1);0.5];
            n = length(x);

            if n < nw
                x(nw) = 0;
                n = nw;
            end

            x = x(:);

            k = fix((n - noverlap) / (nw - noverlap)); % no. of windows
            index = 1:nw;
            kmu = k * norm(w)^2; % Normalizing scale factor
            y = zeros(nfft,1);

            for i = 1:k
%                xw = w.*detrend(x(index),'linear');
                xw = w .* x(index);
                index = index + (nw - noverlap);
                Xx = abs(fft(xw,nfft)).^2;
                y = y + Xx;
            end

            y  = y * (1 / kmu); % normalize

            n2 = floor(nfft / 2);
            y  = y(1:n2);
            f  = 1 ./ n * (0:(n2-1));

        end % spectrum

    end % MCMC


%   Miscellaneous
    methods (Static = true, Access = 'public')

        function y = GammarDistribution(m, n, a, b)
%-----------------------------------------------------------------------------------------
% GammarDistrib random deviates from gamma distribution
%
%  GAMMAR_MT(M,N,A,B) returns a M*N matrix of random deviates from the Gamma
%  distribution with shape parameter A and scale parameter B:
%
%  p(x|A,B) = B^-A/gamma(A)*x^(A-1)*exp(-x/B)
%
%  Uses method of Marsaglia and Tsang (2000)
%
% G. Marsaglia and W. W. Tsang:
% A Simple Method for Generating Gamma Variables,
% ACM Transactions on Mathematical Software, Vol. 26, No. 3, September 2000, 363-372.
%-----------------------------------------------------------------------------------------


            if nargin < 4, b = 1; end

            y = zeros(m, n);
            for j = 1:n
                for i=1: m
                    y(i, j) = Gammar(a, b);
                end
            end

            function y = Gammar(a, b)

                if a < 1

                    y = Gammar(1+a, b) * rand(1) ^ (1/a);

                else

                    d = a - 1/3;
                    c = 1 / sqrt(9*d);

                    while(1)

                        while(1)
                            x = randn(1);
                            v = 1 + c*x;

                            if v > 0, break, end

                        end

                        v = v^3; u = rand(1);

                        if u < 1 - 0.0331*x^4, break, end

                        if log(u) < 0.5 * x^2 + d * (1-v+log(v)), break, end

                    end

                    y = b * d * v;

                end

            end % Gammar

        end % GammerDistribution

        function xLog = pTransfer(str, x)
%------------------------------------------------------------------------------
% If log-scale is enabled
%
% Parameters:
%       - str:
%           + 'exp'. X = exp(x)
%           + 'log'. X = log(x)
%       - x. Parameter set
%------------------------------------------------------------------------------


            if OptAlgo.logScale
                if strcmp('exp', str)
                    xLog = exp(x);
                elseif strcmp('log', str)
                    xLog = log(x);
                else
                    error('OptAlgo.LogTransfer \n');
                end
            else
                xLog = x;
            end

        end % logTransfer

        function prior = priorPDF(points)
%------------------------------------------------------------------------------
% This is a routine used for constructe the prior distribution for MCMC
%
% Parameters:
%       points. The estimated parameters
%
% Return:
%       prior. The prior possibility, prior = f(point).
%------------------------------------------------------------------------------


            if isempty(OptAlgo.prior)
                prior = 1;
                return;
            end

            % Load prior data file
            if OptAlgo.logScale
                data = OptAlgo.pTransfer('exp', OptAlgo.prior(:, 1:end-1));
            else
                data = OptAlgo.prior(:, 1:end-1);
            end
            [~, d] = size(data);

            % Get mesh grid and pdf of the prior
            if exist('tmp.mat', 'file') ~= 2
                OptAlgo.multivariatePrior(data, d);
            end

            proposedPoint = cell(1, d);
            for i = 1:d
                if OptAlgo.logScale
                    proposedPoint{i} = OptAlgo.pTransfer('exp', points(i));
                else
                    proposedPoint{i} = points(i);
                end
            end

            load('tmp.mat');
            % Evaluate the density value of the new proposal
            prior = interpn(fullAxisMesh{:}, pdfVal, proposedPoint{:});
            if isnan(prior), prior = 1; end

        end % priorPDF

        function multivariatePrior(data, d)
%------------------------------------------------------------------------------
% Build mesh grid and invoke mvkde (multivariate kernel density estimator) routine
%------------------------------------------------------------------------------


            % Preallocation
            temp = [];
            meshSize = 20;
            axisMesh = cell(1, d);
            fullAxisMesh = cell(1, d);

            maxVector = max(data, [], 1);
            minVector = min(data, [], 1);
            rangeVec  = maxVector - minVector;

            % Axial mesh
            for i = 1:d
                axisMesh{i} = minVector(i) : rangeVec(i)/(meshSize-1) : maxVector(i);
            end

            % Generate n-demensional grid
            [fullAxisMesh{:}] = ndgrid(axisMesh{:});

            % Reshape grid to call mvkde routine
            for i = 1:d
                temp = [temp, fullAxisMesh{i}(:)];
            end
            grids = reshape(temp, meshSize^d, d);

            % Invoke multivariate kernel density estimator
            pdfVal = OptAlgo.mvkde(data, grids);
            % Inverse reshape
            pdfVal = reshape(pdfVal, size(fullAxisMesh{1}));

            % Store the mesh and pdf information for further use
            save('tmp.mat', 'axisMesh', 'fullAxisMesh', 'pdfVal');

        end % multivariatePrior

        function FigurePlot(Population, opt)
%------------------------------------------------------------------------------
% Plot the histgram and scatter figures of the last population
%------------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgo.FigurePlot: There are no enough input arguments \n');
            end

            figure(1);clf
            for i = 1: opt.Nparams

                subplot(round(opt.Nparams/2),2,i,'Parent', figure(1));

%                histfit(OptAlgo.pTransfer('exp', Population(:,i)), 50, 'kernel');
                [counts, centers] = hist( Population(:,i), 25 );
                bar(centers, counts, 'r');

                xlabel(sprintf('$\\ln(x_%d)$', i), 'FontSize', 16, 'Interpreter', 'latex');
                ylabel(sprintf('Frequency'), 'FontSize', 14, 'Interpreter', 'latex');
                set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
%                OptAlgo.tickLabelFormat(gca, 'x', '%0.2e');
%                OptAlgo.tickLabelFormat(gca, 'x', []);
%                set(gca, 'XTickLabel', num2str(get(gca, 'xTick')', '%g'));
%                OptAlgo.xtickLabelRotate([], 15, [], 'FontSize', 20, 'FontName', 'Times New Roman');
                set(gca, 'ygrid', 'on');

            end

            figure(2);clf
            for i = 1: opt.Nparams-1
                for j = i: opt.Nparams-1

                    subplot(opt.Nparams-1, opt.Nparams-1, j+(i-1)*(opt.Nparams-1), 'Parent', figure(2));

                    scatter(Population(:,j+1), Population(:,i), 8, 'o', 'MarkerEdgeColor', [0, 0.7, 0.7]);

                    xlabel(sprintf('$\\ln(x_%d)$', j+1), 'FontName', 'Times New Roman', 'FontSize', 16, 'Interpreter', 'latex');
                    ylabel(sprintf('$\\ln(x_%d)$', i), 'FontName', 'Times New Roman', 'FontSize', 16, 'Interpreter', 'latex');
                    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
%                    OptAlgo.tickLabelFormat(gca, 'x', '%0.2e');
%                    set(gca, 'XTickLabel', num2str(get(gca, 'xTick')', '%g'));
%                    OptAlgo.xtickLabelRotate([], 15, [], 'FontSize', 20, 'FontName', 'Times New Roman');
                    grid on;

                end
            end

            figure(3);clf
            for i = 1: opt.Nparams

                subplot(round(opt.Nparams/2),2,i,'Parent', figure(3));

%                [y, x] = ksdensity(OptAlgo.pTransfer('exp', Population(:,i))); area(x, y, 'FaceColor', 'g');
                [y, x] = ksdensity( Population(:,i));
                area(x, y, 'FaceColor', 'g');

                xlabel(sprintf('$\\ln(x_%d)$', i), 'FontSize', 16, 'Interpreter', 'latex');
                ylabel(sprintf('$p(\\theta)$'), 'FontSize', 16, 'Interpreter', 'latex');
                set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
                set(gca, 'xgrid', 'on');
                set(gca, 'ygrid', 'on');

            end

        end % FigurePlot

        function [pdf, X1, X2] = mvkde(X, grid, gam)
%------------------------------------------------------------------------------
% adaptive kernel density estimation in high dimensions;
%
% INPUTS:   X  - data as a 'n' by 'd' vector;
%
%         grid - 'm' points of dimension 'd' over which pdf is computed;
%                default provided only for 2-dimensional data;
%                see example on how to construct it in higher dimensions
%
%          gam - cost/accuracy tradeoff parameter, where gam<n;
%                default value is gam=ceil(n^(1/2)); larger values
%                may result in better accuracy, but always reduce speed;
%                to speedup the code, reduce the value of "gam";
%
% OUTPUT: pdf   - the value of the estimated density at 'grid'
%         X1,X2 - grid only for 2 dimensional data
%
%
%  Reference:
%  Kernel density estimation via diffusion
%  Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
%  Annals of Statistics, Volume 38, Number 5, pages 2916-2957.
%------------------------------------------------------------------------------


            [n, d] = size(X);

            % begin scaling preprocessing
            MAX = max(X, [], 1);
            MIN = min(X, [], 1);
            scaling = MAX - MIN;

            MAX = MAX + scaling/10;
            MIN = MIN - scaling/10;
            scaling = MAX - MIN;

            X = bsxfun(@minus, X, MIN);
            X = bsxfun(@rdivide, X, scaling);

            % failing to provide grid
            if (nargin < 2) || isempty(grid)

                warning('Assuming data is 2 dimensional. For higher dimensions, provide a grid as in example.')

                % create meshgrid in 2-dimensions
                [X1, X2] = meshgrid( MIN(1):scaling(1)/(2^7-1):MAX(1), MIN(2):scaling(2)/(2^7-1):MAX(2) );

                % create grid for plotting
                grid=reshape([X1(:),X2(:)], 2^14, d);

            end

            mesh = bsxfun(@minus, grid, MIN);
            mesh = bsxfun(@rdivide, mesh, scaling);

            % failing to provide speed/accuracy tradeoff
            if nargin < 3
                gam = ceil(n^(1/2));
            end

            % algorithm initialization
            del = 0.1 / n^(d/(d+4));
            perm = randperm(n);
            mu = X(perm(1:gam), :);
            w = rand(1, gam);
            w = w / sum(w);
            Sig = bsxfun(@times, rand(d,d,gam), eye(d)*del);
            ent = -Inf;

            % begin algorithm
            for iter = 1:1500
                Eold = ent;

                % update parameters
                [w, mu, Sig, del, ent] = OptAlgo.regEM(w, mu, Sig, del, X);

                % stopping condition
                err = abs( (ent-Eold) / ent );
                if (err < 10^(-4)) || iter > 200, break, end
            end

            % now output density values at grid
            pdf = OptAlgo.probfun(mesh, w, mu, Sig) / prod(scaling); % evaluate density

            % adjust bandwidth for scaling
            del = del * scaling;

        end % mvkde

        function pdf = probfun(x, w, mu, Sig)
%------------------------------------------------------------------------------
%
%------------------------------------------------------------------------------


            [gam, d] = size(mu);

            pdf = 0;

            for k = 1:gam

                L = chol(Sig(:, :, k));
                s = diag(L);

                logpdf = -0.5 * sum( (bsxfun(@minus, x, mu(k, :)) / L).^2, 2 ) + log(w(k)) - ...
                    sum(log(s)) - d * log(2*pi) / 2;

                pdf = pdf + exp(logpdf);

            end

        end % probfun

        function [w, mu, Sig, del, ent] = regEM(w, mu, Sig, del, X)
%------------------------------------------------------------------------------
%
%------------------------------------------------------------------------------


            [gam, d] = size(mu);
            [n, d] = size(X);

            log_lh = zeros(n, gam);
            log_sig = log_lh;

            for i = 1:gam

                L = chol(Sig(:, :, i));

                Xcentered = bsxfun(@minus, X, mu(i,:));

                xRinv = Xcentered / L; xSig = sum((xRinv / L').^2,2) + eps;

                log_lh(:, i) =-0.5 * sum(xRinv.^2, 2) - sum(log(diag(L))) + ...
                    log(w(i)) - d * log(2*pi) / 2 - 0.5 * del^2 * trace((eye(d)/L)/L');

                log_sig(:, i) = log_lh(:, i) + log(xSig);

            end

            maxll = max (log_lh, [], 2);
            maxlsig = max (log_sig, [], 2);

            p = exp(bsxfun(@minus, log_lh, maxll));
            psig = exp(bsxfun(@minus, log_sig, maxlsig));

            density = sum(p, 2);
            psigd = sum(psig, 2);

            logpdf = log(density) + maxll;
            logpsigd = log(psigd) + maxlsig;

            p = bsxfun(@rdivide, p, density);

            ent = sum(logpdf);

            w = sum(p, 1);

            for i = find(w > 0)
                %compute mu's
                mu(i, :) = p(:, i)' * X / w(i);

                Xcentered = bsxfun(@minus, X, mu(i,:));
                Xcentered = bsxfun(@times, sqrt(p(:, i)), Xcentered);

                % compute sigmas
                Sig(:, :, i) = Xcentered' * Xcentered / w(i) + del^2 * eye(d);
            end

            % estimate curvature
            w = w / sum(w);

            curv = mean( exp(logpsigd - logpdf) );

            del = 1 / (4 * n * (4*pi)^(d/2) *curv)^(1 / (d + 2));

        end % regEM

        function tickLabelFormat(hAxes, axName, format)
%------------------------------------------------------------------------------
% Sets the format of the tick labels
%
% Syntax:
%    ticklabelformat(hAxes, axName, format)
%
% Input Parameters:
%    hAxes  - handle to the modified axes, such as returned by the gca function
%    axName - name(s) of axles to modify: 'x','y','z' or combination (e.g. 'xy')
%    format - format of the tick labels in sprintf format (e.g. '%.1f V') or a
%             function handle that will be called whenever labels need to be updated
%
%    Note: Calling TICKLABELFORMAT again with an empty ([] or '') format will revert
%    ^^^^  to Matlab's normal tick labels display behavior
%
% Examples:
%    ticklabelformat(gca,'y','%.6g V') - sets y axis on current axes to display 6 significant digits
%    ticklabelformat(gca,'xy','%.2f')  - sets x & y axes on current axes to display 2 decimal digits
%    ticklabelformat(gca,'z',@myCbFcn) - sets a function to update the Z tick labels on current axes
%    ticklabelformat(gca,'z',{@myCbFcn,extraData}) - sets an update function as above, with extra data
%
% Warning:
%    This code heavily relies on undocumented and unsupported Matlab functionality.
%    It works on Matlab 7+, but use at your own risk!
%
% Technical description and more details:
%    http://UndocumentedMatlab.com/blog/setting-axes-tick-labels-format
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany@gmail.com)
%------------------------------------------------------------------------------


            % Check # of args (we now have narginchk but this is not available on older Matlab releases)
            if nargin < 1
                help(mfilename)
                return;
            elseif nargin < 3
                error('OptAlgo.ticklabelformat: Not enough input arguments \n');
            end

            % Check input args
            if ~ishandle(hAxes) || (~isa(handle(hAxes),'axes') && ~isa(handle(hAxes),'matlab.graphics.axis.Axes'))
                error('OptAlgo.ticklabelformat: hAxes input argument must be a valid axes handle \n');
            elseif ~ischar(axName)
                error('OptAlgo.ticklabelformat: axName input argument must be a string \n');
            elseif ~isempty(format) && ~ischar(format) && ~isa(format,'function_handle') && ~iscell(format)
                error('OptAlgo.ticklabelformat: format input argument must be a string or function handle \n');
            end

            % normalize axes name(s) to lowercase
            axName = lower(axName);

            if strfind(axName,'x')
                install_adjust_ticklbl(hAxes,'X',format)
            elseif strfind(axName,'y')
                install_adjust_ticklbl(hAxes,'Y',format)
            elseif strfind(axName,'z')
                install_adjust_ticklbl(hAxes,'Z',format)
            end

            function install_adjust_ticklbl(hAxes,axName,format)
            % Install the new tick labels for the specified axes

                % If empty format was specified
                if isempty(format)

                    % Remove the current format (revert to default Matlab format)
                    set(hAxes,[axName 'TickLabelMode'],'auto')
                    setappdata(hAxes,[axName 'TickListener'],[])
                    return

                end

                % Determine whether to use the specified format as a
                % sprintf format or a user-specified callback function
                if ischar(format)
                    cb = {@adjust_ticklbl axName format};
                else
                    cb = format;
                end

                % Now install axis tick listeners to adjust tick labels
                % (use undocumented feature for adjustments)
                ha = handle(hAxes);
                propName = [axName 'Tick'];
                hp = findprop(ha,propName);

                try
                    % R2014a or older
                    hl = handle.listener(ha,hp,'PropertyPostSet',cb);

                    % Adjust tick labels now
                    % eventData.AffectedObject = hAxes;
                    % adjust_ticklbl([],eventData,axName,format)
                    set(hAxes,[axName 'TickLabelMode'],'manual')
                    set(hAxes,[axName 'TickLabelMode'],'auto')
                catch
                    % R2014b or newer
                    if iscell(cb)
                        cb = @(h,e) feval(cb{1},h,e,cb{2:end});
                    end
                    % hl(1) = addlistener(ha,propName,'PostSet',cb);
                    hl(1) = event.proplistener(ha,hp,'PostSet',cb);

                    % *Tick properties don't trigger PostSet events when updated automatically in R2014b - need to use *Lim
                    % addlistener(ha,[axName 'Lim'],'PostSet',cb);
                    hRuler = get(ha,[axName 'Ruler']);
                    % hl(2) = addlistener(hRuler,'MarkedClean',cb);
                    hl(2) = event.listener(hRuler,'MarkedClean',cb);

                    % Adjust tick labels now
                    eventData.AffectedObject = hAxes;
                    if ischar(format)
                        adjust_ticklbl([],eventData,axName,format)
                    else
                        hgfeval(format);
                    end
                    set(hAxes,[axName 'TickLabelMode'],'manual')
                    % set(hAxes,[axName 'TickLabelMode'],'auto')  % causes labels not to be updated in R2014b!
                end

                setappdata(hAxes,[axName 'TickListener'],hl)
                % drawnow;

                function adjust_ticklbl(hProp,eventData,axName,format)
                % Default tick labels update callback function (used if user did not specify their own function)

                    try
                        hAxes = eventData.AffectedObject;
                    catch
                        hAxes = ancestor(eventData.Source,'Axes');
                    end
                    tickValues = get(hAxes,[axName 'Tick']);
                    tickLabels = arrayfun(@(x)(sprintf(format,x)),tickValues,'UniformOutput',false);
                    set(hAxes,[axName 'TickLabel'],tickLabels);

                end

            end % install_adjust_ticklbl

        end % tickLabelFormat

        function hText = xtickLabelRotate(XTick, rot, varargin)
%------------------------------------------------------------------------------
% xtick label rotate
%
% Parameter:
%       - XTick. vector array of XTick positions & values (numeric) uses current
%           XTick values or XTickLabel cell array by default (if empty)
%       - rot. angle of rotation in degrees, 90° by default
%       - XTickLabel: cell array of label strings
%       - [var]. Optional. "Property-value" pairs passed to text generator
%           ex: 'interpreter','none', 'Color','m','Fontweight','bold'
%
% Return:
%       - hText. handle vector to text labels
%
% Example 1:  Rotate existing XTickLabels at their current position by 90°
%    xticklabel_rotate
%
% Example 2:  Rotate existing XTickLabels at their current position by 45° and change
% font size
%    xticklabel_rotate([],45,[],'FontSize',14)
%
% Example 3:  Set the positions of the XTicks and rotate them 90°
%    figure;  plot([1960:2004],randn(45,1)); xlim([1960 2004]);
%    xticklabel_rotate([1960:2:2004]);
%
% Example 4:  Use text labels at XTick positions rotated 45° without tex interpreter
%    xticklabel_rotate(XTick,45,NameFields,'interpreter','none');
%
% Example 5:  Use text labels rotated 90° at current positions
%    xticklabel_rotate([],90,NameFields);
%
% Example 6:  Multiline labels
%    figure;plot([1:4],[1:4])
%    axis([0.5 4.5 1 4])
%    xticklabel_rotate([1:4],45,{{'aaa' 'AA'};{'bbb' 'AA'};{'ccc' 'BB'};{'ddd' 'BB'}})
%
% This is a modified version of xticklabel_rotate90 by Denis Gilbert
% Modifications include Text labels (in the form of cell array)
%       Arbitrary angle rotation
%       Output of text handles
%       Resizing of axes and title/xlabel/ylabel positions to maintain same overall size
%          and keep text on plot
%          (handles small window resizing after, but not well due to proportional placement with
%           fixed font size. To fix this would require a serious resize function)
%       Uses current XTick by default
%       Uses current XTickLabel is different from XTick values (meaning has been already defined)
%
% Author: Brian FG Katz, bfgkatz@hotmail.com
%------------------------------------------------------------------------------


            % check to see if xticklabel_rotate has already been here (no other reason for this to happen)
            if isempty(get(gca,'XTickLabel'))
                error('OptAlgo.xtickLabelRotate: can not process, either xticklabel_rotate has already been run or XTickLabel field has been erased \n');
            end

            % Modified with forum comment by "Nathan Pust" allow the current text labels to be used and property value pairs to be changed for those labels
            if (nargin < 3 || isempty(varargin{1})) && (~exist('XTick') || isempty(XTick))

                xTickLabels = get(gca,'XTickLabel');

                if ~iscell(xTickLabels)
                    % remove trailing spaces if exist (typical with auto generated XTickLabel)
                    temp1 = num2cell(xTickLabels,2);
                    for loop = 1:length(temp1),
                        temp1{loop} = deblank(temp1{loop});
                    end
                    xTickLabels = temp1;
                end

                varargin = varargin(2:length(varargin));

            end

            % if no XTick is defined use the current XTick
            if (~exist('XTick') || isempty(XTick))
                XTick = get(gca,'XTick');
            end

            % Make XTick a column vector
            XTick = XTick(:);

            if ~exist('xTickLabels')
                % Define the xtickLabels
                % If XtickLabel is passed as a cell array then use the text
                if ~isempty(varargin) && (iscell(varargin{1}))
                    xTickLabels = varargin{1};
                    varargin = varargin(2:length(varargin));
                else
                    xTickLabels = num2str(XTick);
                end
            end

            if length(XTick) ~= length(xTickLabels)
                error('OptAlgo.xtickLabelRotate: must have same number of elements in "XTick" and "XTickLabel" \n');
            end

            % Set the Xtick locations and set XTicklabel to an empty string
            set(gca,'XTick',XTick,'XTickLabel','');

            if nargin < 2
                rot = 90 ;
            end

            % Determine the location of the labels based on the position of the xlabel
            hxLabel = get(gca,'XLabel');
            xLabelString = get(hxLabel,'String');

            set(hxLabel,'Units','data');
            xLabelPosition = get(hxLabel,'Position');
            y = xLabelPosition(2);

            % CODE below was modified following suggestions from Urs Schwarz
            y = repmat(y,size(XTick,1),1);
            % retrieve current axis' fontsize
            fs = get(gca,'fontsize');

            if ~iscell(xTickLabels)
                % Place the new xTickLabels by creating TEXT objects
                hText = text(XTick, y, xTickLabels,'fontsize',fs);
            else
                % Place multi-line text approximately where tick labels belong
                for cnt=1:length(XTick)
                    hText(cnt) = text(XTick(cnt),y(cnt),xTickLabels{cnt}, ...
                        'VerticalAlignment','top', 'UserData','xtick');
                end
            end

            % Rotate the text objects by ROT degrees
            % Modified with modified forum comment by "Korey Y" to deal with labels at top
            % Further edits added for axis position
            xAxisLocation = get(gca, 'XAxisLocation');
            if strcmp(xAxisLocation,'bottom')
                set(hText,'Rotation',rot,'HorizontalAlignment','right',varargin{:});
            else
                set(hText,'Rotation',rot,'HorizontalAlignment','left',varargin{:});
            end

            % Adjust the size of the axis to accomodate for longest label (like if they are text ones)
            % This approach keeps the top of the graph at the same place and tries to keep xlabel at the same place
            % This approach keeps the right side of the graph at the same place

            set(get(gca,'xlabel'),'units','data');
                labxorigpos_data = get(get(gca,'xlabel'),'position');
            set(get(gca,'ylabel'),'units','data');
                labyorigpos_data = get(get(gca,'ylabel'),'position');
            set(get(gca,'title'),'units','data');
                labtorigpos_data = get(get(gca,'title'),'position');

            set(gca,'units','pixel');
            set(hText,'units','pixel');
            set(get(gca,'xlabel'),'units','pixel');
            set(get(gca,'ylabel'),'units','pixel');

            origpos = get(gca,'position');

            % Modified with forum comment from "Peter Pan" to deal with case when only one XTickLabelName is given.
            x = get( hText, 'extent' );
            if iscell( x ) == true
                textsizes = cell2mat( x );
            else
                textsizes = x;
            end

            largest =  max(textsizes(:,3));
            longest =  max(textsizes(:,4));

            laborigext = get(get(gca,'xlabel'),'extent');
            laborigpos = get(get(gca,'xlabel'),'position');

            labyorigext = get(get(gca,'ylabel'),'extent');
            labyorigpos = get(get(gca,'ylabel'),'position');
            leftlabdist = labyorigpos(1) + labyorigext(1);

            % assume first entry is the farthest left
            leftpos = get(hText(1),'position');
            leftext = get(hText(1),'extent');
            leftdist = leftpos(1) + leftext(1);
            if leftdist > 0, leftdist = 0; end

            % Modified to allow for top axis labels and to minimize axis resizing
            if strcmp(xAxisLocation,'bottom')
                newpos = [origpos(1)-(min(leftdist,labyorigpos(1)))+labyorigpos(1) ...
                    origpos(2)+((longest+laborigpos(2))-get(gca,'FontSize')) ...
                    origpos(3)-(min(leftdist,labyorigpos(1)))+labyorigpos(1)-largest ...
                    origpos(4)-((longest+laborigpos(2))-get(gca,'FontSize'))];
            else
                newpos = [origpos(1)-(min(leftdist,labyorigpos(1)))+labyorigpos(1) ...
                    origpos(2) ...
                    origpos(3)-(min(leftdist,labyorigpos(1)))+labyorigpos(1)-largest ...
                    origpos(4)-(longest)+get(gca,'FontSize')];
            end
            set(gca,'position',newpos);

            % readjust position of text labels after resize of plot
            set(hText,'units','data');
            for loop= 1:length(hText)
                set(hText(loop),'position',[XTick(loop), y(loop)]);
            end

            % adjust position of xlabel and ylabel
            laborigpos = get(get(gca,'xlabel'),'position');
            set(get(gca,'xlabel'),'position',[laborigpos(1) laborigpos(2)-longest 0]);

            % switch to data coord and fix it all
            set(get(gca,'ylabel'),'units','data');
            set(get(gca,'ylabel'),'position',labyorigpos_data);
            set(get(gca,'title'),'position',labtorigpos_data);

            set(get(gca,'xlabel'),'units','data');
                labxorigpos_data_new = get(get(gca,'xlabel'),'position');
            set(get(gca,'xlabel'),'position',[labxorigpos_data(1) labxorigpos_data_new(2)]);

            % Reset all units to normalized to allow future resizing
            set(get(gca,'xlabel'),'units','normalized');
            set(get(gca,'ylabel'),'units','normalized');
            set(get(gca,'title'),'units','normalized');
            set(hText,'units','normalized');
            set(gca,'units','normalized');

            if nargout < 1
                clear hText
            end

        end % xtickLabelRotate

        function res = crashSaver(x, e)

            save('crashParams.dat', 'x', '-ascii', '-append');
            fprintf('%s The parameter set results in this crash has been stored in crashParams.dat\n', e.message);
            res = 68106800;

        end % crashSaver

    end % misc

end % classdef
% =============================================================================
%              The MATLAB library for optimization case studies
%
%      Copyright © 2015-2017: Qiaole He
%
%      Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
