function fitLysozyme_AIO
%==============================================================================
% Fit the breakthrough curves of lysozyme
% In the estimation of parameters due to the datasets (e.g., different elution
% gradient) from one experiment, the effcts of datasets are lumped together to
% determine parameter distributions
%
% fitting parameters: theta = {D_ax, K_f, D_p, e_c, e_p, k_eq, nu, sigma}
%==============================================================================


    % Create fit object
    pf = ParameterFit();

% Experiment one
%------------------------------------------------------------------------------
    % Collect data of experiment in cell array (one cell per observed compoenent)
    % Note that time points of the measurements are given in sim.solutionTimes
    temp = load('data_lyz_5cv.dat'); % data = [time, conc]
    data = cell(1,1);
    data{1} = temp(:,2);

    sim = createModel(temp(:,1), '5');

    % Set model parameters and enable sensitivities
    params = cell(8, 1);
    params{1} = makeSensitivity([0], {'COL_DISPERSION'},[-1], [-1], [-1], [-1]);
    params{2} = makeSensitivity([0], {'FILM_DIFFUSION'}, [0], [-1], [-1], [-1]);
    params{3} = makeSensitivity([0], {'PAR_DIFFUSION'},  [0], [-1], [-1], [-1]);
    params{4} = makeSensitivity([0], {'COL_POROSITY'},  [-1], [-1], [-1], [-1]);
    params{5} = makeSensitivity([0], {'PAR_POROSITY'},  [-1], [-1], [-1], [-1]);
    params{6} = makeSensitivity([0], {'SMA_KA'},         [1], [-1], [0], [-1]);
    params{7} = makeSensitivity([0], {'SMA_NU'},         [1], [-1], [0], [-1]);
    params{8} = makeSensitivity([0], {'SMA_SIGMA'},      [1], [-1], [0], [-1]);
    sim.setParameters(params, false(8, 1));

    % Specify which components are observed in which factor for each observation / wavelength
    idxComp = cell(1,1);
    idxComp{1} = [0, 1];

    % Add the experiment to the fit
    pf.addExperiment(data, sim, [0], idxComp, [], [], [], [], 'Lysozyme', {'5cv'});

    % Enable logarithmic parameter transformation for all parameters
    pf.parameterTransform = LogParameterTransformation(ones(1, 8), true(1, 8));


% Experiment two
%------------------------------------------------------------------------------
    % Collect data of experiment in cell array (one cell per observed compoenent)
    % Note that time points of the measurements are given in sim.solutionTimes
    temp = load('data_lyz_10cv.dat'); % data = [time, conc]
    data = cell(1,1);
    data{1} = temp(:,2);

    sim = createModel(temp(:,1), '10');

    % Set model parameters and enable sensitivities
    params = cell(8, 1);
    params{1} = makeSensitivity([0], {'COL_DISPERSION'},[-1], [-1], [-1], [-1]);
    params{2} = makeSensitivity([0], {'FILM_DIFFUSION'}, [0], [-1], [-1], [-1]);
    params{3} = makeSensitivity([0], {'PAR_DIFFUSION'},  [0], [-1], [-1], [-1]);
    params{4} = makeSensitivity([0], {'COL_POROSITY'},  [-1], [-1], [-1], [-1]);
    params{5} = makeSensitivity([0], {'PAR_POROSITY'},  [-1], [-1], [-1], [-1]);
    params{6} = makeSensitivity([0], {'SMA_KA'},         [1], [-1], [0], [-1]);
    params{7} = makeSensitivity([0], {'SMA_NU'},         [1], [-1], [0], [-1]);
    params{8} = makeSensitivity([0], {'SMA_SIGMA'},      [1], [-1], [0], [-1]);
    sim.setParameters(params, false(8, 1));

    % Specify which components are observed in which factor for each observation / wavelength
    idxComp = cell(1,1);
    idxComp{1} = [0, 1];

    % Add the experiment to the fit
    pf.addExperiment(data, sim, [0], idxComp, [], [], [], [], 'Lysozyme', {'10cv'});

    % Enable logarithmic parameter transformation for all parameters
    pf.parameterTransform = LogParameterTransformation(ones(1, 8), true(1, 8));


% Experiment three
%------------------------------------------------------------------------------
    % Collect data of experiment in cell array (one cell per observed compoenent)
    % Note that time points of the measurements are given in sim.solutionTimes
    temp = load('data_lyz_30cv.dat'); % data = [time, conc]
    data = cell(1,1);
    data{1} = temp(:,2);

    sim = createModel(temp(:,1), '30');

    % Set model parameters and enable sensitivities
    params = cell(8, 1);
    params{1} = makeSensitivity([0], {'COL_DISPERSION'},[-1], [-1], [-1], [-1]);
    params{2} = makeSensitivity([0], {'FILM_DIFFUSION'}, [0], [-1], [-1], [-1]);
    params{3} = makeSensitivity([0], {'PAR_DIFFUSION'},  [0], [-1], [-1], [-1]);
    params{4} = makeSensitivity([0], {'COL_POROSITY'},  [-1], [-1], [-1], [-1]);
    params{5} = makeSensitivity([0], {'PAR_POROSITY'},  [-1], [-1], [-1], [-1]);
    params{6} = makeSensitivity([0], {'SMA_KA'},         [1], [-1], [0], [-1]);
    params{7} = makeSensitivity([0], {'SMA_NU'},         [1], [-1], [0], [-1]);
    params{8} = makeSensitivity([0], {'SMA_SIGMA'},      [1], [-1], [0], [-1]);
    sim.setParameters(params, false(8, 1));

    % Specify which components are observed in which factor for each observation / wavelength
    idxComp = cell(1,1);
    idxComp{1} = [0, 1];

    % Add the experiment to the fit
    pf.addExperiment(data, sim, [0], idxComp, [], [], [], [], 'Lysozyme', {'30cv'});

    % Enable logarithmic parameter transformation for all parameters
    pf.parameterTransform = LogParameterTransformation(ones(1, 8), true(1, 8));


% Experiment four
%------------------------------------------------------------------------------
    % Collect data of experiment in cell array (one cell per observed compoenent)
    % Note that time points of the measurements are given in sim.solutionTimes
    temp = load('data_lyz_60cv.dat'); % data = [time, conc]
    data = cell(1,1);
    data{1} = temp(:,2);

    sim = createModel(temp(:,1), '60');

    % Set model parameters and enable sensitivities
    params = cell(8, 1);
    params{1} = makeSensitivity([0], {'COL_DISPERSION'},[-1], [-1], [-1], [-1]);
    params{2} = makeSensitivity([0], {'FILM_DIFFUSION'}, [0], [-1], [-1], [-1]);
    params{3} = makeSensitivity([0], {'PAR_DIFFUSION'},  [0], [-1], [-1], [-1]);
    params{4} = makeSensitivity([0], {'COL_POROSITY'},  [-1], [-1], [-1], [-1]);
    params{5} = makeSensitivity([0], {'PAR_POROSITY'},  [-1], [-1], [-1], [-1]);
    params{6} = makeSensitivity([0], {'SMA_KA'},         [1], [-1], [0], [-1]);
    params{7} = makeSensitivity([0], {'SMA_NU'},         [1], [-1], [0], [-1]);
    params{8} = makeSensitivity([0], {'SMA_SIGMA'},      [1], [-1], [0], [-1]);
    sim.setParameters(params, false(8, 1));

    % Specify which components are observed in which factor for each observation / wavelength
    idxComp = cell(1,1);
    idxComp{1} = [0, 1];

    % Add the experiment to the fit
    pf.addExperiment(data, sim, [0], idxComp, [], [], [], [], 'Lysozyme', {'60cv'});

    % Enable logarithmic parameter transformation for all parameters
    pf.parameterTransform = LogParameterTransformation(ones(1, 8), true(1, 8));


% Experiment five
%------------------------------------------------------------------------------
    % Collect data of experiment in cell array (one cell per observed compoenent)
    % Note that time points of the measurements are given in sim.solutionTimes
    temp = load('data_lyz_120cv.dat'); % data = [time, conc]
    data = cell(1,1);
    data{1} = temp(:,2);

    sim = createModel(temp(:,1), '120');

    % Set model parameters and enable sensitivities
    params = cell(8, 1);
    params{1} = makeSensitivity([0], {'COL_DISPERSION'},[-1], [-1], [-1], [-1]);
    params{2} = makeSensitivity([0], {'FILM_DIFFUSION'}, [0], [-1], [-1], [-1]);
    params{3} = makeSensitivity([0], {'PAR_DIFFUSION'},  [0], [-1], [-1], [-1]);
    params{4} = makeSensitivity([0], {'COL_POROSITY'},  [-1], [-1], [-1], [-1]);
    params{5} = makeSensitivity([0], {'PAR_POROSITY'},  [-1], [-1], [-1], [-1]);
    params{6} = makeSensitivity([0], {'SMA_KA'},         [1], [-1], [0], [-1]);
    params{7} = makeSensitivity([0], {'SMA_NU'},         [1], [-1], [0], [-1]);
    params{8} = makeSensitivity([0], {'SMA_SIGMA'},      [1], [-1], [0], [-1]);
    sim.setParameters(params, false(8, 1));

    % Specify which components are observed in which factor for each observation / wavelength
    idxComp = cell(1,1);
    idxComp{1} = [0, 1];

    % Add the experiment to the fit
    pf.addExperiment(data, sim, [0], idxComp, [], [], [], [], 'Lysozyme', {'120cv'});

    % Enable logarithmic parameter transformation for all parameters
    pf.parameterTransform = LogParameterTransformation(ones(1, 8), true(1, 8));

%------------------------------------------------------------------------------

    % The intent here is that I am trying to find a parameter set, which fits
    % all the datasets (e.g., varying elution gradients) here. Thus I link all
    % the parameters in the fashion of ([indexExp], [indexParam])
    pf.addLink([1, 2, 3, 4, 5], [1, 1, 1, 1, 1]);
    pf.addLink([1, 2, 3, 4, 5], [2, 2, 2, 2, 2]);
    pf.addLink([1, 2, 3, 4, 5], [3, 3, 3, 3, 3]);
    pf.addLink([1, 2, 3, 4, 5], [4, 4, 4, 4, 4]);
    pf.addLink([1, 2, 3, 4, 5], [5, 5, 5, 5, 5]);
    pf.addLink([1, 2, 3, 4, 5], [6, 6, 6, 6, 6]);
    pf.addLink([1, 2, 3, 4, 5], [7, 7, 7, 7, 7]);
    pf.addLink([1, 2, 3, 4, 5], [8, 8, 8, 8, 8]);

    % This variable serves as storage for the plot handles returned by the plot function of parameterFit
    plotHd = [];

    opt.params = params;
    % Specify the searching domain boundary for the algorithm
    opt.paramBound = [2.29e-11, 1.38e-10; 8.29e-6, 3.04e-5; 1.7e-10, 5.08e-10; 0.23, 0.245; 0.65, 0.71; 1e-1 1; 1 10; 10 35];

    % Call algorithms in OptAlgorithms to implement parameter estimation
    OptAlgo.Markov_Chain_Monte_Carlo(@residualSumOfSquares, opt);


    function res = residualSumOfSquares(x)

        % Calculate the residual between meas data and simulated data
        res = pf.residualVector(x);
        % Calculate the sum of square: res = R' * R
        res = sum(res.^2);
        % Visualization
        plotHd = pf.plot(plotHd, 0, [], false, false, false);

    end

end

function sim = createModel(tOut, str)
%------------------------------------------------------------------------------
% Create the model
%------------------------------------------------------------------------------


    mGrm = SingleGRM();

    % components
    mGrm.nComponents = 2;
    % Discretization
    mGrm.nCellsColumn = 40;
    mGrm.nCellsParticle = 3;
    mGrm.nBoundStates = ones(mGrm.nComponents, 1);

    % Initial conditions
    mGrm.initialBulk = [20.0, 0.0];
    mGrm.initialSolid = [962, 0.0];

    % Transport
    mGrm.dispersionColumn          = 3.3746e-8;
    mGrm.filmDiffusion             = [1.3209e-5, 1.3209e-5];
    mGrm.diffusionParticle         = [6.009e-10, 6.009e-10];
    mGrm.diffusionParticleSurface  = [0.0, 0.0];

    % Geometry
    mGrm.columnLength      = 0.025;
    mGrm.particleRadius    = 4.5e-5;
    mGrm.porosityColumn    = 0.242;
    mGrm.porosityParticle  = 0.6877;

    volumeFlow = 0.5e-6/60; % 0.5 ml/min --> l/s
    columnRadius = 0.35e-2; % m
    mGrm.interstitialVelocity = volumeFlow / (pi * columnRadius^2 * mGrm.porosityColumn); % m/s

    % Adsorption
    mSma = StericMassActionBinding();
    mSma.kineticBinding = false;
    mSma.lambda  = 962; % mM
    mSma.kA      = [0.0, 1.14];
    mSma.kD      = [0.0, 1];
    mSma.nu      = [1, 4];
    mSma.sigma   = [1, 16];
    mGrm.bindingModel = mSma;

    % Inlet
    nInletSections = 4;

    mGrm.constant   = zeros(nInletSections, mGrm.nComponents);
    mGrm.linear     = zeros(nInletSections, mGrm.nComponents);
    mGrm.quadratic  = zeros(nInletSections, mGrm.nComponents);
    mGrm.cubic      = zeros(nInletSections, mGrm.nComponents);

    % section 1
    mGrm.constant(1,1) = 20.0; % component 1 is salt
    mGrm.constant(1,2) = 0.2;  % component 2 is lysozyme

    % section 2
    mGrm.constant(2,1) = 20.0; % component 1 salt

    % section 3
    mGrm.constant(3,1) = 20.0; % component 1 salt
    switch str
        case '5'
            mGrm.linear(3,1)   = 0.73121; %% 5cv
        case '10'
            mGrm.linear(3,1)   = 0.42022; %% 10cv
        case '30'
            mGrm.linear(3,1)   = 0.14607; %% 30cv
        case '60'
            mGrm.linear(3,1)   = 0.07305; %% 60cv
        case '120'
            mGrm.linear(3,1)   = 0.03524; %% 120cv
    end

    % section 4
    mGrm.constant(4,1) = 526;


    % Configure simulator
    sim = Simulator.create();
    sim.solutionTimes = tOut;
    switch str
        case '5'
            sim.sectionTimes  = [0.0, 12.6, 591, 1284, tOut(end)]; %% 5cv
        case '10'
            sim.sectionTimes  = [0.0, 12.6, 591, 1795, tOut(end)]; %% 10cv
        case '30'
            sim.sectionTimes  = [0.0, 12.6, 591, 4055, tOut(end)]; %% 30cv
        case '60'
            sim.sectionTimes  = [0.0, 12.6, 591, 7518, tOut(end)]; %% 60cv
        case '120'
            sim.sectionTimes  = [0.0, 12.6, 591, 14950, tOut(end)]; %% 120cv
    end
    sim.sectionContinuity = false(nInletSections-1, 1);
    sim.nThreads = 4;

    % Assign model
    sim.model = mGrm;

end
% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%
%  Copyright © 2008-2014: Eric von Lieres¹, Joel Andersson,
%                         Andreas Puettmann¹, Sebastian Schnittert¹,
%                         Samuel Leweke¹
%
%    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
