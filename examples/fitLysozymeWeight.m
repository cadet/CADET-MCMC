function fitLysozyme
%==============================================================================
% Fit the breakthrough curves of lysozyme
%
% fitting parameters: theta = {D_ax, K_f, D_p, e_c, e_p, k_eq, nu, sigma}
%==============================================================================


    % Collect data of experiment in cell array (one cell per observed compoenent)
    % Note that time points of the measurements are given in sim.solutionTimes
    temp = load('data_lyz.dat'); % data = [time, conc]
    data = cell(1,1);
    data{1} = temp(:,2);

    sim = createModel(temp(:,1));

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

    % Create fit object
    pf = ParameterFit();


    % Specify which components are observed in which factor for each observation / wavelength
    idxComp = cell(1,1);
    idxComp{1} = [0, 1];

    % Add weight to dataset. Because different experiments have different concentration unit.
    % The order of magnitude will affect the Metropolis probability
    dataWeight = cell(1,1);
    dataWeight{1} = ones(length(temp), 1) .* 100;

    % If you have several proteins, you can also take mass of chromatograms into account
%    mass = cellfun(@(x) trapz(temp(:,1), x), data);

    % Add the experiment to the fit
%    pf.addExperiment(data, sim, [0], idxComp, [], [], 1./mass, dataWeight, 'Lysozyme', {'Lysozyme'});
    pf.addExperiment(data, sim, [0], idxComp, [], [], [], dataWeight, 'Lysozyme', {'Lysozyme'});

    % Enable logarithmic parameter transformation for all parameters
    pf.parameterTransform = LogParameterTransformation(ones(1, 8), true(1, 8));


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

function sim = createModel(tOut)
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
    mGrm.linear(3,1)   = 0.14607; %% 30cv

    % section 4
    mGrm.constant(4,1) = 526;


    % Configure simulator
    sim = Simulator.create();
    sim.solutionTimes = tOut;
    sim.sectionTimes  = [0.0, 12.6, 591, 4055, tOut(end)]; %% 30cv
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
