function fitAcetoneWeight
%==============================================================================
% Fit the single peak injection curves
%
% fitting parameters: theta = {D_ax, K_f, D_p, e_c, e_p}
%==============================================================================


    % Collect data of experiment in cell array (one cell per observed compoenent)
    % Note that time points of the measurements are given in sim.solutionTimes
    temp = load('data_ace.dat');
    data = cell(1,1);
    data{1} = temp(:,2);

    sim = createModel(temp(:,1));

    % Set model parameters and enable sensitivities
    params = cell(5, 1);
    params{1} = makeSensitivity([0], {'COL_DISPERSION'},[-1], [-1], [-1], [-1]);
    params{2} = makeSensitivity([0], {'FILM_DIFFUSION'}, [0], [-1], [-1], [-1]);
    params{3} = makeSensitivity([0], {'PAR_DIFFUSION'},  [0], [-1], [-1], [-1]);
    params{4} = makeSensitivity([0], {'COL_POROSITY'},  [-1], [-1], [-1], [-1]);
    params{5} = makeSensitivity([0], {'PAR_POROSITY'},  [-1], [-1], [-1], [-1]);
    sim.setParameters(params, false(5, 1));

    % Create fit object
    pf = ParameterFit();


    % Specify which components are observed in which factor for each observation / wavelength
    idxComp = cell(1,1);
    idxComp{1} = [1];

    % Add weight to dataset. Because different experiments have different concentration unit.
    % The order of magnitude will affect the Metropolis probability 
    dataWeight = cell(1,1);
    dataWeight{1} = ones(length(temp), 1) .* 1;

    % Add the experiment to the fit
    pf.addExperiment(data, sim, [0], idxComp, [], [], [], dataWeight, 'Acetone', {'Acetone'});

    % Enable logarithmic parameter transformation for all parameters
    pf.parameterTransform = LogParameterTransformation(ones(1, 5), true(1, 5));


    % This variable serves as storage for the plot handles returned by the plot function of parameterFit
    plotHd = [];

    opt.params = params;
    % Specify the searching domain boundary for the algorithm
    opt.paramBound = [1e-11, 1e-5; 1e-8, 1e-4; 1e-11, 1e-6; 0.1, 0.9; 0.01, 0.9];

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
    mGrm.nComponents = 1;
    % Discretization
    mGrm.nCellsColumn = 40;
    mGrm.nCellsParticle = 3;
    mGrm.nBoundStates = ones(mGrm.nComponents, 1);

    % Initial conditions
    mGrm.initialBulk = [0.0];
    mGrm.initialSolid = [0.0];

    % Transport
    mGrm.dispersionColumn          = 1.36e-7;
    mGrm.filmDiffusion             = 4.73e-5;
    mGrm.diffusionParticle         = 7.66e-9;
    mGrm.diffusionParticleSurface  = 0.0;

    % Geometry
    mGrm.columnLength      = 0.025;
    mGrm.particleRadius    = 4.5e-5;
    mGrm.porosityColumn    = 0.261;
    mGrm.porosityParticle  = 0.839;

    volumeFlow = 0.5e-6/60; % 0.5 ml/min --> l/s
    columnRadius = 0.35e-2; % m
    mGrm.interstitialVelocity = volumeFlow / (pi * columnRadius^2 * mGrm.porosityColumn); % m/s

    % Adsorption
    mLinear = LinearBinding();
    mLinear.kineticBinding = false;
    mLinear.kA      = 0;
    mLinear.kD      = 1;
    mGrm.bindingModel = mLinear;

    % Inlet
    temp = load('inlet_ace.dat');
    mIn = PiecewiseCubicPolyProfile.fromUniformData(temp(:,1), temp(:,2));

    % Nonnegative limitation
    mIn.constant(mIn.constant < 0) = 0;
    mIn.linear(mIn.linear < 0) = 0;
    mIn.quadratic(mIn.quadratic < 0) = 0;
    mIn.cubic(mIn.cubic < 0) = 0;

    mGrm.inlet = mIn;

    % Simulator
    sim = Simulator.create();
    sim.solutionTimes = tOut;
    sim.sectionTimes = mIn.breaks;
    sim.sectionContinuity = mIn.continuity;

    sim.nThreads = 4;
    sim.returnLastState = true;
    sim.returnLastStateSens = false;

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
