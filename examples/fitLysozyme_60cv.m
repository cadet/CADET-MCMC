function fitLysozyme
%==============================================================================
% Fit the 60 CV lysozyme peak curve
%==============================================================================

    % Number of optimized parameters
    opt.params = 11;

    % Specify the searching domain boundary for the algorithm
    % Length of DPFR, axial dispersion of DPFR, volume of CSTR
    % axial dispersion of column, column porosity
    % film mass transfer, pore diffusion, particle porosity
    % equilibrium constant, nu, sigma
    opt.paramBound = [0.01 1; 1e-11 1e-8; 1e-11 1e-7;...
        1e-11 1e-9; 0.01 0.90;...
        1e-9 1e-4; 1e-11 1e-7; 0.40 0.99; ...
        1e-3 1.0; 1.0 30; 5.0 80];

    % Call the algorithm in OptAlgo to implement parameter estimation
    OptAlgo.Markov_Chain_Monte_Carlo(@objectiveFunc, opt);

end % fitLysozyme

function res = objectiveFunc(params)

    str = '60cv';
    enablePlot = false;

    % Load dataset
    dataset = load(['columnLysozyme' str '.dat']);

    % Generate chromatographic models and invoke the simulator
    try
        solution = createModelcolumnLysozyme(dataset(1:38000, 1), params, str);
    catch e
        res = 106800;
        save('crashSave.dat', 'params', '-ascii', '-append');
        return
    end

    % Calculate the least square residuals
    res = [(solution(1:38000, 2) - dataset(1:38000, 2))];
    res = sum(res.^2);

    % Visualization if enabled
    if enablePlot
        figure(01); clf
        plot(dataset(14570:38000,1),  dataset(14570:38000,2), 'b:'); hold on
        plot(solution(14570:38000,1), solution(14570:38000,2), 'r'); hold off
        legend('Meas', 'Sim');
        grid on
    end

end % objectiveFunc

function solution = createModelcolumnLysozyme(tOut, params, str)


    nComp = 2; % set global number of components

    % Inlet unit operation
    mIn = PiecewiseCubicPolyInlet();
    mIn.nComponents = nComp;

    % Reserve space: nSections x nComponents (a section can be thought of being a
    % step in the process, see below)
    mIn.constant    = zeros(4, nComp);
    mIn.linear      = zeros(4, nComp);
    mIn.quadratic   = zeros(4, nComp);
    mIn.cubic       = zeros(4, nComp);

    % Section 1: load
    mIn.constant(1, 1)  = 20.0; % salt
    mIn.constant(1, 2)  = 0.2;  % lysozyme
    % Section 2: wash
    mIn.constant(2, 1)  = 20.0;
    % Section 3: elute
    mIn.constant(3, 1)  = 20.0; % salt
    switch str
        case '0cv'
            mIn.constant(3, 1) = 526; % 0cv
        case '5cv'
            mIn.linear(3, 1)   = 0.73121; % 5cv
        case '10cv'
            mIn.linear(3, 1)   = 0.42022; % 10cv
        case '30cv'
            mIn.linear(3, 1)   = 0.14607; % 30cv
        case '60cv'
            mIn.linear(3, 1)   = 0.07305; % 60cv
        case '120cv'
            mIn.linear(3, 1)   = 0.03524; % 120cv
    end
    % Section 4: strip
    mIn.constant(4, 1) = 526; % salt

    % Construct DPFR unit operation
    dpfr1 = LumpedRateModelWithoutPores();
    dpfr1.nComponents = nComp;
    dpfr1.columnLength = params(1);
    dpfr1.porosity     = 1;
    dpfr1.nCellsColumn = 100;
    dpfr1.dispersionColumn      = params(2);
    dpfr1.interstitialVelocity  = 1;
    dpfr1.crossSectionArea      = pi * (5e-4/2)^2; % cross sectional area
    dpfr1.bindingModel = [];
    dpfr1.nBoundStates = [0, 0];
    dpfr1.initialBulk  = [20, 0];

    % Symmetric DPFR unit
    dpfr2 = LumpedRateModelWithoutPores();
    dpfr2.nComponents = nComp;
    dpfr2.columnLength = dpfr1.columnLength;
    dpfr2.porosity     = dpfr1.porosity;
    dpfr2.nCellsColumn = dpfr1.nCellsColumn;
    dpfr2.dispersionColumn      = dpfr1.dispersionColumn;
    dpfr2.interstitialVelocity  = dpfr1.interstitialVelocity;
    dpfr2.crossSectionArea      = dpfr1.crossSectionArea;
    dpfr2.bindingModel = dpfr1.bindingModel;
    dpfr2.nBoundStates = dpfr1.nBoundStates;
    dpfr2.initialBulk  = dpfr1.initialBulk;

    % Construct CSTR unit operation
    cstr1 = StirredTankModel();
    cstr1.nComponents = nComp;
    cstr1.porosity = 1;
    cstr1.nBoundStates = [0, 0];
    cstr1.initialConcentration = [20.0, 0];
    cstr1.initialVolume = params(3);

    % Symmetric CSTR unit
    cstr2 = StirredTankModel();
    cstr2.nComponents = nComp;
    cstr2.porosity = cstr1.porosity;
    cstr2.nBoundStates = cstr1.nBoundStates;
    cstr2.initialConcentration =cstr1.initialConcentration;
    cstr2.initialVolume = cstr1.initialVolume;

    % Construct GRM column unit
    column = GeneralRateModel();
    column.nComponents = nComp;
    column.nCellsColumn = 80;
    column.nCellsParticle = 8;
    column.dispersionColumn          = params(4);
    column.filmDiffusion             = [params(6), params(6)];
    column.diffusionParticle         = [params(7), params(7)];
    column.diffusionParticleSurface  = [0, 0];
    column.interstitialVelocity      = 1;
    column.crossSectionArea          = pi * 0.35e-2^2;
    column.columnLength        = 0.025;
    column.particleRadius      = 4.5e-5;
    column.porosityColumn      = params(5);
    column.porosityParticle    = params(8);
    column.nBoundStates = ones(column.nComponents, 1);
    column.initialBulk = [20.0, 0];
    column.initialSolid = [ionCapacity(column.porosityColumn, column.porosityParticle), 0];

    % Adsorption kinetics
    mSma = StericMassActionBinding();
    mSma.kineticBinding = true;
    mSma.lambda     = ionCapacity(column.porosityColumn, column.porosityParticle);
    mSma.kA         = [0.0 params(9)];
    mSma.kD         = [0.0 1];
    mSma.nu         = [0.0 params(10)];
    mSma.sigma      = [0.0 params(11)];
    column.bindingModel = mSma;

    % Outlet unit operation
    mOut = OutletModel();
    mOut.nComponents = nComp;

    % Assemble system of unit operations
    % Construct ModelSystem and assign unit operations (order determines IDs)
    mSys = ModelSystem();
    mSys.models = [mIn, dpfr1, cstr1, column, cstr2, dpfr2, mOut];
    mSys.connectionStartSection = [0];
    mSys.connections = {[0, 1, -1, -1, 0.5e-6/60; ...
                         1, 2, -1, -1, 0.5e-6/60; ...
                         2, 3, -1, -1, 0.5e-6/60; ...
                         3, 4, -1, -1, 0.5e-6/60; ...
                         4, 5, -1, -1, 0.5e-6/60; ...
                         5, 6, -1, -1, 0.5e-6/60]};


    % Create simulator and configure it
    sim = Simulator.create();
    sim.solutionTimes = tOut; % [s], time points at which solution is computed

    % sectionTimes holds the sections and sectionContinuity indicates whether
    % the transition between two adjacent sections is continuous
    switch str
        case '0cv'
            sim.sectionTimes = [0.0, 12.6, 591, 1168, tOut(end)];
        case '5cv'
            sim.sectionTimes = [0.0, 12.6, 591, 1284, tOut(end)];
        case '10cv'
            sim.sectionTimes = [0.0, 12.6, 591, 1795, tOut(end)];
        case '30cv'
            sim.sectionTimes = [0.0, 12.6, 591, 4055, tOut(end)];
        case '60cv'
            sim.sectionTimes = [0.0, 12.6, 591, 7518, tOut(end)];
        case '120cv'
            sim.sectionTimes = [0.0, 12.6, 591, 14950, tOut(end)];
    end
    sim.sectionContinuity = false(length(sim.sectionTimes)-2,1);

    sim.nThreads = 4;
    sim.returnLastState = false;
    sim.returnLastStateSens = false;

    % Hand model over to simulator
    sim.model = mSys;

    % Run the models
    result = sim.run();
    solution = [tOut, result.solution.outlet{length(mSys.models),1}(:,2)];

end % createModelColumnLysozyme

function capacity = ionCapacity(ec, ep)

    c = 10; % [mol/m^3]
    V = 19.252e-6; % [m^3]
    Vc = 9.6211e-7; % [m^3]
    et = ec + (1 - ec) * ep;

    capacity = c * V / (Vc * (1 - et));

end
% =============================================================================
%  CADET - The Chromatography Analysis and Design Toolkit
%
%  Copyright © 2008-2019: Eric von Lieres¹, Joel Andersson,
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
