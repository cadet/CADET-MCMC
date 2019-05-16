function fitBypassDextran
%==============================================================================
% Fit the triple bypass dextran peak curves
%==============================================================================

    % Number of optimized parameters
    opt.params = 3;

    % Specify the searching domain boundary for the algorithm
    % Length of DPFR, axial dispersion of DPFR, volume of CSTR
    opt.paramBound = [0.01, 1; 1e-11, 1e-8; 1e-11, 1e-7]; % [m, m^2/s, m^3]

    % Call the algorithm in OptAlgo to implement parameter estimation
    OptAlgo.Markov_Chain_Monte_Carlo(@objectiveFunc, opt);

end % fitBypassDextran

function res = objectiveFunc(params)

    % Load the triple experimental datasets
    dataset = cell(3, 1);
    dataset{1} = load('bypassDextran1.dat');
    dataset{2} = load('bypassDextran2.dat');
    dataset{3} = load('bypassDextran3.dat');

    % Generate chromatographic models and invoke the simulator
    solution = createModelBypassDextran(dataset{3}(1:201, 1), params);

    % Calculate the least square residuals
    res = [];
    for i = 1:3
        res = [res; (solution(1:201, 2) - dataset{i}(1:201, 2))];
    end

    res = sum(res.^2);

    % Visualization
    figure(01); clf
    plot(dataset{1}(1:201,1), dataset{1}(1:201,2), 'b:'); hold on
    plot(dataset{2}(1:201,1), dataset{1}(1:201,2), 'b:'); hold on
    plot(dataset{3}(1:201,1), dataset{1}(1:201,2), 'b:'); hold on
    plot(solution(1:201,1), solution(1:201,2), 'r'); hold off
    legend('Meas1', 'Meas2', 'Meas3', 'Sim');
    grid on

end % objectiveFunc

function solution = createModelBypassDextran(time, params)


    nComp = 1; % set global number of components

    % Inlet unit operation
    mIn = PiecewiseCubicPolyInlet();
    mIn.nComponents = nComp;

    % Reserve space: nSections x nComponents (a section can be thought of being a
    % step in the process, see below)
    mIn.constant    = zeros(2, nComp);
    mIn.linear      = zeros(2, nComp);
    mIn.quadratic   = zeros(2, nComp);
    mIn.cubic       = zeros(2, nComp);

    % Section 1: load
    mIn.constant(1,1)  = 2.70e-3; % Dextran solution [mol/m^3]
    % Section 2: transport
    mIn.constant(2,1)  = 0;  % [mol/m^3]

    % Construct DPFR unit operation
    dpfr1 = LumpedRateModelWithoutPores();
    dpfr1.nComponents = nComp;
    dpfr1.columnLength = params(1); % [m]
    dpfr1.porosity     = 1;
    dpfr1.nCellsColumn = 100;
    dpfr1.dispersionColumn      = params(2); % [m^2/s]
    dpfr1.interstitialVelocity  = 0.5e-6/60; % [m/s]
    dpfr1.crossSectionArea      = pi * (5e-4/2)^2; % fixed
    dpfr1.bindingModel = [];
    dpfr1.nBoundStates = [0];
    dpfr1.initialBulk  = [0]; % [mol/m^3]

    % Symmetric DPFR unit
    dpfr2 = LumpedRateModelWithoutPores();
    dpfr2.nComponents = nComp;
    dpfr2.columnLength = dpfr1.columnLength; % [m]
    dpfr2.porosity     = dpfr1.porosity; % [-]
    dpfr2.nCellsColumn = dpfr1.nCellsColumn; %
    dpfr2.dispersionColumn      = dpfr1.dispersionColumn; % [m^2/s]
    dpfr2.interstitialVelocity  = dpfr1.interstitialVelocity; % [m/s]
    dpfr2.crossSectionArea      = dpfr1.crossSectionArea; %
    dpfr2.bindingModel = [];
    dpfr2.nBoundStates = [0];
    dpfr2.initialBulk  = [0]; % [mol/m^3]

    % Construct CSTR unit operation
    cstr1 = StirredTankModel();
    cstr1.nComponents = nComp;
    cstr1.porosity = 1;
    cstr1.nBoundStates = [0];
    cstr1.initialConcentration = [0.0];
    cstr1.initialVolume = params(3);

    % Symmetric CSTR unit
    cstr2 = StirredTankModel();
    cstr2.nComponents = nComp;
    cstr2.porosity = cstr1.porosity;
    cstr2.nBoundStates = cstr1.nBoundStates;
    cstr2.initialConcentration = cstr1.initialConcentration;
    cstr2.initialVolume = cstr1.initialVolume;


    % Outlet unit operation
    mOut = OutletModel();
    mOut.nComponents = nComp;

    % Assemble system of unit operations
    % Construct ModelSystem and assign unit operations (order determines IDs)
    mSys = ModelSystem();
    mSys.models = [mIn, dpfr1, cstr1, cstr2, dpfr2, mOut]; % here index of unit operations are set here: unitOperarion 0 is inlet
    mSys.connectionStartSection = [0];
    mSys.connections = {[0, 1, -1, -1, 0.5e-6/60; ...
                         1, 2, -1, -1, 0.5e-6/60; ...
                         2, 3, -1, -1, 0.5e-6/60; ...
                         3, 4, -1, -1, 0.5e-6/60; ...
                         4, 5, -1, -1, 0.5e-6/60]}; %connection: [from unit, to unit, vol. flow rate]



    % Create simulator and configure it
    sim = Simulator.create();
    sim.solutionTimes = time; % [s], time points at which solution is computed

    % sectionTimes holds the sections and sectionContinuity indicates whether
    % the transition between two adjacent sections is continuous
    sim.sectionTimes = [0.0 12 sim.solutionTimes(end)]; % [s]
    sim.sectionContinuity = false(length(sim.sectionTimes)-2,1);

    sim.nThreads = 4;
    sim.returnLastState = true;
    sim.returnLastStateSens = false;

    % Hand model over to simulator
    sim.model = mSys;

    % Run the model
    result = sim.run();
    solution = [time, result.solution.outlet{length(mSys.models)}];

end
% =============================================================================
%  CADET-MCMC
%
%  Copyright © 2018-2019: Qiaole He¹, Eric von Lieres¹
%
%    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
%
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
