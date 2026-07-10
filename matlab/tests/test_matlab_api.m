function tests = test_matlab_api
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testFile = mfilename('fullpath');
testCase.TestData.repo = fileparts(fileparts(fileparts(testFile)));
addpath(fullfile(testCase.TestData.repo, 'matlab'));
end

function testParametersBpxAndUpdate(testCase)
parameters = slide.ParameterValues('Chen2020');
verifyEqual(testCase, get(parameters, 'Nominal cell capacity [A.h]'), 5.0);
parameters('Initial state-of-charge') = 0.7;
update(parameters, {'Contact resistance [Ohm]', 2e-3});
verifyEqual(testCase, parameters('Initial state-of-charge'), 0.7);
verifyEqual(testCase, parameters('Contact resistance [Ohm]'), 2e-3);

copied = copy(parameters);
copied('Initial state-of-charge') = 0.6;
verifyEqual(testCase, parameters('Initial state-of-charge'), 0.7);

fixture = fullfile(testCase.TestData.repo, 'matlab', 'tests', 'fixture.bpx.json');
bpx = slide.ParameterValues(fixture);
verifyEqual(testCase, bpx('Nominal cell capacity [A.h]'), 5.0);
simulation = slide.Simulation(slide.SPM(struct('nch', 12)), ...
    experiment=slide.Experiment('Rest for 1 second', period='1 second'), ...
    parameter_values=bpx);
solution = solve(simulation);
verifyTrue(testCase, all(isfinite(variable(solution, 'Voltage [V]').entries)));
end

function testPythonWheelSampleParity(testCase)
reference = readmatrix(fullfile(testCase.TestData.repo, 'tests', 'reference', ...
    'slide_python_chen2020_1c_600.csv'), 'NumHeaderLines', 1);
simulation = slide.Simulation(slide.SPM(struct('nch', 12)), ...
    experiment=slide.Experiment('Discharge at 1 C for 600 seconds', ...
                                period='10 seconds'));
solution = solve(simulation, initial_soc=0.8);
voltage = variable(solution, 'Voltage [V]');
current = variable(solution, 'Current [A]');
verifyEqual(testCase, solution.t, reference(:, 1));
verifyEqual(testCase, current.entries, reference(:, 2));
maximumError = max(abs(voltage.entries - reference(:, 3)));
fprintf('P8-G1 MATLAB/Python maximum voltage error: %.17g V\n', maximumError);
verifyLessThanOrEqual(testCase, maximumError, 2e-12);
end

function testTutorialFiveParity(testCase)
steps = {
    'Discharge at C/10 for 10 hours or until 3.3 V'
    'Rest for 1 hour'
    'Charge at 1 A until 4.1 V'
    'Hold at 4.1 V until 50 mA'
    'Rest for 1 hour'
};
reference = readmatrix(fullfile(testCase.TestData.repo, 'tests', 'reference', ...
    'pybamm_chen2020_tutorial5.csv'), 'NumHeaderLines', 1);
simulation = slide.Simulation(slide.lithium_ion.SPM(struct('nch', 12)), ...
    experiment=slide.Experiment(steps, period='60 seconds'));
solution = solve(simulation, initial_soc=1.0);
voltage = variable(solution, 'Terminal voltage [V]').entries;
errors = [];
for segment = 1:numel(steps)
    referenceMask = reference(:, 1) == segment - 1;
    slideMask = solution.sample_segment == segment;
    if segment == 1
        start = 0;
    else
        previous = solution.t(solution.sample_segment == segment - 1);
        start = previous(end);
    end
    localTime = solution.t(slideMask) - start;
    comparable = referenceMask ...
        & reference(:, 3) >= localTime(1) - 1e-9 ...
        & reference(:, 3) <= localTime(end) + 1e-9;
    interpolated = interp1(localTime, voltage(slideMask), ...
        reference(comparable, 3), 'linear', 'extrap');
    errors = [errors; interpolated - reference(comparable, 5)]; %#ok<AGROW>
end
verifyFalse(testCase, any(isnan(errors)));
maximumError = max(abs(errors));
rmsError = sqrt(mean(errors .^ 2));
fprintf('P8-G1 Tutorial-5 maximum/RMS error: %.6g/%.6g V\n', ...
    maximumError, rmsError);
verifyLessThanOrEqual(testCase, maximumError, 15e-3);
verifyLessThanOrEqual(testCase, rmsError, 8e-3);
end

function testProcessedVariablePlotSaveAndVaried(testCase)
parameters = slide.ParameterValues('Chen2020');
set(parameters, 'Initial state-of-charge', slide.varied([0.4, 0.6, 0.8]));
simulation = slide.Simulation(slide.SPM(struct('nch', 12)), ...
    experiment=slide.Experiment('Discharge at 1 C for 60 seconds', ...
                                period='10 seconds'), ...
    parameter_values=parameters);
solution = solve(simulation);
voltage = variable(solution, 'Voltage [V]');
verifySize(testCase, voltage.entries, [7, 3]);
verifyEqual(testCase, voltage(15), interp1(solution.t, voltage.entries, 15));

single = solve(slide.Simulation(slide.SPM(), ...
    experiment=slide.Experiment('Rest for 2 seconds', period='1 second')));
csvPath = [tempname, '.csv'];
matPath = [tempname, '.mat'];
cleanup = onCleanup(@() deleteIfPresent({csvPath, matPath})); %#ok<NASGU>
saveData(single, csvPath, 'csv');
saveData(single, matPath, 'mat');
verifyTrue(testCase, isfile(csvPath));
verifyTrue(testCase, isfile(matPath));
figureHandle = figure('Visible', 'off');
figureCleanup = onCleanup(@() close(figureHandle)); %#ok<NASGU>
plot(variable(single, 'Voltage [V]'));
end

function testStableErrorsRecoverAndRepeat(testCase)
badOptions = slide.Simulation(slide.SPM(struct('nch', 7)), ...
    experiment=slide.Experiment('Rest for 1 second', period='1 second'));
verifyError(testCase, @() solve(badOptions), 'slide:Options');

missing = fullfile(testCase.TestData.repo, 'matlab', 'tests', 'missing.bpx.json');
verifyError(testCase, @() slide.ParameterValues(missing), 'slide:BPX');

badStep = slide.Simulation(slide.SPM(), ...
    experiment=slide.Experiment('This is not an experiment step'));
verifyError(testCase, @() solve(badStep), 'slide:Experiment');

for repetition = 1:50
    simulation = slide.Simulation(slide.SPM(), ...
        experiment=slide.Experiment('Rest for 1 second', period='1 second'));
    solution = solve(simulation, initial_soc=0.55);
    verifyEqual(testCase, solution.t, [0; 1]);
    verifyTrue(testCase, all(isfinite(variable(solution, 'Voltage [V]').entries)));
    clear simulation solution
end
end

function testDevicePreflight(testCase)
devices = slide.availableDevices();
verifyEqual(testCase, devices(1), "cpu");
verifyError(testCase, @() slide.Simulation(slide.SPM(), device="gpu"), ...
    'slide:Device');
end

function deleteIfPresent(paths)
for index = 1:numel(paths)
    if isfile(paths{index})
        delete(paths{index});
    end
end
end
