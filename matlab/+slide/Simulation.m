classdef Simulation < handle
    properties (SetAccess = private)
        Model
        Experiment
        ParameterValues
        Device (1,1) string
        Solution
    end

    methods
        function obj = Simulation(model, options)
            arguments
                model (1,1) slide.SPM = slide.SPM()
                options.experiment (1,1) slide.Experiment = slide.Experiment("Rest for 1 second")
                options.parameter_values (1,1) slide.ParameterValues = slide.ParameterValues("Chen2020")
                options.device (1,1) string = "cpu"
            end
            devices = slide.availableDevices();
            device = lower(options.device);
            if ~any(devices == device)
                error('slide:Device', 'device=%s is unavailable; devices: %s', ...
                    device, strjoin(devices, ', '));
            end
            obj.Model = model;
            obj.Experiment = options.experiment;
            obj.ParameterValues = options.parameter_values;
            obj.Device = device;
            obj.Solution = [];
        end

        function solution = solve(obj, options)
            arguments
                obj
                options.initial_soc double = []
            end
            [names, values, variation_names, variation_values] = ...
                marshal(obj.ParameterValues);
            if ~isempty(options.initial_soc)
                validateattributes(options.initial_soc, {'double'}, ...
                    {'scalar', 'real', 'finite', '>=', 0, '<=', 1});
                names{end + 1, 1} = 'Initial state-of-charge';
                values(end + 1, 1) = options.initial_soc;
            end
            [option_names, option_values] = marshal(obj.Model);
            native = slide_mex('solve', char(obj.ParameterValues.Source), ...
                names, values, option_names, option_values, ...
                obj.Experiment.Steps, obj.Experiment.Period, ...
                char(obj.Device), variation_names, variation_values);
            solution = slide.Solution(native);
            obj.Solution = solution;
        end
    end
end
