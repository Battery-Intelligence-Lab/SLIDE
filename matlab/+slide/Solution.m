classdef Solution
    properties (SetAccess = private)
        t (:,1) double
        termination (1,1) string
        termination_detail (1,1) string
        segment (1,1) double
        sample_segment (:,1) double
        n_lanes (1,1) double
    end

    properties (Access = private)
        Fields
    end

    methods
        function obj = Solution(native)
            obj.t = native.time(:);
            obj.termination = string(native.termination);
            obj.termination_detail = string(native.termination_detail);
            obj.segment = native.segment;
            obj.sample_segment = native.sample_segment(:);
            obj.n_lanes = native.n_lanes;
            obj.Fields = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.Fields('Time [s]') = obj.t;
            obj.Fields('Terminal voltage [V]') = native.voltage;
            obj.Fields('Voltage [V]') = native.voltage;
            obj.Fields('Current [A]') = native.current;
        end

        function result = variable(obj, name)
            key = char(string(name));
            if ~isKey(obj.Fields, key)
                error('slide:Solution', 'variable was not recorded: %s', key);
            end
            result = slide.ProcessedVariable(key, obj.t, obj.Fields(key));
        end

        function handles = plot(obj, variables)
            if nargin < 2
                variables = "Terminal voltage [V]";
            end
            variables = string(variables);
            handles = gobjects(numel(variables), 1);
            tiledlayout(numel(variables), 1);
            for index = 1:numel(variables)
                nexttile;
                field = variable(obj, variables(index));
                handles(index) = plot(field);
            end
        end

        function saveData(obj, filename, format)
            if nargin < 3
                [~, ~, extension] = fileparts(filename);
                format = erase(lower(string(extension)), '.');
            end
            time = obj.t; %#ok<NASGU>
            voltage = obj.Fields('Terminal voltage [V]'); %#ok<NASGU>
            current = obj.Fields('Current [A]'); %#ok<NASGU>
            if any(strcmpi(format, {'mat', 'matlab'}))
                save(filename, 'time', 'voltage', 'current');
            elseif strcmpi(format, 'csv')
                if size(voltage, 2) ~= 1
                    error('slide:Solution', 'CSV save currently requires one lane');
                end
                table_data = table(time, voltage, current, ...
                    'VariableNames', {'Time_s', 'Terminal_voltage_V', 'Current_A'});
                writetable(table_data, filename);
            else
                error('slide:Solution', 'format must be csv or mat');
            end
        end

    end
end
