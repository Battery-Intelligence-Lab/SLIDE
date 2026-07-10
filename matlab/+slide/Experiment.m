classdef Experiment
    properties (SetAccess = private)
        Steps (:,1) cell
        Period (1,1) double
    end

    methods
        function obj = Experiment(steps, options)
            arguments
                steps
                options.period = "60 seconds"
            end
            if ischar(steps) || (isstring(steps) && isscalar(steps))
                obj.Steps = {char(string(steps))};
            elseif isstring(steps)
                obj.Steps = cellstr(steps(:));
            elseif iscell(steps) && all(cellfun(@(x) ischar(x) || ...
                    (isstring(x) && isscalar(x)), steps(:)))
                obj.Steps = cellfun(@(x) char(string(x)), steps(:), ...
                    'UniformOutput', false);
            else
                error('slide:Experiment', 'steps must be text or a cell/string array');
            end
            if isempty(obj.Steps)
                error('slide:Experiment', 'an experiment needs at least one step');
            end
            obj.Period = slide.Experiment.durationSeconds(options.period);
        end
    end

    methods (Static, Access = private)
        function seconds = durationSeconds(value)
            if isnumeric(value)
                validateattributes(value, {'double'}, {'scalar', 'real', 'finite', 'positive'});
                seconds = double(value);
                return
            end
            token = regexp(lower(strtrim(char(string(value)))), ...
                '^([0-9]+(?:\.[0-9]+)?)\s*([a-z]+)$', 'tokens', 'once');
            if isempty(token)
                error('slide:Experiment', 'invalid sample period: %s', string(value));
            end
            magnitude = str2double(token{1});
            units = token{2};
            if ismember(units, {'s', 'sec', 'second', 'seconds'})
                scale = 1;
            elseif ismember(units, {'min', 'minute', 'minutes'})
                scale = 60;
            elseif ismember(units, {'h', 'hr', 'hour', 'hours'})
                scale = 3600;
            else
                error('slide:Experiment', 'unsupported sample-period unit: %s', units);
            end
            seconds = magnitude * scale;
            if ~(isfinite(seconds) && seconds > 0)
                error('slide:Experiment', 'sample period must be positive and finite');
            end
        end
    end
end
