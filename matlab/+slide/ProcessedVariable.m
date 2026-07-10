classdef ProcessedVariable
    properties (SetAccess = private)
        Name (1,1) string
        t (:,1) double
        entries double
        data double
    end

    methods
        function obj = ProcessedVariable(name, time, entries)
            obj.Name = string(name);
            obj.t = time(:);
            obj.entries = entries;
            obj.data = entries;
        end

        function values = interpolate(obj, query)
            if size(obj.entries, 2) == 1
                values = interp1(obj.t, obj.entries, query, 'linear');
            else
                values = zeros(numel(query), size(obj.entries, 2));
                for lane = 1:size(obj.entries, 2)
                    values(:, lane) = interp1(obj.t, obj.entries(:, lane), ...
                        query, 'linear');
                end
            end
        end

        function handle = plot(obj, varargin)
            handle = plot(obj.t, obj.entries, varargin{:});
            xlabel('Time [s]');
            ylabel(obj.Name);
        end

        function value = subsref(obj, indexing)
            if strcmp(indexing(1).type, '()') && numel(indexing(1).subs) == 1
                value = obj.interpolate(indexing(1).subs{1});
                if isscalar(value)
                    value = double(value);
                end
                if numel(indexing) > 1
                    value = subsref(value, indexing(2:end));
                end
            else
                value = builtin('subsref', obj, indexing);
            end
        end
    end
end
