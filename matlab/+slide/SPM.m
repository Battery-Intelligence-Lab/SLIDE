classdef SPM
    properties (SetAccess = private)
        Options (1,1) struct
    end

    methods
        function obj = SPM(options)
            arguments
                options (1,1) struct = struct()
            end
            defaults = struct('nch', 8, 'thermal', 0, 'sei_mask', 0, ...
                'sei_porosity', 0, 'crack_mask', 0, 'crack_diffusivity', 0, ...
                'lam_mask', 0, 'plating', 0);
            names = fieldnames(options);
            for index = 1:numel(names)
                if ~isfield(defaults, names{index})
                    error('slide:Options', 'unknown SPM option: %s', names{index});
                end
                defaults.(names{index}) = options.(names{index});
            end
            obj.Options = defaults;
        end

        function [names, values] = marshal(obj)
            names = fieldnames(obj.Options);
            values = zeros(numel(names), 1);
            for index = 1:numel(names)
                value = obj.Options.(names{index});
                validateattributes(value, {'numeric', 'logical'}, ...
                    {'scalar', 'real', 'finite'});
                values(index) = double(value);
            end
        end
    end
end
