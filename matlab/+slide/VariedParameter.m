classdef VariedParameter
    properties (SetAccess = private)
        Values (1,:) double
    end

    methods
        function obj = VariedParameter(values)
            arguments
                values (1,:) double
            end
            if isempty(values) || any(~isfinite(values))
                error('slide:Variation', 'varied values must be non-empty and finite');
            end
            obj.Values = values;
        end
    end
end
