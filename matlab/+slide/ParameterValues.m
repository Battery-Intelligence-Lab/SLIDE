classdef ParameterValues < handle
    properties (SetAccess = private)
        Source (1,1) string
    end

    properties (Access = private)
        Base
        Overrides
        Variations
    end

    methods
        function obj = ParameterValues(source)
            if nargin == 0
                source = "Chen2020";
            end
            obj.Source = string(source);
            [names, values] = slide_mex('parameters', char(obj.Source));
            obj.Base = containers.Map(names, values, 'UniformValues', false);
            obj.Overrides = containers.Map('KeyType', 'char', 'ValueType', 'double');
            obj.Variations = containers.Map('KeyType', 'char', 'ValueType', 'any');
        end

        function value = get(obj, name)
            key = char(string(name));
            if isKey(obj.Variations, key)
                value = slide.VariedParameter(obj.Variations(key));
            elseif isKey(obj.Overrides, key)
                value = obj.Overrides(key);
            elseif isKey(obj.Base, key)
                value = obj.Base(key);
            else
                error('slide:ParameterValues', 'unknown parameter: %s', key);
            end
        end

        function set(obj, name, value)
            key = char(string(name));
            if isa(value, 'slide.VariedParameter')
                obj.Variations(key) = value.Values;
                if isKey(obj.Overrides, key)
                    remove(obj.Overrides, key);
                end
                return
            end
            validateattributes(value, {'double'}, {'scalar', 'real', 'finite'});
            obj.Overrides(key) = double(value);
            if isKey(obj.Variations, key)
                remove(obj.Variations, key);
            end
        end

        function update(obj, values)
            if isa(values, 'containers.Map')
                names = keys(values);
                for index = 1:numel(names)
                    obj.set(names{index}, values(names{index}));
                end
            elseif iscell(values) && size(values, 2) == 2
                for index = 1:size(values, 1)
                    obj.set(values{index, 1}, values{index, 2});
                end
            else
                error('slide:ParameterValues', ...
                    'update expects a containers.Map or an N-by-2 cell array');
            end
        end

        function result = copy(obj)
            result = slide.ParameterValues(obj.Source);
            override_names = keys(obj.Overrides);
            for index = 1:numel(override_names)
                set(result, override_names{index}, obj.Overrides(override_names{index}));
            end
            variation_names = keys(obj.Variations);
            for index = 1:numel(variation_names)
                set(result, variation_names{index}, ...
                    slide.VariedParameter(obj.Variations(variation_names{index})));
            end
        end

        function [names, values, variation_names, variation_values] = marshal(obj)
            names = reshape(keys(obj.Overrides), [], 1);
            values = zeros(numel(names), 1);
            for index = 1:numel(names)
                values(index) = obj.Overrides(names{index});
            end
            variation_names = reshape(keys(obj.Variations), [], 1);
            variation_values = cell(numel(variation_names), 1);
            for index = 1:numel(variation_names)
                variation_values{index} = obj.Variations(variation_names{index});
            end
        end

        function value = subsref(obj, indexing)
            if strcmp(indexing(1).type, '()') && numel(indexing(1).subs) == 1
                value = obj.get(indexing(1).subs{1});
                if numel(indexing) > 1
                    value = subsref(value, indexing(2:end));
                end
            else
                value = builtin('subsref', obj, indexing);
            end
        end

        function obj = subsasgn(obj, indexing, value)
            if strcmp(indexing(1).type, '()') && numel(indexing(1).subs) == 1
                if numel(indexing) ~= 1
                    error('slide:ParameterValues', 'nested parameter assignment is unsupported');
                end
                obj.set(indexing(1).subs{1}, value);
            else
                obj = builtin('subsasgn', obj, indexing, value);
            end
        end
    end
end
