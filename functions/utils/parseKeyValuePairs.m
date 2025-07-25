function parsedInputs = parseKeyValuePairs(defaults, varargin)
    % initialize a struct to save the key-value pairs of input parameters
    parsedInputs = struct();

    % set default value
    for i = 1:2:length(defaults)
        parsedInputs.(defaults{i}) = defaults{i+1};
    end

    % ensure the number of parameters to be even 
    if mod(length(varargin), 2) ~= 0
        error('The number of parameters must be even (key-value pairs).');
    end

    % read key-value pairs
    for i = 1:2:length(varargin)
        key = varargin{i};
        value = varargin{i+1};

        if ~ischar(key) && ~isstring(key)
            error('Input key must be string or char.');
        end

        parsedInputs.(key) = value;
    end
end