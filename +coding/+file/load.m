function spikes = load(animal, index, varargin)
% function load(index, varargin)
%
% loads my results either one day or many index

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParameter('appendTo', []);
ip.addParameter('combine', true); % combine invidual day results (if false, looks for already combined)
ip.addParameter('variables', {}); % which variables to load; if empty, all
ip.addParameter('squeeze', []);
ip.addParameter('shuffle', 'cache'); % whether to cache or load the shuffle file
ip.addParameter('filename_full', []);
ip.parse(varargin{:})
Opt = ip.Results;

if isempty(Opt.squeeze) 

    if numel(index) == 1
        Opt.squeeze = true;
    else
        Opt.squeeze = false;
    end

end


if Opt.combine

    % Iteratively load index
    for day = index
        filename = sprintf("checkpoint_index=%s", join(string(day),","));
        if ~isempty(Opt.filename_full)
            filename_full = Opt.filename_full;
        else
            filename_full = coding.file.filename(animal, index);
        end
        spikes(day) = load(filename_full, Opt.variables{:});
    end

    % Combine across index
    if numel(index) > 1
        sp = spikes(index(1));
        fieldList = util.struct.childEnumerate(sp); % enumerates full set of possible string based addresses of child field
        for day = index(2:end)

            % Per day, per type of analysis, invoke the concatonation method
            for name = string(fieldnames(sp))'
                lambda =  eval(['@coding.file.combineMethods.' name]);
                sp.(name) = lambda(sp.(name), spikes(day).(name));
            end

        end
    end

else
    filename_full = coding.file.filename(animal, index, varargin{:});
    error("Not implemented")
end

if Opt.squeeze
    spikes = nd.removeEmpty(spikes);
end
if ~isempty(Opt.appendTo)
    spikes = util.struct.update(Opt.appendTo, spikes);
end


if ~isempty(Opt.shuffle)
    switch char(Opt.shuffle)
        case 'cache'
            spikes.shuffle = coding.file.shufflematfile(animal, index);
        case 'load'
            spikes.shuffle = coding.file.load(...
                coding.file.shufflefilename(animal, index));
        otherwise
            error("Unrecognized option")
    end
end

