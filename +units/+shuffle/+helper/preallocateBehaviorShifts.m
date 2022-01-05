function newtimes = preallocateBehaviorShifts(shifts, beh, groups, varargin)
% Creates our shifted behavior times per shuffle, per neuron

ip = inputParser;
ip.addParameter('maxSize_timeVector', 1000);
ip.parse(varargin{:})
Opt = ip.Results;

[S, N, G] = size(shifts);
newtimes = repmat(shiftdim(beh.time,-2), [S, N]); % make an Shifts x Neuron x Time new record of shifted times
tic
for s = progress(1:S, 'Title', 'Preallocating shuffle times')
    for g = 1:G
        groupselect = groups.time.groups==g;
        scale = mean(diff(beh.time(groupselect)));
        vector_of_timeshifts = round(shifts(s,:,g)/scale)';
        if sum(groupselect) < Opt.maxSize_timeVector
            newtimes(s,:,groupselect) = util.vector.circshift(beh.time(groupselect), vector_of_timeshifts);
        else
            % Break the problem into digenstible chunks
            indices = find(groupselect);
            for ii = 1:Opt.maxSize_timeVector:length(indices)
                subindices = ii:min(ii+Opt.maxSize_timeVector-1, length(indices));
                newtimes(s,:,subindices) = util.vector.circshift(beh.time(subindices), vector_of_timeshifts);
            end
        end
    end
end
toc
