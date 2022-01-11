function newtimes = preallocateBehaviorShifts(shifts, beh, groups, varargin)
% Creates our shifted behavior times per shuffle, per neuron
%
% TODO add the yarsev shift type

ip = inputParser;
ip.addParameter('maxSize_timeVector', 1000);
ip.addParameter('shifttype', 'circshift');
ip.parse(varargin{:})
Opt = ip.Results;

[S, N, G] = size(shifts);
newtimes = repmat(shiftdim(beh.time,-2), [S, N]); % make an Shifts x Neuron x Time new record of shifted times
tic

iscircshift = strcmp(Opt.shifttype, 'circshift');
disp("Using shift type = " + Opt.shifttype);

for s = progress(1:S, 'Title', 'Preallocating shuffle times')
    for g = 1:G
        groupselect = groups.time.groups==g;
        scale = mean(diff(beh.time(groupselect)));
        vector_of_timeshifts = round(shifts(s,:,g)/scale)';

        if iscircshift
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
        else % LINEAR INDEX MOVEMENT
            keyboard; % untested!
            newtimes(s,:,groupselect) = beh.time(groupselect + vector_of_timeshifts);

            
            %NOTE : I do not have implented the shuffle in yatrsev dotson
            %where indices shift up to the edge of a boundary and disappear
        end

    end
end
toc
