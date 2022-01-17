function answer = checkpointActive(checkpoint, iS)

if nargin  == 0
    answer = false;
elseif nargin < 2 % standard checkpoint
    answer = checkpoint;
else % within loop checkpoint
    answer = checkpoint && ~islogical(checkpoint) && mod(iS, checkpoint) == 0;
end

