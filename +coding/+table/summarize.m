function out = summarize(T, index)
% function out = summaraze(T, index)
% Summarizes info table entries for an index

index_label_order = ["$day", "$epoch", "$tetrode", "$cell"];
index_label_order_wo_dollar = ["day", "epoch", "tetrode", "cell"];

good_index = index > 0;
index_label_order = index_label_order(good_index);
index = index(good_index);
index = string(index);

% Create the filterstring
filtstring = index_label_order(:) + " == " + index(:);
filtstring = join(filtstring," & ");

% Acquire the entries
filtstring = replace(filtstring, '$','T.');
evalstring = sprintf("T(%s, :)", filtstring);
%disp("Filtering with:")
%disp(evalstring);
%disp('');
out = eval(evalstring);

% If entries not unique, initiate summmary proces
if height(out) > 1
    Tnew = out(1, :);
    for field  = string(fieldnames(out))'
        if field == "Properties" || field == "Row" || field == "Variables"
            continue
        end
        if ismember(field, index_label_order_wo_dollar)
            if ~util.isunique(out.(field))
                Tnew.(field) = out(1,:).(field);
            else
                Tnew.(field) = join(string(out.(field)), ", ");
            end
        elseif contains(field, "mean")
            Tnew.(field) = nanmean(out.(field));
        elseif isnumeric(out.(field))
            Tnew.(field) = nansum(out.(field));
        else
            Tnew.(field) = join(unique(out.(field)), " + ");
        end
    end
    out = Tnew;
end

drop_list = index_label_order_wo_dollar(~good_index);
for field = drop_list(:)'
    out.(field) = [];
end
