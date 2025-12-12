function resultsToCSV(results, filename)
% resultsToCSV Export detected microtubule anchor points into CSV format
%
% resultsToCSV(results, filename)
%   results   : detection results (struct, same as used in resultsToMask)
%   filename  : output CSV filename
%
% Output CSV format:
%   IDs,X [A],Y [A],Z [A]

    % Open file
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open file %s for writing.', filename);
    end
    
    % Write header
    fprintf(fid, 'IDs,X [A],Y [A],Z [A]\n');

    % Number of objects
    nObj = 0;
    if ~isempty(results) && isfield(results,'data') && ~isempty(results.data)
        nObj = max(nObj, numel(results.data));
    end

    if ~isempty(results) && isfield(results,'points') && ~isempty(results.points)
        nObj = max(nObj, numel(results.points));
    end

    if isfield(results,'center_x') && ~isempty(results.center_x)
        nObj = max(nObj, numel(results.center_x));
    end

    % Loop over all objects
    for k = 1:nObj
        pts = [];

        % Preferred: interpolated curve
        if isfield(results,'data') && ~isempty(results.data) && k <= numel(results.data)
            D = results.data{k};
            if ~isempty(D) && size(D,2) >= 2
                pts = [double(D(:,1)), double(D(:,2))];
            end
        end

        % Fallback: fitted points
        if isempty(pts) && isfield(results,'points') && ~isempty(results.points) && k <= numel(results.points)
            P = results.points{k};
            if ~isempty(P)
                try
                    xk = arrayfun(@(q) double(q.x(1).value), P);
                    yk = arrayfun(@(q) double(q.x(2).value), P);
                catch
                    xk = arrayfun(@(q) double(q.x(1)), P);
                    yk = arrayfun(@(q) double(q.x(2)), P);
                end
                pts = [xk(:), yk(:)];
            end
        end

        % Fallback: single center
        if isempty(pts)
            if isfield(results,'center_x') && numel(results.center_x) >= k ...
               && isfield(results,'center_y') && numel(results.center_y) >= k
                pts = [double(results.center_x(1,k)), double(results.center_y(1,k))];
            elseif isfield(results,'com_x') && size(results.com_x,2) >= k ...
                   && isfield(results,'com_y') && size(results.com_y,2) >= k
                pts = [double(results.com_x(1,k)), double(results.com_y(1,k))];
            else
                continue
            end
        end

        % Write points to CSV (Z=0 for all points)
        for p = 1:size(pts,1)
            fprintf(fid, '%.1f,%.15g,%.15g,0.0\n', k, pts(p,1)-1, pts(p,2)-1);
        end
    end

    fclose(fid);

    %fprintf('Results saved to: %s', filename)
end
