% HBIC A biclustering algorithm for heterogeneous and missing data.
%   HBIC(DATA, VARS) finds coherent submatrices in mixed-type datasets.
%
%   DATA is an N-by-P data matrix of data-type TABLE with one row per 
%   observation and one column per variable. Each column on DATA can be 
%   of different data-type such as numeric, binary or categorical. The
%   datatype of each variable is indicated in VARS.
%
%   VARS is 1-by-P row vector of type string with different datatypes 
%   such as 'numeric', 'binary', and 'categorical'. Each values corresponds
%   to one column, variable in DATA. The available datatypes are:
%       'numeric'       - Euclidean distance (the default).
%       'binary'        - Normalized Euclidean distance.
%       'categorical'   - Cosine similarity.
%   
%   HBIC returns an structure (object) with the following properties:
%   	nbicluster,B     - the number of B candidate bicluster
%       RowxNum          - is NxB logical matrix with the row bicluster's positions
%       NumxCol          - is BxP logical matrix with the column bicluster's positions
%       Bic()            - an structure of size B with row and columns with 
%                           the positions of rows and columns, respectively.
%
%
%   Examples:
%   -------
%   see demo_heterogeneous_data.m; % STAGE-I: Generation of Candidate Biclusters
%   see demo_numerical_data.m;     % STAGE-I: Generation of Candidate Biclusters
%
% -------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2023, A. Jose-Garcia (adan.jose-garcia@univ-lille.fr)
%   November 2023
% -------------------------------------------------------------------------
function [biclusters] = hbic_algorithm(data,vars)
[nrows, ncols] = size(data);

n_mincols = max(round((ncols*3)/100),2);
n_minrows = max(round((nrows*3)/100),2);
var_dtype = vars;

nbicluster = 1;
biclusters_temp = [];
for icol = 1:ncols
    
    disp(['Iteration: ' num2str(icol)]);
    
    c_values = data(:,icol).Variables; % from table to array (cell or double)
    c_datatype = var_dtype(icol);
    cunique_vals = unique(c_values,'stable');
    for ival = 1:numel(cunique_vals)
        
        col_value = cunique_vals(ival);
        if isa(col_value,'cell')
            idx_rows = strcmpi(col_value,c_values);
        else
            idx_rows = (c_values==col_value);
        end
        
        if sum(idx_rows) < n_minrows
            continue;
        end

        idx_cols = false(ncols,1);
        idx_cols(icol) = true; % WARNING

        bic_size = sum(idx_rows)*sum(idx_cols);
        mybic = struct('lrows',idx_rows,'lcols',idx_cols,'size',bic_size);
        STOP = false;

        while not(STOP)

            % Step 1: Finding the best column for the bicluster
            H = zeros(1,ncols);
            H_values = cell(1,ncols);
            iter_cols = transpose(find(~mybic.lcols)); % WARNING: what if empty?
            temp_data = data(mybic.lrows,:);
            for jcol = iter_cols
                
                temp_cdtype = var_dtype(jcol);

                temp_cvalues = temp_data(:,jcol).Variables;
                [C,~,iunique] = unique(temp_cvalues);
                counts = accumarray(iunique,1);     % count the unique values 
                [jc_max,pmax] = max(counts);        % position of the max count

                if isa(C(pmax),'cell')
                	H_values(jcol)=C(pmax);
                else
                    H_values{jcol}=C(pmax);
                end
                
                H(jcol)=jc_max;
            end
            tbic = mybic;
            % Adding the best column to the bicluster
            [~,cbest_idx] = max(H);
            tbic.lcols(cbest_idx) = true; 

            % Step: Getting a homogeneous bicluster
            cbest_value = H_values{cbest_idx}; % num (double) or cat (cell)

            cbest_values = data(:,cbest_idx).Variables;
            if isa(cbest_values,'cell')
                cbest_rows = strcmpi(cbest_value,cbest_values); 
            else
                cbest_rows = (cbest_values==cbest_value);
            end
            
            % Selecting the row indices matching the new bicluster
            tbic.lrows = tbic.lrows & cbest_rows;

            % Compute bicluster size
            tbic.size = sum(tbic.lrows) * sum(tbic.lcols);

            if (tbic.size > mybic.size)
                mybic = tbic;
            else
                STOP = true;
                %disp('WARNING STOP');
            end
        end
        % Validate the generated biclusters
        if sum(mybic.lrows)>n_minrows && sum(mybic.lcols)>n_mincols
            disp(['Bicluster (' num2str(nbicluster) '): ' num2str(sum(mybic.lrows)) 'x' num2str(sum(mybic.lcols))])
            nbicluster = nbicluster + 1;
            % ADD the new bicluster
            biclusters_temp = [biclusters_temp mybic];
        else
            % disp(['WARNING (' num2str(icol) ') bicluster of size: ' num2str(sum(mybic.lrows)) 'x' num2str(sum(mybic.lcols))])
        end
    end
end

if (nbicluster-1) > 0
    % Remove repeated biclusters
    biclusters_temp = remove_bics(biclusters_temp);
end
% Save unique biclusters in a mat struct
if numel(biclusters_temp) > 0
    biclusters.nbicluster = numel(biclusters_temp);
    biclusters.RowxNum = [biclusters_temp.lrows];
    biclusters.NumxCol = transpose([biclusters_temp.lcols]);
    %res_biclusters.Size = transpose([biclusters.size]);
    for ibic=1:biclusters.nbicluster
        rows = find(biclusters.RowxNum(:,ibic) > 0);
        cols =find(biclusters.NumxCol(ibic,:) > 0);
        biclusters.Bic(ibic) = struct('rows',rows,'cols',cols);
    end
else
    disp('WARNING: No biclusters found')
    biclusters = nan;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ubiclusters] = remove_bics(biclusters)
    nbics = numel(biclusters);
    uniquebics = false(nbics,1);
    uniquebics(1) = true;
    
    RowxNum = [biclusters.lrows];
    NumxCol = transpose([biclusters.lcols]);
    for ibic = 2:nbics
        if isuniqueb(ibic,RowxNum,NumxCol)
            uniquebics(ibic) = true;
        end
    end
    ubiclusters = biclusters(uniquebics);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uniqueb = isuniqueb(ibic,RxN,NxC)
    uniqueb = true;
    for jbic = 1:(ibic-1)
        if all(RxN(:,jbic)==RxN(:,ibic))
            if all(NxC(jbic,:)==NxC(ibic,:))
                uniqueb = false;
                break;
            end
        end
            
    end
end