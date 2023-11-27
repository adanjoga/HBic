% HBIC_DISCRETIZATION Group numeric attributes into bins or categories.
%   MDATA = HBIC(XDATA, VARS, NBINS, DATASET_TYPE) is function that groups 
%   the numeric attributes into bins.
%
%   XDATA is an N-by-P data matrix of data-type TABLE with one row per 
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
%   NBINS is a scalar integer, divides the range a numeric attribute into N uniform bins, 
%
%   DATASET_TYPE is an string indicating the overall data-typesin XDATA
%       'ds_numeric'    - All attributes in XDATA are numeric.
%       'ds_mixed'      - Attributes in XDATA are if diffrent data-type.
%
%   DISCRETIZATION return MDATA, a table matrix with numeric attributes
%   discretized according to NBINS.
%
% -------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2023, A. Jose-Garcia (adan.jose-garcia@univ-lille.fr)
%   November 2023
% -------------------------------------------------------------------------

function [data] = hbic_discretization(data,vars,nbins,dataset_type)

var_dtype = vars;

if nargin < 3
    nbins = 5;
    dataset_type = 'ds_numeric';
end

if strcmpi(dataset_type,'ds_numeric')
    data_values = table2array(data);
    data_values = discretize(data_values,nbins);
    data(:,:) = array2table(data_values);
    figure(2); imagesc(table2array(data));
    
elseif strcmpi(dataset_type,'ds_mixed')
    
    for i_att = 1:size(var_dtype,1)
        att_values = data(:,i_att).Variables; % from table to array

        if (strcmpi(var_dtype{i_att},'Numeric') || strcmpi(var_dtype{i_att},'Integer') || strcmpi(var_dtype{i_att},'Double')) 
            att_values = normalize(att_values,'range');
            att_values = discretize(att_values,nbins);
            data(:,i_att) = table(att_values);
        end
    end
    %figure(2); imagesc(table2array(data));
end

%figure(2); imagesc(table2array(data));
end
