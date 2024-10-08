% HBIC_QUALITY the heterogeneous intra-bicluster variance uused in HBIC
%   HBIC_QUALITY(BICLUSTERS,XDATA,VARS) computes the quality of the
%   heterogeneous biclusters. In the HBIC manuscript is is defined as
%   the heterogeneous intra-bicluster variance (HIV).
%
%   BICLUSTERS is an structure with B biclusters generated by HBIC with 
%   the following properties:
%   	nbicluster,B     - the number of B candidate bicluster
%       RowxNum          - is NxB logical matrix with the row bicluster's positions
%       NumxCol          - is BxP logical matrix with the column bicluster's positions
%       Bic()            - an structure of size B with row and columns with 
%                           the positions of rows and columns, respectively.
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
%   HBIC_QUALITY returns the same input structure (BICLUSTERS) and adds the
%   following property:
%   	quality         - the heterogeneous intra-bicluster variance for
%                         each bicluster in BICLUSTERS
%
%   Examples:
%   -------
%   see demo_heterogeneous_data.m; % STAGE-II: Computation of HIV 
%   see demo_numerical_data.m;     % STAGE-II: Computation of HIV 
%
% -------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2023, A. Jose-Garcia (adan.jose-garcia@univ-lille.fr)
%   November 2023
% -------------------------------------------------------------------------
function [biclusters] = hbic_quality(biclusters,Xdata,vars)
%QUALITY_BICS Computation of the quality for heterogeneous biclusters
%   BICS is data structure with the following attributes:
%       nbicluster, RowxNum, NumxCol, Bic

var_dtype = vars;
nbics = biclusters.nbicluster;

quality = zeros(nbics,1);
sizeRxC = [];

var_cols = variance_cols(Xdata,var_dtype);

for ibic = 1:nbics
    Bmatrix = Xdata(biclusters.RowxNum(:,ibic),biclusters.NumxCol(ibic,:)); 
    Bvar_dtype = var_dtype(biclusters.NumxCol(ibic,:));
    Bvar_cols = var_cols(biclusters.NumxCol(ibic,:));
    
    quality(ibic) = quality_hiv(Bmatrix,Bvar_dtype,Bvar_cols);
    
    %[nr,nc] = size(ibic);
    sizeRxC = [sizeRxC; size(Bmatrix)];
end

biclusters.quality = quality;
biclusters.sizeRC  = sizeRxC; % optional
biclusters.sizeRxC = sizeRxC(:,1).* sizeRxC(:,2); % optional
end

%##########################################################################
% Bicluster criterion for heterogeneous data
%##########################################################################
function value = quality_hiv(Bmatrix,Bvar_dtype,Bvar_cols)

[nrow,ncol] = size(Bmatrix);
nvars_num = 0; nvars_cat = 0;

sum_variance = [];
sum_homegeneity = [];
for ivar = 1:ncol
    if (strcmpi(Bvar_dtype{ivar},'numeric') || strcmpi(Bvar_dtype{ivar},'integer') || strcmpi(Bvar_dtype{ivar},'double')) 
        ivar_ = var(Bmatrix(:,ivar).Variables);
        ivar_ = ivar_ / Bvar_cols(ivar);
        sum_variance = [sum_variance ivar_];
        nvars_num =  nvars_num + 1;
    elseif (strcmpi(Bvar_dtype{ivar},'binary') || strcmpi(Bvar_dtype{ivar},'categorical'))
        
        [~,~,iunique] = unique(Bmatrix(:,ivar).Variables);
        counts = accumarray(iunique,1); % count for unique values 
        mostfrecuent = max(counts);
        ihomo = 1 - (mostfrecuent/nrow);
        sum_homegeneity = [sum_homegeneity ihomo];
        nvars_cat = nvars_cat + 1;
    else
        error('unrecognized data type');
    end
end

% Should be minimized
if isempty(sum_variance) 
    value =  sum(sum_homegeneity)/nvars_cat;
elseif isempty(sum_homegeneity)
    value =  sum(sum_variance)/nvars_num;
else
    value =  (sum(sum_variance)/nvars_num) + (sum(sum_homegeneity)/nvars_cat) ;
end

end

function cvars = variance_cols(dmatrix,var_dtype)

ncol = size(dmatrix,2);
cvars = NaN(1,ncol);

for ivar = 1:ncol
    if (strcmpi(var_dtype{ivar},'integer') || strcmpi(var_dtype{ivar},'double') || strcmpi(var_dtype{ivar},'numeric')) 
        cvars(ivar) = var(dmatrix(:,ivar).Variables);
    end
end

end