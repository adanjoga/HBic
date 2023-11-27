%% DEMO: Evaluation of HBIC on numerical datasets
% -------------------------------------------------------------------------
% HBIC is biclustering algorithm for heterogeneous and missing data.
% HBIC handle mixed-type data such as numerical, binary and categorical.
% HBIC HBic automatically identifies the number of biclusters or 
% takes this parameter as input if this knowledge if available. 
%
% -------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2023, A. Jose-Garcia (adan.jose-garcia@univ-lille.fr)
%   November 2023
% -------------------------------------------------------------------------
%% Reading dataset and reference biclustering solution
% -------------------------------------------------------------------------
clearvars; clc; close all;

addpath([pwd '/datasets/numeric-padilla-data']);
addpath([pwd '/hbic']);
addpath([pwd '/metrics']);
addpath([pwd '/utils']);

nbiclusters = 5; % posible options: {1,2,3,4,5}

% reading the dataset and creating the numeric data-type  vector
Xdata = readtable(['number_' num2str(nbiclusters) '_1.txt']);
Xdata = Xdata(:,2:end);
var_dtype = repmat("numeric",[1 size(Xdata,2)]);
[nr,nc] = size(Xdata);

% Creation of the reference biclustering solution
ref_biclusters = reference_biclustering(nbiclusters);
T = ref_biclusters.Bic;

%% STAGE-I: Generation of Candidate Biclusters
% -------------------------------------------------------------------------
nbins = 5;      % number of bins for the discretization function
Mdata           = hbic_discretization(Xdata,var_dtype,nbins,'data_numeric');
Bics_candidate  = hbic_algorithm(Mdata,var_dtype);

% evaluating the performance of the candidate biclusters
B = Bics_candidate.Bic;
metric_recovery     = external_biclustering_indices('prec', B, T, nr, nc);
metric_relevance    = external_biclustering_indices('prel', B, T, nr, nc);

%% STAGE-II: Bicluster Model Selection (OPTIONAL)
% -------------------------------------------------------------------------
% computation of HIV  (heterogeneous intra-bicluster variance)
Bics_candidate  = hbic_quality(Bics_candidate,Xdata,var_dtype); 

% selection types: {'all','top_b','tree_b','merge_b','tree_auto','merge_auto'}
selection_type  = 'top_b'; 
Biclusters      = hbic_selection(selection_type, Bics_candidate, ref_biclusters.nbicluster);

% evaluating the performance of the final-selected biclusters
B = Biclusters.Bic;
metric_recovery2     = external_biclustering_indices('prec', B, T, nr, nc);
metric_relevance2    = external_biclustering_indices('prel', B, T, nr, nc);

