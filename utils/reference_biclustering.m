function [ref_biclusters] = reference_biclustering(nbiclusters)
%REFERENCE_BICLUSTERING return a structure with NBICLUSTERS according to
%the Padilla-Contant-number dataset see details in https://github.com/padilha/biclustlib

if nbiclusters == 1
ref_biclusters.Bic(1).rows = transpose(1:50);       ref_biclusters.Bic(1).cols = (1:50);
elseif nbiclusters == 2
ref_biclusters.Bic(1).rows = transpose(1:50);       ref_biclusters.Bic(1).cols = (1:50);
ref_biclusters.Bic(2).rows = transpose(51:100);     ref_biclusters.Bic(2).cols = (51:100);    
elseif nbiclusters == 3
ref_biclusters.Bic(1).rows = transpose(1:50);       ref_biclusters.Bic(1).cols = (1:50);
ref_biclusters.Bic(2).rows = transpose(51:100);     ref_biclusters.Bic(2).cols = (51:100);
ref_biclusters.Bic(3).rows = transpose(101:150);    ref_biclusters.Bic(3).cols = (101:150);    
elseif nbiclusters == 4
ref_biclusters.Bic(1).rows = transpose(1:50);       ref_biclusters.Bic(1).cols = (1:50);
ref_biclusters.Bic(2).rows = transpose(51:100);     ref_biclusters.Bic(2).cols = (51:100);
ref_biclusters.Bic(3).rows = transpose(101:150);    ref_biclusters.Bic(3).cols = (101:150);
ref_biclusters.Bic(4).rows = transpose(151:200);    ref_biclusters.Bic(4).cols = (151:200);
elseif nbiclusters == 5    
ref_biclusters.Bic(1).rows = transpose(1:50);       ref_biclusters.Bic(1).cols = (1:50);
ref_biclusters.Bic(2).rows = transpose(51:100);     ref_biclusters.Bic(2).cols = (51:100);
ref_biclusters.Bic(3).rows = transpose(101:150);    ref_biclusters.Bic(3).cols = (101:150);
ref_biclusters.Bic(4).rows = transpose(151:200);    ref_biclusters.Bic(4).cols = (151:200);
ref_biclusters.Bic(5).rows = transpose(201:250);    ref_biclusters.Bic(5).cols = (201:250);
end
ref_biclusters.nbicluster = nbiclusters;

end

