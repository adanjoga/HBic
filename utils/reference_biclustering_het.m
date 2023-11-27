function [ref_biclusters] = reference_biclustering_het(nbiclusters)
%REFERENCE_BICLUSTERING_HET return a structure with NBICLUSTERS according to
%the heteregneous syntetic datasets reported in our manuscript

if nbiclusters == 1
ref_biclusters.Bic(1).rows = transpose(1:50); ref_biclusters.Bic(1).cols = (1:30);
elseif nbiclusters == 3
ref_biclusters.Bic(1).rows = transpose(1:50); ref_biclusters.Bic(1).cols = (1:30);
ref_biclusters.Bic(2).rows = transpose(51:100); ref_biclusters.Bic(2).cols = (31:60);
ref_biclusters.Bic(3).rows = transpose(101:150); ref_biclusters.Bic(3).cols = (61:90);
elseif nbiclusters == 5    
ref_biclusters.Bic(1).rows = transpose(1:50); ref_biclusters.Bic(1).cols = (1:30);
ref_biclusters.Bic(2).rows = transpose(51:100); ref_biclusters.Bic(2).cols = (31:60);
ref_biclusters.Bic(3).rows = transpose(101:150); ref_biclusters.Bic(3).cols = (61:90);
ref_biclusters.Bic(4).rows = transpose(151:200); ref_biclusters.Bic(4).cols = (91:120);
ref_biclusters.Bic(5).rows = transpose(201:250); ref_biclusters.Bic(5).cols = (121:150);
elseif nbiclusters == 8
ref_biclusters.Bic(1).rows = transpose(1:50); ref_biclusters.Bic(1).cols = (1:30);
ref_biclusters.Bic(2).rows = transpose(51:100); ref_biclusters.Bic(2).cols = (31:60);
ref_biclusters.Bic(3).rows = transpose(101:150); ref_biclusters.Bic(3).cols = (61:90);
ref_biclusters.Bic(4).rows = transpose(151:200); ref_biclusters.Bic(4).cols = (91:120);
ref_biclusters.Bic(5).rows = transpose(201:250); ref_biclusters.Bic(5).cols = (121:150);
ref_biclusters.Bic(6).rows = transpose(251:300); ref_biclusters.Bic(6).cols = (151:180);
ref_biclusters.Bic(7).rows = transpose(301:350); ref_biclusters.Bic(7).cols = (181:210);
ref_biclusters.Bic(8).rows = transpose(351:400); ref_biclusters.Bic(8).cols = (211:240);
elseif nbiclusters == 10
ref_biclusters.Bic(1).rows = transpose(1:50);    ref_biclusters.Bic(1).cols = (1:30);
ref_biclusters.Bic(2).rows = transpose(51:100);  ref_biclusters.Bic(2).cols = (31:60);
ref_biclusters.Bic(3).rows = transpose(101:150); ref_biclusters.Bic(3).cols = (61:90);
ref_biclusters.Bic(4).rows = transpose(151:200); ref_biclusters.Bic(4).cols = (91:120);
ref_biclusters.Bic(5).rows = transpose(201:250); ref_biclusters.Bic(5).cols = (121:150);
ref_biclusters.Bic(6).rows = transpose(251:300); ref_biclusters.Bic(6).cols = (151:180);
ref_biclusters.Bic(7).rows = transpose(301:350); ref_biclusters.Bic(7).cols = (181:210);
ref_biclusters.Bic(8).rows = transpose(351:400); ref_biclusters.Bic(8).cols = (211:240);
ref_biclusters.Bic(9).rows = transpose(401:450); ref_biclusters.Bic(9).cols = (241:270);
ref_biclusters.Bic(10).rows = transpose(451:500); ref_biclusters.Bic(10).cols = (271:300);    
end
ref_biclusters.nbicluster = nbiclusters;

end

