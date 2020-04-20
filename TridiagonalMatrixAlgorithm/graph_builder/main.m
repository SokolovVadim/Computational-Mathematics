% graph builder

% input values
fileID = fopen('data.txt','r');
formatSpec = '%f %f %f';
sizeA = [3 Inf];

A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
A = A';

y = Graph_builder(A);



