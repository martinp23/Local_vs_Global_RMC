function [Ax,Ay,Az] = area_finder_local(xyz,conn,geomTXT,ring_size);
% This function estimates normalised
% net cross-section areas of the RCM
% relatively to the x,y,z axes.
% Only used for sinle porphyrin subunits
%count = 1
%conn = conn_local_sub{count,1}
conn_temp = conn(:,2:3);
idx_temp = unique(conn_temp(:));

%idx_temp2 = find(conn(:,1) == 1);
%conn(idx_temp2,:)
%return

if ring_size == 20 || ring_size == 28;

distMat = squareform(pdist(xyz(idx_temp,:))); % distance matrix
adjMat = distMat < 1.75;
adjMat = adjMat - diag(diag(adjMat));
G = graph(adjMat);
cycles = allcycles(G);
size_graphs = zeros(size(cycles,1),1);
for i = 1:size(cycles,1)
    size_graphs(i,1) = size(cycles{i,1},2);
end

Idx20 = find(size_graphs == 20);

idx_ring = idx_temp(cycles{Idx20,1}');
xyz_temp = xyz(idx_ring,:);

idx_store = [];
for i = 1:size(conn,1)

idx_store(i,1) = ...
    ismember(conn(i,2),idx_ring) & ismember(conn(i,3),idx_ring) ;
end
idx_temp3 = find(idx_store == 1);

conn20 = conn(idx_temp3,:);

reor_conn20 = zeros(size(conn20));
reor_conn20(1,:) = conn20(1,:);

num_1 = reor_conn20(1,3);
count = 2;

while num_1 ~= reor_conn20(1,2)
    num_1 = reor_conn20(count-1,3);
    num_1 = find(conn20(:,2) == num_1);
    reor_conn20(count,:) = conn20(num_1,:);
    count = count + 1;
    num_1 = reor_conn20(count-1,3);
end

%scatter3(xyz_temp(:,1),xyz_temp(:,2),xyz_temp(:,3))



X = xyz(reor_conn20(:,2),1);
Y = xyz(reor_conn20(:,2),2);
Z = xyz(reor_conn20(:,2),3);


elseif ring_size == 16

distMat = squareform(pdist(xyz(idx_temp,:))); % distance matrix
adjMat = distMat < 1.75;
adjMat = adjMat - diag(diag(adjMat));
G = graph(adjMat);
cycles = allcycles(G);
size_graphs = zeros(size(cycles,1),1);
for i = 1:size(cycles,1)
    size_graphs(i,1) = size(cycles{i,1},2);
end

Idx16= find(size_graphs == 16);

idx_ring = idx_temp(cycles{Idx16,1}');
xyz_temp = xyz(idx_ring,:);

idx_store = [];
for i = 1:size(conn,1)

idx_store(i,1) = ...
    ismember(conn(i,2),idx_ring) & ismember(conn(i,3),idx_ring) ;
end
idx_temp3 = find(idx_store == 1);

conn16 = conn(idx_temp3,:);

reor_conn16 = zeros(size(conn16));
reor_conn16(1,:) = conn16(1,:);

num_1 = reor_conn16(1,3);
count = 2;

while num_1 ~= reor_conn16(1,2)
    num_1 = reor_conn16(count-1,3);
    num_1 = find(conn16(:,2) == num_1);
    reor_conn16(count,:) = conn16(num_1,:);
    count = count + 1;
    num_1 = reor_conn16(count-1,3);
end

X = xyz(reor_conn16(:,2),1);
Y = xyz(reor_conn16(:,2),2);
Z = xyz(reor_conn16(:,2),3);

end
% disp(ring_size)

Ax = polyarea(Y,Z);
Ay = polyarea(X,Z);
Az = polyarea(X,Y);


if ~ispolycw(Z,Y)
    Ax = -1*Ax;
end

if ~ispolycw(X,Z)
    Ay = -1*Ay;
end

if ~ispolycw(Y,X)
    Az = -1*Az;
end

scale_factor = sqrt(Ax^2+Ay^2+Az^2);

Ax = Ax/scale_factor;
Ay = Ay/scale_factor;
Az = Az/scale_factor;

end


