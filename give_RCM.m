function M_spec = give_RCM(xyz,conn,separation,special_rotation)
% Calculation of the RCGFx,y,z for the selected atoms.

m = size(conn,1);
offset_cell = cell(size(conn,1),2);

if separation == 0;
M_spec = zeros(m,7);
for i = 1:m
    M_spec(i,1:7) = [conn(i,1),...
        xyz(conn(i,2),1:3),...
        xyz(conn(i,3),1:3)];
end

else
M_spec = zeros(2*m,7);
for i=1:m
    for j=1:2
        depo_okoli = [];
        count = 1;
        r_temp = xyz(conn(i,j+1),1:3);
        for k=1:size(xyz,1)
            r_scan = xyz(k,1:3);
            if norm(r_temp-r_scan) < 2.35 % atoms in vicinity
                depo_okoli(count,1:3) = r_scan;
                count = count + 1;
            else
            end
        end
        % calculates in which direction should be
        % the RCM offset
        A = depo_okoli;
        [~,~,V]=svd(A-mean(A,1),0);
        n0=V(:,end);
        fktr = separation/norm(n0);
        n = n0'*fktr; % n is an offset vector. 
        % arbitrary orientation (if up or down)

%         if 0%size(depo_okoli,1) <= 3
%             n = n*0.0001;
%         end

        if special_rotation
           n_tot = (r_temp-mean(xyz));
            n = separation*n_tot/norm(n_tot);
            
        end

        offset_cell{i,j} = n;

    end
    
    % decide orientation within the segment
    X1 = xyz(conn(i,j),1:3);
    X2 = xyz(conn(i,j+1),1:3);
    off_1 = offset_cell{i,j-1}; 
    off_2 = offset_cell{i,j}; 
    dis1 = norm(X1+off_1-X2-off_2);
    dis2 = norm(X1+off_1-X2+off_2);
    if dis1 > dis2
        offset_cell{i,j} = -1*n;
    elseif dis1 == dis2
        disp('error, the path cannot be found unambigously');
    end
end

M_spec = zeros(2*m,7);
for i = 1:m
    M_spec(i,1:7) = [0.5*conn(i,1),...
        xyz(conn(i,2),1:3)+offset_cell{i,1},...
        xyz(conn(i,3),1:3)+offset_cell{i,2}];
    M_spec(i+m,1:7) = [0.5*conn(i,1),...
        xyz(conn(i,2),1:3)-offset_cell{i,1},...
        xyz(conn(i,3),1:3)-offset_cell{i,2}];
end



end





end





