function x=RCM_anal(conn,idx,geomTXT,geomNO,...
    separation,special_rotation,conn_file)


xyz=geomNO;

if contains(conn_file,'local')
    
step =size(conn,1)/6;
disp(step)
count = 1;
for i = 1:6
    conn_temp = conn((1+(i-1)*step):((i)*step),:);
[write_storex,write_storey,write_storez] = ...
    RCM(xyz,conn_temp,idx,separation,special_rotation); % returns RCGF(x,y,z)

[Ax,Ay,Az] = area_finder_local(xyz,conn_temp,geomTXT,step);
RCM_local{count,1} = [write_storex,write_storey,write_storez] ; 
RCM_local{count,2} = [Ax,Ay,Az];
count = count + 1;
end

x_loc = zeros(size(RCM_local{1,1},1),6);
for j = 1:6
for i = 1:size(RCM_local{j,1},1)
    x_loc(i,j) = ...
        [RCM_local{j,2}(1,1)*mean(RCM_local{j,1}{i,1},1)+...
        RCM_local{j,2}(1,2)*mean(RCM_local{j,1}{i,2},1)+...
        RCM_local{j,2}(1,3)*mean(RCM_local{j,1}{i,3},1)]/3;
end
end

x_loc_sum = sum(x_loc,2);

x=x_loc_sum;
% disp(x_loc)
else

    [write_storex,write_storey,write_storez] = ...
        RCM(xyz,conn,idx,separation,special_rotation); % returns RCGF(x,y,z)
    [Ax,Ay,Az] = area_finder(xyz,conn); % returns cross-section areas

    RCM_global{1,1} = [write_storex,write_storey,write_storez] ;
    RCM_global{1,2} = [Ax,Ay,Az];

    x_glob = zeros(size(RCM_global{1,1},1),1);
    for i = 1:size(RCM_global{1,1},1)
        x_glob(i,1) = ...
            [RCM_global{1,2}(1,1)*mean(RCM_global{1,1}{i,1},1)+...
            RCM_global{1,2}(1,2)*mean(RCM_global{1,1}{i,2},1)+...
            RCM_global{1,2}(1,3)*mean(RCM_global{1,1}{i,3},1)]/3;
    end
    x=x_glob;
    
end

% geometry_plot(geomTXT,geomNO,xyz_file,idx,...
%     conn,show_ring_path,show_spectator_atoms); % geometry visualisation

% [continue_fitting,x] = output_box(Ax,Ay,Az,...
%     write_storex,write_storey,write_storez);


% if 0% continue_fitting
%     fit_and_plot(Ax,Ay,Az,...
%         write_storex,write_storey,write_storez);
% else
% end

end


