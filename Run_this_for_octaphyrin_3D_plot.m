close all
clear

xyz_path='.\octaphyrin\';
xyz_file='calc.xyz';

manual_setup = 1;
special_rotation=0;
separation = 0.7; 

dim = [30,30,30,0.5];

for get_init_data = 1;

    if ~manual_setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    msg1 = 'Select xyz file';
    disp(msg1);
    [xyz_file,xyz_path] = uigetfile('*.xyz',msg1);

    if isequal(xyz_file,0)
        disp('No xyz file provided');
        return
    else
        %     disp(['User selected ', fullfile(xyz_path,xyz_file)]);
        xyz_crude = readcell( (fullfile(xyz_path,xyz_file)),...
            'FileType', 'text');
        if isnumeric(xyz_crude{2,2})
            init_idx = 2;
        elseif isnumeric(xyz_crude{3,2})
            init_idx = 3;
        else
            disp('hey, problem')
            return
        end

        if isnumeric(xyz_crude{init_idx,1})
            for i = init_idx:size(xyz_crude,1)
                if xyz_crude{i,1} == 1
                    xyz_crude{i,1} = "H";
                elseif xyz_crude{i,1} == 6
                    xyz_crude{i,1} = "C";
                elseif xyz_crude{i,1} == 7
                    xyz_crude{i,1} = "N";
                elseif xyz_crude{i,1} == 30
                    xyz_crude{i,1} = "Zn";
                else
                    disp('heeeelp')
                end
            end
        end
        geomTXT = string(xyz_crude(init_idx:end,1));
        geomNO = cell2mat(xyz_crude(init_idx:end,2:4));
        xyz=geomNO;
    end

    else
        xyz_crude = readcell( (fullfile(xyz_path,xyz_file)),...
            'FileType', 'text');
        if isnumeric(xyz_crude{2,2})
            init_idx = 2;
        elseif isnumeric(xyz_crude{3,2})
            init_idx = 3;
        else
            disp('hey, problem')
            return
        end
        if isnumeric(xyz_crude{init_idx,1})
            for i = init_idx:size(xyz_crude,1)
                if xyz_crude{i,1} == 1
                    xyz_crude{i,1} = "H";
                elseif xyz_crude{i,1} == 6
                    xyz_crude{i,1} = "C";
                elseif xyz_crude{i,1} == 7
                    xyz_crude{i,1} = "N";
                elseif xyz_crude{i,1} == 30
                    xyz_crude{i,1} = "Zn";
                else
                    disp('heeeelp')
                end
            end
        end
        geomTXT = string(xyz_crude(init_idx:end,1));
        geomNO = cell2mat(xyz_crude(init_idx:end,2:4));
        xyz=geomNO;

    end

    atoms_list = csvread([xyz_path,xyz_file(1:end-4),'_atoms.csv']);
    conn_glob = csvread([xyz_path,xyz_file(1:end-4),'_connectivity_global.csv']);
%     conn_loc_S = csvread([xyz_path,xyz_file(1:end-4),'_connectivity_local_S.csv']);
    conn_loc_SS = csvread([xyz_path,xyz_file(1:end-4),'_connectivity_local_SS.csv']);
    conn_glob_ave = csvread([xyz_path,xyz_file(1:end-4),'_connectivity_global_ave.csv']);
    conn_loc_SS_ave = csvread([xyz_path,xyz_file(1:end-4),'_connectivity_local_SS_ave.csv']);

end

% conn = conn_glob_ave;
% conn = conn_glob;
% conn = conn_loc_SS;
conn = conn_loc_SS;
%%%% here you swap what you want to visualise



list=zeros((dim(1)/dim(4)+1)*(dim(2)/dim(4)+1)*(dim(3)/dim(4)+1),3);
count=0;
for k=-(dim(3)/2):dim(4):(dim(3)/2)
    for j=-(dim(2)/2):dim(4):(dim(2)/2)
        for i=-(dim(1)/2):dim(4):(dim(1)/2)
            count=count+1;
            list(count,1)=i;
            list(count,2)=j;
            list(count,3)=k;
        end
    end
end

[write_storex,write_storey,write_storez] = ...
    RCM_for_3d(xyz,conn,list,separation,special_rotation); % returns RCGF(x,y,z)

[Ax,Ay,Az] = area_finder(xyz,conn); % returns cross-section areas

x = zeros(size(write_storez,1),1);
for i = 1:size(write_storez,1)
    x(i,1) = ...
        [Ax*mean(write_storex(i,1),1)+...
        Ay*mean(write_storey(i,1),1)+...
        Az*mean(write_storez(i,1),1)]/3;
end

x=-x;

vel = round(size(x,1)^(1/3));
V = zeros(vel,vel,vel);
b=x;

for page = 0:(vel-1);
    c = b((page*(vel^2)+1):((page+1)*(vel^2)));
    %     V(:,:,page+1) = -reshape(c,vel,vel);
    V(:,:,page+1) = rot90(rot90(rot90(flipud(reshape(c,vel,vel)))));
%     V(:,:,page+1) = reshape(c,vel,vel);
end


roz = -12.5:dim(4):12.5;
roz = -dim(1)/2:dim(4):dim(1)/2;

scale_faktor = 2.5;

% this adjusts pop-up window in the centre of the screen
temp_handle = get(0,'ScreenSize');
fig = figure('units','centimeters','position',...
    [temp_handle(3:4)/2-4.5,12*scale_faktor,9*scale_faktor]);
movegui('center');
view([0 90]);

meze = 10;

barvicky = [-meze,-meze/2,0,meze/2,meze];

vel_list = size(barvicky,2);


if mod(vel_list,2) == 0
    points = vel_list+1;
else
    points = vel_list;
end

% points=3;
newmap=zeros(points,3);
for i=0:points-1
    newmap(i+1,:)=diverging_map(i/(points-1),[33 102 172]./255, [176 24 43]./255);
end
newmap(ceil(points/2),:) = [0.7 0.7 0.7];
% end

show_0 = 0;

count = 1;
for i = barvicky %-meze:(meze/(points-3)):meze
    if show_0 || i~=0
        s_iso = patch(isosurface(roz,roz,roz,V,i*1e-2));
        s_iso.FaceColor = newmap(count,:);
        s_iso.EdgeColor = 'none';
        s_iso.FaceAlpha = 0.5;
        disp(newmap(count,:))
    end
    count = count + 1;
    hold on
end

axis 'equal';
% camlight
% lighting gouraud
% view(45,30)

% plot 3d stuff

BT = 0.9*scale_faktor; % bond thickness
BT2 = 0.7*scale_faktor; % C-H bond thickness
SS = 17*scale_faktor; % sizes of atom
SS2 = 100*scale_faktor; % sizes of spectator atom
ringT = 6*scale_faktor; % ring current thickness
trp = 0.5; % transparency
spec_trp = 0.5; % spectator atoms transparency
text_font_size = 8*scale_faktor; % font size for the text

% defining the xyz coordinates of individual atoms
c = geomNO(find(geomTXT== "C"),:);
h = geomNO(find(geomTXT== "H"),:);
n = geomNO(find(geomTXT== "N"),:);
s = geomNO(find(geomTXT== "S"),:);
f = geomNO(find(geomTXT== "F"),:);
o = geomNO(find(geomTXT== "O"),:);
zn = geomNO(find(geomTXT== "Zn"),:);
at_rest = geomNO(find(geomTXT ~= "C" & geomTXT ~= "H"...
    & geomTXT ~= "N" & geomTXT ~= "S" ...
    & geomTXT ~= "F" & geomTXT ~= "O" ...
    & geomTXT ~= "Zn"),:);

distMat = squareform(pdist(geomNO)); % distance matrix
adjMat = distMat < 1.75; % cutoff for drawing bonds
adjMatlong = distMat < 2.25; % cutoff for drawing bonds
for i=1:size(adjMat,1)
    for j=1:size(adjMat,2)
        if i ~= j && adjMat(i,j) == 1
            hold on
            g = plot3([geomNO(i,1);geomNO(j,1)],...
                [geomNO(i,2);geomNO(j,2)],...
                [geomNO(i,3);geomNO(j,3)],'k');
            set(g,'LineStyle','-');
            if geomTXT(i,1) == "H" & ...
                    geomTXT(j,1) == "C" || ...
                    geomTXT(i,1) == "C" & ...
                    geomTXT(j,1) == "H"
                set(g,'LineWidth',BT2);
                hold on
            else
                set(g,'LineWidth',BT);
                hold on
            end
        else
            % do nothing
        end
    end
end

for i=1:size(adjMatlong,1)
    for j=1:size(adjMatlong,2)
        if i ~= j && adjMatlong(i,j) == 1

            if geomTXT(i,1) == "S" || ...
                    geomTXT(j,1) == "S"
                hold on
                g = plot3([geomNO(i,1);geomNO(j,1)],...
                    [geomNO(i,2);geomNO(j,2)],...
                    [geomNO(i,3);geomNO(j,3)],'k');
                set(g,'LineWidth',BT2);
                hold on
            end
        end
    end
end

for degen=1 % drawing atoms
    s1 = scatter3(c(:,1),c(:,2),c(:,3),SS,'filled');
    s1.MarkerFaceColor = [0 0 0];
    s1.MarkerEdgeColor = [0 0 0];
    hold on
    s2 = scatter3(s(:,1),s(:,2),s(:,3),SS,'filled');
    s2.MarkerFaceColor = [1 165/255 0];
    s2.MarkerEdgeColor = [0 0 0];
    hold on
    s3 = scatter3(h(:,1),h(:,2),h(:,3),SS,'filled');
    s3.MarkerFaceColor = 0.8*[1 1 1];
    s3.MarkerEdgeColor = [0 0 0];
    hold on
    s4 = scatter3(n(:,1),n(:,2),n(:,3),SS,'filled');
    s4.MarkerFaceColor = [0 0 1];
    s4.MarkerEdgeColor = [0 0 0];
    hold on
    s5 = scatter3(f(:,1),f(:,2),f(:,3),SS,'filled');
    s5.MarkerFaceColor = [178 102 255]/255;
    s5.MarkerEdgeColor = [0 0 0];
    hold on
    s6 = scatter3(o(:,1),o(:,2),o(:,3),SS,'filled');
    s6.MarkerFaceColor = [255 0 0]/255;
    s6.MarkerEdgeColor = [0 0 0];
    hold on
    s7 = scatter3(zn(:,1),zn(:,2),zn(:,3),SS,'filled');
    s7.MarkerFaceColor = [0 155 155]/255;
    s7.MarkerEdgeColor = [0 0 0];
    hold on
    s8 = scatter3(at_rest(:,1),at_rest(:,2),at_rest(:,3),SS,'filled');
    s8.MarkerFaceColor = [0.5 0.5 0.5];
    s8.MarkerEdgeColor = [0 0 0];
    hold on
end

hold off
axis 'equal'
axis tight
axis vis3d
set(gca,'xtick',[],'ytick',[])
% xlim([-15 15])
% ylim([-15 15])
xlim([-dim(1)/2 dim(1)/2])
ylim([-dim(1)/2 dim(1)/2])

axis off
set(gca, 'Position', ...
    get(gca, 'OuterPosition') - ...
    get(gca,'TightInset')...
    *[-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

