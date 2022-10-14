clear
close all

run_RCGF = 1;
show_fig = 1;
show_dif_colors = 1;
special_rotation = 0;

separation=0.7;

% for local vs global swap change later the conn = conn_glob to conn = conn_loc_SS
% setup of the script
show_connectivity = 1;
show_heads = 1;
show_numbering = 0;
show_hydrogens = 1;
print_figure = 0;

xyz_file = 'calc.xyz' ;
xyz_path='.\octaphyrin\';

xyz_crude = readcell( (fullfile(xyz_path,xyz_file)),...
    'FileType', 'text');

if isnumeric(xyz_crude{2,2})
    init_idx = 2;
elseif isnumeric(xyz_crude{3,2})
    init_idx = 3;
else
    disp('problem')
    return
end

geomTXT = string(xyz_crude(init_idx:end,1));
geomNO = cell2mat(xyz_crude(init_idx:end,2:4));

atoms_list = csvread(strcat(xyz_path,xyz_file(1:end-4),'_atoms.csv'));
conn_glob = csvread(strcat(xyz_path,xyz_file(1:end-4),'_connectivity_global.csv'));
% conn_loc_S = csvread(strcat(xyz_path,xyz_file(1:end-4),'_connectivity_local_S.csv'));
conn_loc_SS = csvread(strcat(xyz_path,xyz_file(1:end-4),'_connectivity_local_SS.csv'));
conn_glob_ave = csvread(strcat(xyz_path,xyz_file(1:end-4),'_connectivity_global_ave.csv'));
conn_loc_SS_ave = csvread(strcat(xyz_path,xyz_file(1:end-4),'_connectivity_local_SS_ave.csv'));

% conn = conn_loc_SS_ave;
conn = conn_glob; % global ring
% conn = conn_loc_SS; % local ring


for setting_par = 1;

    max_special = ceil(max(max(geomNO)))+0.5;
    min_special = floor(min(min(geomNO)))-0.5;
    max_lims = ceil(max(geomNO))+0.5;
    min_lims = floor(min(geomNO))-0.5;

    figsize = 2;
    size_faktor = 15/max_special;
    % this adjusts pop-up window in the centre of the screen
    if show_fig;
        temp_handle = get(0,'ScreenSize');
        fig = figure('units','centimeters','position',...
            [temp_handle(3:4)/2-9*figsize/2,9*figsize,9*figsize]);
        ax1 = axes;
        movegui('center');
        view([83 21]);
    end

    BT = 0.4*figsize*size_faktor; % bond thickness
    BT2 = 0.3*figsize*size_faktor; % C-H bond thickness
    SS = 17*figsize*size_faktor; % sizes of atom
    SS2 = 100*figsize*size_faktor; % sizes of spectator atom
    ringT = 4*figsize*size_faktor; % ring current thickness
    trp = 0.5; % transparency
    spec_trp = 0.5; % spectator atoms transparency
    text_font_size = 8*figsize*size_faktor; % font size for the text

    % defining the xyz coordinates of individual atoms
    c = geomNO(find(geomTXT== "C"),:);
    non_h = geomNO(find(geomTXT~= "H"),:);
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

end
distMat = squareform(pdist(geomNO)); % distance matrix
adjMat = distMat < 2.0; % cutoff for drawing bonds
if show_fig
    for i=1:size(adjMat,1)
        for j=1:size(adjMat,2)
            if i ~= j && adjMat(i,j) == 1
                if geomTXT(i) == "H" && geomTXT(j) == "H"
                else
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
                end
            else
                % do nothing
            end
        end
    end
end

xyz=geomNO;
M_spec = give_RCM(xyz,conn,separation,special_rotation); % returns RCGF(x,y,z)

if show_fig;
    for scat1 = 1;
        % drawing atoms
        s1 = scatter3(c(:,1),c(:,2),c(:,3),SS,'filled');
        s1.MarkerFaceColor = [0 0 0];
        s1.MarkerEdgeColor = [0 0 0];
        hold on
        s2 = scatter3(s(:,1),s(:,2),s(:,3),SS,'filled');
        s2.MarkerFaceColor = [1 165/255 0];
        s2.MarkerEdgeColor = [0 0 0];
        hold on
        if show_hydrogens
            s3 = scatter3(h(:,1),h(:,2),h(:,3),SS,'filled');
            s3.MarkerFaceColor = 0.8*[1 1 1];
            s3.MarkerEdgeColor = [0 0 0];
            hold on
        else
        end

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
    ring_list = [1,0,0;
        0,0.5,0;
        0,0,1];
    conn_color_map = [...
        252,12,12;...
        0,17,255;...
        204,0,255;...
        95,252,12;...
        102,102,102]/255;
    j=1;
    for i = 1:size(M_spec,1)
        r(i) = plot3(...
            [M_spec(i,2);M_spec(i,5)],...
            [M_spec(i,3);M_spec(i,6)],...
            [M_spec(i,4);M_spec(i,7)]...
            );

        ring_current_weights = flipud(unique(M_spec(:,1)));
        %     for i = 1:size(conn,1)
        if show_dif_colors
            r(i).Color = conn_color_map(...
                find(ring_current_weights ==  M_spec(i,1)),:);
        else
            r(i).Color = ring_list(j,:);
        end
        r(i).LineStyle = '-';
        r(i).LineWidth = 2;
        r(i).Color(4)=1;
        hold on
    end
    for conn_set = 1;
        if show_connectivity % show path (the path is only single loop)
            conn_color_map = [...
                252,12,12;...
                0,17,255;...
                204,0,255;...
                95,252,12;...
                102,102,102]/255;
            ring_current_weights = flipud(unique(conn(:,1)));

            if 0
                for i = 1:size(conn,1)
                    r(i) = plot3(...
                        [geomNO(conn(i,2),1);geomNO(conn(i,3),1)],...
                        [geomNO(conn(i,2),2);geomNO(conn(i,3),2)],...
                        [geomNO(conn(i,2),3);geomNO(conn(i,3),3)]);
                    r(i).Color = conn_color_map(...
                        find(ring_current_weights ==  conn(i,1)),:);
                    r(i).LineStyle = '-';
                    r(i).LineWidth = abs(conn(i,1))*ringT;
                    r(i).Color(4)=trp;
                    hold on
                end
            end

            if ~show_heads
                leg = legend([r(1:size(ring_current_weights,1))]',...
                    num2str(ring_current_weights,5));
                leg.Position(1:2) = leg.Position(1:2) - 0.05;
            end

            if show_heads
                start = M_spec(:,2:4);
                stop = M_spec(:,5:7);

                %                 start = geomNO(conn(:,2),1:3);
                %                 stop = geomNO(conn(:,3),1:3);
                head_clr = zeros(size(r,2),3);
                for i = 1:size(r,2)
                    head_clr(i,1:3) = r(i).Color;
                end
                dvec=stop-start;
                dis=sqrt(sum(dvec.^2,2));
                hv=min(dis)*0.35;
                cosrang=acos(dvec(:,3)./dis)*180/pi;
                nvec=[-dvec(:,2) dvec(:,1) zeros(size(dis))];
                hheads=[];
                hhgrd=[];
                pv=dis-hv;
                for i=1:length(dis)
                    [xi,yi,zi] = cylinder([tan(30/180*pi),0],10);
                    xi=xi*hv;yi=yi*hv;zi=zi*hv+pv(i);
                    [rx,ry,rz] = rotatedata(xi,yi,zi,nvec(i,:),cosrang(i),[0,0,0]);
                    cx=start(i,1)+rx;cy=start(i,2)+ry;cz=start(i,3)+rz;
                    hheads(i)=surf(cx,cy,cz,'edgecolor','none','facecolor',head_clr(i,:)*0.6);
                end
            end
        end

    end

    for fig_set = 1;

        set(gca,'xtick',[],'ytick',[],'ztick',[])
        axis equal
        view([-8 24])
        box off
        axis tight
        set(gca,'visible','off')
        hold off
        axis 'equal'
        axis tight
        set(gca,'xtick',[],'ytick',[])

        ax = gca;
        ax.XLim = [min_lims(1) max_lims(1)];
        ax.YLim = [min_lims(2) max_lims(2)];
        ax.ZLim = [min_lims(3) max_lims(3)];

        axis vis3d
        axis off
        ax.Position = ax.OuterPosition - ax.TightInset...
            *[-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1];

        view([0 90])

    end
end

if run_RCGF
    idx = atoms_list;
    [write_storex,write_storey,write_storez] = ...
        RCM(xyz,conn,idx,separation,special_rotation); % returns RCGF(x,y,z)

    [Ax,Ay,Az] = area_finder(xyz,conn); % returns cross-section areas

    x = zeros(size(write_storez,1),1);
    for i = 1:size(write_storez,1)
        x(i,1) = ...
            [Ax*mean(write_storex{i,1},1)+...
            Ay*mean(write_storey{i,1},1)+...
            Az*mean(write_storez{i,1},1)]/3;
    end

    answer_list = "The averaged RCGF for the atoms are:";
    answer_list_console = answer_list;

    for i = 1:size(x,1)
        answer_list_console = [answer_list_console;
            strcat("atoms no. ",num2str(i),...
            " = ", string(x(i,1))," ppm T / nA")];
    end

    opts.Interpreter = 'tex';
    opts.Default = 'OK and exit';

    disp(answer_list_console);

end


