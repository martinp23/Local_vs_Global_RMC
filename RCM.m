function [write_storex,write_storey,write_storez] = ...
    RCM(xyz,conn,idx,separation,special_rotation)
% Calculation of the RCGFx,y,z for the selected atoms.

M_spec = give_RCM(xyz,conn,separation,special_rotation);

% split the connectivity matrix into two
% loops above/below the molecule. M_spec
% stores connectivity matrix of xyz  
% coordinates definintg the RCM.

count = 1;
store = [];

write_storex = cell(size(idx,2),1);
write_storey = cell(size(idx,2),1);
write_storez = cell(size(idx,2),1);

for nuclei = 1:size(idx,2)
    a = idx(:,nuclei);
    a = a(find(a~=0));
    Mato=[];
    Matox=[];
    Matoy=[];
    cnt=1;
    for i = 1:size(a,1)
        x = xyz(a(i,1),1);
        y = xyz(a(i,1),2);
        z = xyz(a(i,1),3);
        Bz = zeros(size(M_spec,1),1);
        count = 1;
        Bx = zeros(size(M_spec,1),1);
        By = zeros(size(M_spec,1),1);
        Bz = zeros(size(M_spec,1),1);

        % This section sums contribution from every
        % single RCM block.
        for j = 1:size(M_spec,1);
            a1 = M_spec(j,5)-M_spec(j,2);
            a2 = M_spec(j,2);
            b1 = M_spec(j,6)-M_spec(j,3);
            b2 = M_spec(j,3);
            c1 = M_spec(j,7)-M_spec(j,4);
            c2 = M_spec(j,4);
            
            t = 1;
            [Bxt,Byt,Bzt] = unit3d(a1,a2,b1,b2,c1,c2,t,x,y,z);
            A_Bx = M_spec(j,1)*Bxt;
            A_By = M_spec(j,1)*Byt;
            A_Bz = M_spec(j,1)*Bzt;
            t = 0;
            [Bxt,Byt,Bzt] = unit3d(a1,a2,b1,b2,c1,c2,t,x,y,z);
            B_Bx = M_spec(j,1)*Bxt;
            B_By = M_spec(j,1)*Byt;
            B_Bz = M_spec(j,1)*Bzt;
            
            Bx(count,1) = Bx(count,1) + A_Bx - B_Bx;
            By(count,1) = By(count,1) + A_By - B_By;
            Bz(count,1) = Bz(count,1) + A_Bz - B_Bz;
            count=count+1;
        end
        Mato(cnt,1) = sum(Bz);
        Matox(cnt,1) = sum(Bx);
        Matoy(cnt,1) = sum(By);
        cnt=cnt+1;
    end
    write_storex{nuclei,1} = Matox;
    write_storey{nuclei,1} = Matoy;
    write_storez{nuclei,1} = Mato;
end
end












