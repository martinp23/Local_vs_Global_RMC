function [write_storex,write_storey,write_storez] = ...
    RCM_for_3d(xyz,conn,list,separation,special_rotation)
% Calculation of the RCGFx,y,z for the selected atoms.

M_spec = give_RCM(xyz,conn,separation,special_rotation);

% split the connectivity matrix into two
% loops above/below the molecule. M_spec
% stores connectivity matrix of xyz
% coordinates definintg the RCM.

write_storex = zeros(size(list,1),1);
write_storey = zeros(size(list,1),1);
write_storez = zeros(size(list,1),1);

cnt=1;
for i = 1:size(list,1)
    x=list(i,1);
    y=list(i,2);
    z=list(i,3);

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
    write_storez(cnt,1) = sum(Bz);
    write_storex(cnt,1) = sum(Bx);
    write_storey(cnt,1) = sum(By);
    cnt=cnt+1;
end

end












