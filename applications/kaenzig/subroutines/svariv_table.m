SVARIV_FVD_med = min(1,max(0,SVARIV_FVD.estim));
SVARIV_FVD_lb  = min(1,max(0,SVARIV_FVD.ci.lower));
SVARIV_FVD_ub  = min(1,max(0,SVARIV_FVD.ci.upper));

nvar = size(SVARIV_FVD_med,2);

nvarbase = 5;
FID = fopen('tables/svariv_kaenzig.tex','w');
fprintf(FID, strcat('\\begin{tabular}{l',repmat('c',1,6),'}\\toprule\\midrule  \n'));  
fprintf(FID,'\\multicolumn{6}{l}{\\textit{Global variables and exchange rates:}} \\\\ \n');  
nameString = [ ];
for ii = 1:nvarbase-1
    nameString = [nameString ' %s &'];
end
nameString = [nameString ' %s '];
fprintf(FID,strcat(' &  ',nameString,'  \\\\ \\midrule  \n'),data.endo_vars{1:nvarbase});
numString = ' %1.0f & ';
for ii = 1:nvarbase-1
    numString = [numString '%8.2f & '];
end
numString = [numString '%8.2f'];
numStringBands = [];
for ii = 1:nvarbase-1
    numStringBands = [numStringBands ' [%1.2f,%8.2f] &'];
end
numStringBands = [numStringBands ' [%1.2f,%8.2f]'];
for kk=[0 12 24 48] %48
    kk2=kk+1;
    if kk<36
        fprintf(FID,strcat(numString,'  \\\\  \n'),kk,SVARIV_FVD_med(kk2,1:nvarbase));
        fprintf(FID, strcat(' & ',numStringBands,' \\\\  \n'),vec([SVARIV_FVD_lb(kk2,1:nvarbase)' SVARIV_FVD_ub(kk2,1:nvarbase)']'));
    else
        fprintf(FID,strcat(numString,'  \\\\  \n'),kk,SVARIV_FVD_med(kk2,1:nvarbase));
        fprintf(FID, strcat(' & ',numStringBands,' \\\\[2ex] \\midrule  \n'),vec([SVARIV_FVD_lb(kk2,1:nvarbase)' SVARIV_FVD_ub(kk2,1:nvarbase)']'));
    end
end
fprintf(FID,'\\multicolumn{6}{l}{\\textit{U.S. variables:}} \\\\ \n'); 
nameString = [ ];
for ii = 1:nvar-nvarbase-1
    nameString = [nameString ' %s &'];
end
nameString = [nameString ' %s '];
fprintf(FID,strcat(' &  ',nameString,'  \\\\ \\midrule  \n'),data.endo_vars{nvarbase+1:end});
numString = ' %1.0f & ';
for ii = 1:nvar-nvarbase-1
    numString = [numString '%8.2f & '];
end
numString = [numString '%8.2f'];
numStringBands = [];
for ii = 1:nvar-nvarbase-1
    numStringBands = [numStringBands ' [%1.2f,%8.2f] &'];
end
numStringBands = [numStringBands ' [%1.2f,%8.2f]'];
for kk=[0 12 24 48] 
    kk2=kk+1;
    fprintf(FID,strcat(numString,'  \\\\  \n'),kk,SVARIV_FVD_med(kk2,nvarbase+1:end));
    fprintf(FID, strcat(' & ',numStringBands,' \\\\  \n'),vec([SVARIV_FVD_lb(kk2,nvarbase+1:end)' SVARIV_FVD_ub(kk2,nvarbase+1:end)']'));

end
fprintf(FID, '\\midrule\\bottomrule \n');
fprintf(FID, '\\end{tabular}\n');
fclose(FID);

clear ans nvarbase FID ii kk kk2 nameString numString numStringBands nvar SVARIV_FVD_lb SVARIV_FVD_med SVARIV_FVD_ub