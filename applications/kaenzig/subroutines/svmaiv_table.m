SVMAIV_FVD_IS_lb = min(1,max(0,bounds.estim.lower.FVR));
SVMAIV_FVD_IS_ub = min(1,max(0,bounds.estim.upper.FVR));
SVMAIV_FVD_CI_lb = min(1,max(0,bounds.ci.lower.FVR));
SVMAIV_FVD_CI_ub = min(1,max(0,bounds.ci.upper.FVR));
nvar = size(SVMAIV_FVD_IS_lb,2);

nvarbase = 5;
FID = fopen('tables/svmaiv_kaenzig.tex','w');
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
numStringBands1 = [];
for ii = 1:nvarbase-1
    numStringBands1 = [numStringBands1 ' %1.2f,%8.2f &'];
end
numStringBands1 = [numStringBands1 ' %1.2f,%8.2f'];
numStringBands2 = [];
for ii = 1:nvarbase-1
    numStringBands2 = [numStringBands2 ' [%1.2f,%8.2f] &'];
end
numStringBands2 = [numStringBands2 ' [%1.2f,%8.2f]'];
for kk=[0 12 24 48] %48
    kk2=kk+1;
    if kk<36
        fprintf(FID,strcat(' %1.0f '),kk);
        fprintf(FID, strcat(' & ',numStringBands1,' \\\\  \n'),vec([SVMAIV_FVD_IS_lb(kk2,1:nvarbase)' SVMAIV_FVD_IS_ub(kk2,1:nvarbase)']'));
        fprintf(FID, strcat(' & ',numStringBands2,' \\\\  \n'),vec([SVMAIV_FVD_CI_lb(kk2,1:nvarbase)' SVMAIV_FVD_CI_ub(kk2,1:nvarbase)']'));
    else
        fprintf(FID,strcat(' %1.0f '),kk);
        fprintf(FID, strcat(' & ',numStringBands1,' \\\\  \n'),vec([SVMAIV_FVD_IS_lb(kk2,1:nvarbase)' SVMAIV_FVD_IS_ub(kk2,1:nvarbase)']'));
        fprintf(FID, strcat(' & ',numStringBands2,' \\\\[2ex] \\midrule  \n'),vec([SVMAIV_FVD_CI_lb(kk2,1:nvarbase)' SVMAIV_FVD_CI_ub(kk2,1:nvarbase)']'));
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
numStringBands2 = [];
for ii = 1:nvar-nvarbase-1
    numStringBands2 = [numStringBands2 ' [%1.2f,%8.2f] &'];
end
numStringBands2 = [numStringBands2 ' [%1.2f,%8.2f]'];
for kk=[0 12 24 48] 
    kk2=kk+1;
    fprintf(FID,strcat(' %1.0f '),kk);
    fprintf(FID, strcat(' & ',numStringBands1,' \\\\  \n'),vec([SVMAIV_FVD_IS_lb(kk2,nvarbase+1:end)' SVMAIV_FVD_IS_ub(kk2,nvarbase+1:end)']'));
    fprintf(FID, strcat(' & ',numStringBands2,' \\\\  \n'),vec([SVMAIV_FVD_CI_lb(kk2,nvarbase+1:end)' SVMAIV_FVD_CI_ub(kk2,nvarbase+1:end)']'));

end
fprintf(FID, '\\midrule\\bottomrule \n');
fprintf(FID, '\\end{tabular}\n');
fclose(FID);

clear ans nvarbase FID ii kk kk2 nameString numString numStringBands1 numStringBands2 nvar ...
    SVMAIV_FVD_IS_lb SVMAIV_FVD_IS_ub SVMAIV_FVD_CI_lb SVMAIV_FVD_CI_ub