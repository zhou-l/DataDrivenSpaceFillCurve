function LTfull = fullSFC(LTLvls)
% reconstruct full SFC from hierarcical SFC
    LTfull = [];
    for i = 1:length(LTLvls)
         LTfull  = [LTfull; repmat(LTLvls(i,1),LTLvls(i,2)*LTLvls(i,2),1)];
        
%         LTfull  = [LTfull; repmat(LTLvls(i,1),LTLvls(i,2),1)];
    end
    maxLT = max(LTfull);
    disp(maxLT);
end