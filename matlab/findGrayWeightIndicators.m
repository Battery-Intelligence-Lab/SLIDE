function gray = findGrayWeightIndicators(nx,ny,nz)

gray.nx = nx;
gray.ny = ny;
gray.nz = nz;

% gray.nx = data.segm.P.n;
% gray.ny = data.segm.SOC.n;
% gray.nz = data.segm.T.n;

% For 3,4,5 it generates extra bit due to residuals. Find a better way. 

% [~,~,ZZZ] = createRecursiveGrayCode3(4);
% ZZZ = double(ZZZ);
% ZZZ2 = cell(9,1);
% for i =1:9
%     ZZZ2{i} = ZZZ(i,:);
%     
% end


gray.code = createGrayCodeEnum(gray.nx,gray.ny,gray.nz);

gray.nt = length(gray.code{1});

gray.delta         = repmat({false(gray.nx+1,gray.ny+1,gray.nz+1)},gray.nt,1);
gray.one_min_delta = repmat({false(gray.nx+1,gray.ny+1,gray.nz+1)},gray.nt,1);

for iGray =1:(gray.nx+1)
    for jGray =1:(gray.ny+1)
        for kGray =1:(gray.nz+1)
            
            gray.possib = [iGray,       jGray,      kGray;...
                           iGray-1,     jGray,      kGray;...
                           iGray,       jGray-1,    kGray;...
                           iGray,       jGray,      kGray-1;...
                           iGray-1,     jGray-1,    kGray;...
                           iGray,       jGray-1,    kGray-1;...
                           iGray-1,     jGray,      kGray-1;...
                           iGray-1,     jGray-1,    kGray-1];
           
            gray.possib(:,1) = max(min(gray.possib(:,1),gray.nx),1);
            gray.possib(:,2) = max(min(gray.possib(:,2),gray.ny),1);
            gray.possib(:,3) = max(min(gray.possib(:,3),gray.nz),1);

            gray.non_changing = gray.code{gray.possib(1,1),gray.possib(1,2),gray.possib(1,3)};
            
            gray.outcome      = true(size(gray.non_changing));
            
            for iNon = 2:size(gray.possib,1)
                temp = gray.non_changing==gray.code{gray.possib(iNon,1),gray.possib(iNon,2),gray.possib(iNon,3)};
                gray.outcome = temp&gray.outcome;
            end
            
            for iOut = 1:gray.nt
               if(gray.outcome(iOut)) %if it is a non-changing variable. 
                   if(gray.non_changing(iOut)) %then if it is 1 it is delta
                      gray.delta{iOut}(iGray,jGray,kGray) = true;
                   else %else it is 1-delta 
                      gray.one_min_delta{iOut}(iGray,jGray,kGray) = true; 
                   end
               end
            end
        end
    end
end
end