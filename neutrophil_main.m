%% neutrophil ca dynamics
% first define contact area growth
p = struct;
p.t = (-100:5000)';
Amin = 0;
t0 = 20;
% t0 = 20/3;
tAc = p.t;
% Ac = 100*ones(size(t));
% dAcdt = zeros(size(t));

rhoIgGSobol = [24.12, 44.42, 1019, 24891, 100, 300];
popTest = true;
densityTest = false;
sobolTest = false;

if sobolTest
    [txtFile,txtPath] = uigetfile('*.txt','Choose parameter file');
    paramFile = fullfile(txtPath,txtFile);
    fid = fopen(paramFile,'r');
    curLine = fgetl(fid);
    lineNum = 1;
    while ischar(fgetl(fid))
        lineNum = lineNum + 1;
    end
    fclose(fid);
    variations = lineNum;
    fid = fopen(paramFile,'r');
    [saveFile,savePath] = uiputfile('*.txt','Choose place to save output');
    saveFileRef = fullfile(savePath,saveFile);
%     fidSave = fopen(saveFile,'w');
elseif densityTest
%     densityVals = logspace(log10(10),log10(50000),500);
%     variations = length(densityVals);
    paramVals = logspace(log10(.12/1.5),log10(.12*1.5),5);
    variations = length(paramVals);
elseif popTest
    variations = 1;%5000;
    xParamVar = zeros(38,variations);
    sigmaVals = [1.1, 1.1 ,1.1 ,1.5, 1.5, 2, 2, 1.5, 2, 1.5,...
        2, 2, 2, 2, 2, 2, 2, 2, 1.5, 1.5,...
        1.5, 1.5, 2, 1.5, 1.5, 1.5, 2, 1.5, 1.5, 2,...
        2, 1.5, 2, 1.5, 1.5,...
        1.2, 1.2, 1.2]';
    sigmaVals = log(sigmaVals).*ones(38,1);
end

% for m = [1,5,6]%1:length(rhoIgGSobol)
% 
%     if sobolTest
%         fid = fopen(paramFile,'r');
%     end
% 
% 
% %set current save file
% saveFile = strcat(saveFileRef(1:end-4),sprintf('_%.1frho.txt',rhoIgGSobol(m)));
% fidSave = fopen(saveFile,'w');
% p.rhoIgG = rhoIgGSobol(m);

% fidSave = fopen(saveFileRef,'w');

% %rhoIgGSobol = [24.12, 44.42, 1019, 24891, 100, 300];
p.rhoIgG = 24891;
expDensities = [.1, 24.11959, 44.42112, 1019.07155, 24891.35052, 100000];
maxAreaExp =  [0, 203.54056, 203.32, 238.89495, 239.2046, 240];
p.Amax = interp1(log(expDensities'),maxAreaExp',log(p.rhoIgG));
% Amax = 80;


condition = 'Standard';
p = neutrophilParam(p,condition);
p0 = p.pVec;

% solve system
options = odeset('RelTol',1e-5,'MaxStep',2,'NonNegative',1:8);
tSpan = [0 400];
p.startTime = tSpan(1);

auc = zeros(variations,1);
peak = zeros(variations,1);
numPks = zeros(variations,1);
ss = zeros(variations,1);
startFlux = zeros(variations,1);
fluxArea = zeros(variations,1);
auc100 = zeros(variations,1);
maxIncrease = zeros(variations,1);
ySaved = cell(variations,1);


% figure
hold on
tic

for j = 1:variations
    if popTest
        sigmaCur = 0*sigmaVals.*randn([38,1]);
        %     sigmaCur = sigmaVals;
        %     sigmaCur(15) = sigmaCur(15) + log(3);
        %     sigmaCur(18) = sigmaCur(18) + log(0.3);
        %     sigmaCur(21) = sigmaCur(21) + log(0.4);
        xParamVar(:,j) = p0 .* exp(sigmaCur);
        p.pVec = xParamVar(:,j);
    elseif sobolTest
        curParamLine = fgetl(fid);
        curParamVals = sscanf(curParamLine,'%f');
        p.pVec = p0 .* curParamVals;
    elseif densityTest
%         p.rhoIgG = densityVals(j);
%         p.pVec(5) = densityVals(j);
%         curDensity = p.pVec(5);
%         Amax = interp1(log(expDensities'),maxAreaExp',log(curDensity));
%         p.pVec(36) = Amax;
%         p.pVec(38) = 20*(Amax/240);
        p.pVec(22) = paramVals(j);
        % Amax = 80;
    end
    p.Ac = p.pVec(36) ./ (1 + exp(-(tAc-p.pVec(37))/p.pVec(38)));
    p.dAcdt = p.pVec(36)*exp(-(tAc-p.pVec(37))/p.pVec(38))./...
        (p.pVec(38)*(1 + exp(-(tAc-p.pVec(37))/p.pVec(38))).^2);
    
    y0 = neutrophilInit(p);
%     initFlag = any(~isreal(y0));
    if any(~isreal(y0))
        if sobolTest
%             pTest = p;
%             pTest.pVec = p0;
%             y0 = neutrophilInit(pTest);
%             initFlag = true;
              y0 = real(y0);
%             fprintf(fidSave,'%f\t%f\t%f\t%f\t%f\t%f\r\n',nan,nan,nan,nan,nan,nan);
%             continue
        else
            error('imaginary initial conditions computed')
        end
    end

    if any(y0<0)
        fprintf(fidSave,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n',0,0,0,y0(1),tSpan(end),0,0,max(p.Ac));
        continue
    end

%     y0(1) = 0.1;
%     y0(2) = 200;
    if strcmp(condition,'Thaps') && p.rhoIgG > 0
        y0(1) = 0.112;
        y0(2) = 10;
%     elseif strcmp(condition,'Thaps')
%         y0(1) = .08;
%         y0(2) = 100;
    end
    if strcmp(condition, 'Ca Free') && p.rhoIgG > 0
        y0(2) = 160;
    end
    
    if strcmp(condition,'BAPTA')
        y0(1) = .092;
%         y0(2) = 147;
    end
%     y0(1) = .08;
%     y0(2) = 200;
%     y0(8) = .65;
    [t,y] = ode15s(@neutrophil,tSpan,y0,options,p);
%     try [t,y] = ode15s(@neutrophil,tSpan,y0,options,p);
%     catch
%         fprintf(fidSave,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n',nan,nan,nan,nan,nan,nan,nan,nan);
%         continue
%     end
%     AcInterp = interp1(p.t,p.Ac,t);
    
    if popTest || densityTest
        if any(~isreal(y(:,1)))
            y = real(y);
        end
        if max(t) < tSpan(2)
            tFinal = [max(t):tSpan(2),tSpan(2)]';
            t = [t; tFinal];
            finalVals = y(end,:);
            for k = size(y,1)+1:length(t)
                y(k,:) = finalVals;
            end
        end
%         if ~mod(j,100)
%             plot(t,y(:,1)*1000)
%             drawnow
%         end
        auc(j) = trapz(t,y(:,1));
        peak(j) = max(y(:,1));
        numPks(j) = length(findpeaks(y(:,1)));
        ss(j) = y(1,1);
        deriv = diff(y(:,1))./diff(t);
        relDeriv = deriv*100/ss(j);
        derivThresh = 2.5*(peak(j)/ss(j))/5;
        curStart = t(find(relDeriv>derivThresh,1,'first'));
        if isempty(curStart)
            curStart = t(end);
            fluxArea(j) = p.Amax;
            auc100(j) = 100*100;
        else
            fluxArea(j) = interp1(p.t,p.Ac,curStart);
            intStart = find(t > curStart,1,'first'); 
            intEnd = find(t < curStart + 100, 1, 'last');
            intTime = [curStart; t(intStart:intEnd); curStart+100];
            try intCa = 1000*[interp1(t,y(:,1),curStart); y(intStart:intEnd,1); interp1(t,y(:,1),curStart+100)];
            catch % I don't think this will happen
                error('Figure out what to do for this case - simulation doesn''t extend 100 s past start of burst')
            end
            auc100(j) = trapz(intTime,intCa);
        end
        startFlux(j) = curStart;
        maxIncrease(j) = max(deriv);
        ySaved{j} = horzcat(t,y(:,1));
        plot(t,1000*y(:,1),'-')
%         yyaxis right
%         plot(t,y(:,2),'-')
%         hold on
%         yyaxis left
%         if startFlux(j) < 50
%             plot([-50; t - startFlux(j)],1000*[y(1,1); y(:,1)])
%         else
%             plot(t-startFlux(j),1000*y(:,1))
%         end
    else
        if any(~isreal(y(:,1)))
            y = real(y);
        end
        %             for k = 1:size(y,1)
        %                 if any(y(k,:) < 0)
        %                     y = y(1:k-1,:);
        %                     t = t(1:k-1);
        %                     break
        %                 end
        %             end
        if max(t) < tSpan(2)
            tFinal = [max(t):tSpan(2),tSpan(2)]';
            t = [t; tFinal];
            finalVals = y(end,:);
            for k = size(y,1)+1:length(t)
                y(k,:) = finalVals;
            end
        end
        if ~mod(j,500)
            plot(t,y(:,1)*1000)
            drawnow
        end
        auc = trapz(t,y(:,1));
        peak = max(y(:,1));
        numPks = length(findpeaks(y(:,1)));
        ss = y(1,1);
        deriv = diff(y(:,1))./diff(t);
        relDeriv = deriv*100/ss;
        derivThresh = 2.5*(peak/ss)/5;
        curStart = t(find(relDeriv>derivThresh,1,'first'));
        if isempty(curStart)
            curStart = t(end);
            fluxArea = p.Amax;
            auc100 = 100*100;
        else
            fluxArea = interp1(p.t,p.Ac,curStart);
            intStart = find(t > curStart,1,'first'); 
            intEnd = find(t < curStart + 100, 1, 'last');
            intTime = [curStart; t(intStart:intEnd); curStart+100];
            try intCa = 1000*[interp1(t,y(:,1),curStart); y(intStart:intEnd,1); interp1(t,y(:,1),curStart+100)];
                auc100 = trapz(intTime,intCa);
            catch % I don't think this will happen
                auc100 = 100*100;
%                 error('Figure out what to do for this case - simulation doesn''t extend 100 s past start of burst')
            end
        end
        startFlux = curStart;
        maxIncrease = max(deriv);
        fprintf(fidSave,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n',auc,peak,numPks,ss,startFlux,maxIncrease,auc100,fluxArea);
    end
end
toc
fluxes = zeros(length(t),7);
for i = 1:length(t)
    [~,fluxes(i,:)] = neutrophil(t(i),y(i,:),p);
end
numPks(numPks==0) = .1;
startFlux(startFlux==0) = .01;
maxIncrease(maxIncrease==0) = 1e-12;
allOutputs = [auc,peak,numPks,ss,startFlux,maxIncrease];
if sobolTest
    fclose(fidSave);
    fclose(fid);
end

% end


%%
% first plots
figure
yyaxis left
hold on
plot(t,y(:,1),'-')
% ylim([0 .5])
yyaxis right
hold on
plot(t,y(:,2),'-')
ylim([0 210])
legend('c_i','c_{ER}')

%% linear regression for population analysis
xParamVar = xParamVar';
% allOutputs = allOutputs';
X = log(abs(xParamVar));
% X = X(:,[1:26,28:33]);
Y = log(abs(allOutputs));

%% Call the PLS routine
% PLS - nipals algorithm (2003)
[T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,R2x,R2y]=...
    PLS_nipals(X,Y,rank(X));

%% plot correlations
XNum = 18;
YNum = 2;
XCur = X(:,XNum);
YCur = Y(:,YNum);
figure
loglog(exp(XCur),exp(YCur),'*')
prettyGraph

%%
hold on
YBaseline = mean(YCur);
XBaseline = mean(XCur);
XTest = min(XCur):(max(XCur)-min(XCur))/100:max(XCur);
plot(exp(XTest),exp(YBaseline+(XTest-XBaseline)*Bpls(XNum,YNum)))

%% Sobol indices find NaN in output and store assoc parameter values
[outputFile,outputPath] = uigetfile('*.txt','Choose output file to examine');
[paramFile,paramPath] = uigetfile('*.txt','Choose parameter file');
[newSaveFile,newSavePath] = uiputfile('*.txt','Choose trouble param save file');
outputFile = fullfile(outputPath,outputFile);
paramFile = fullfile(paramPath,paramFile);
newSaveFile = fullfile(newSavePath,newSaveFile);
fid = fopen(outputFile,'r');
fidParam = fopen(paramFile,'r');
fidSave = fopen(newSaveFile,'w');
curOutputLine = fgetl(fid);
curParamLine = fgetl(fidParam);
lineNum = 1;
numNaN = 0;
while ischar(curOutputLine)
    try curOutput = sscanf(curOutputLine,'%f');
        if any(isnan(curOutput))
            numNaN = numNaN + 1;
            lineNum
            fprintf(fidSave,'%s\r\n',curParamLine);
        end
    catch
        break
    end
    lineNum = lineNum+1;
    curOutputLine = fgetl(fid);
    curParamLine = fgetl(fidParam);
end
fclose(fid);
fclose(fidParam);
fclose(fidSave);

%% Sobol indices replace NaN in output
[outputFile,outputPath] = uigetfile('*.txt','Choose output file to modify');
[resultsFile,resultsPath] = uigetfile('*.txt','Choose output file with revised results for NaN outputs');
[newSaveFile,newSavePath] = uiputfile('*.txt','Choose new save file');
outputFile = fullfile(outputPath,outputFile);
resultsFile = fullfile(resultsPath,resultsFile);
newSaveFile = fullfile(newSavePath,newSaveFile);
fid = fopen(outputFile,'r');
fidResults = fopen(resultsFile,'r');
fidSave = fopen(newSaveFile,'w');
curOutputLine = fgetl(fid);
lineNum = 1;
numNaN = 0;
while ischar(curOutputLine)
    try curOutput = sscanf(curOutputLine,'%f');
        if any(isnan(curOutput))
            curOutput(isnan(curOutput)) = 0;
            numNaN = numNaN + 1;
            curResultsLine = fgetl(fidResults);
            curOutput = sscanf(curResultsLine,'%f');
        end
        fprintf(fidSave,'%f\t%f\t%f\t%f\t%f\t%f\r\n',curOutput);
    catch
        break
    end
    lineNum = lineNum+1;
    curOutputLine = fgetl(fid);
end
fclose(fid);
fclose(fidResults);
fclose(fidSave);

%% Sobol indices graphs
[sobolFile,sobolPath] = uigetfile('*.txt','Choose file with saved Sobol indices');
sobolData = load(fullfile(sobolPath,sobolFile));
[sobolOrdered,orderIdx] = sort(sobolData(:,3),'descend');
keepIdx = orderIdx(1:5);
tickLabels = paramNames(orderIdx);
% tickLabels = cell(1,5);
% for i = 1:5
%     tickLabels{i} = num2str(keepIdx(i));
% end
% figure
xFirst = (1:5)-.15;
bar(xFirst,sobolData(keepIdx,1),.3,'b')
hold on
errorbar(xFirst,sobolData(keepIdx,1),sobolData(keepIdx,2),'k','LineStyle','none')

xTot = (1:5)+.15;
bar(xTot,sobolData(keepIdx,3),.3,'r')
errorbar(xTot,sobolData(keepIdx,3),sobolData(keepIdx,4),'k','LineStyle','none')
xticks([1 2 3 4 5])
xticklabels(tickLabels)
ylabel('Sobol index')
prettyGraph