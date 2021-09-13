% parse *.dat and *.par file for MultiFilament simulation for processing in
% MATLAB/Python
% started by Simon (sfs14@ic.ac.uk) on 14/03/18)
%
% assumes equal length swimmers
tic
set(0,'DefaultFigureWindowStyle','docked')

SimulationLocation = './';

SimulationName = 'videoCloudIsoBig';





plotDT = 1;
endDt  = 0;
Nworm = 20; %30;
Nsw =  120; %360; %30;
Nsw = 360;

% % 
% plotDT = 1;
% endDt  = 0;
% Nworm = 5; %30;
% % Nworm = 4; %30;
% Nsw = 600; 




% SimulationLocation = '../single_filament_RPY';
% SimulationName = 'output_data.dat';


%% parse *.par file and dimensions of data (based on Matlab's script generator)
parName = [SimulationLocation, SimulationName, '.par'];
delimiter = ' ';
startRow  = 2;
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
fileID = fopen(parName,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end

rawNumericColumns = raw(:, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]);
rawStringColumns = string(raw(:, 23));
PARAMETERS = table;
PARAMETERS.SimulationName = cell2mat(rawNumericColumns(:, 1));
PARAMETERS.SimulationTag = cell2mat(rawNumericColumns(:, 2));
PARAMETERS.Nsw = cell2mat(rawNumericColumns(:, 3));
PARAMETERS.Nworm = cell2mat(rawNumericColumns(:, 4));
PARAMETERS.a = cell2mat(rawNumericColumns(:, 5));
PARAMETERS.ah = cell2mat(rawNumericColumns(:, 6));
PARAMETERS.DL = cell2mat(rawNumericColumns(:, 7));
PARAMETERS.Sp4 = cell2mat(rawNumericColumns(:, 8));
PARAMETERS.Kap = cell2mat(rawNumericColumns(:, 9));
PARAMETERS.C = cell2mat(rawNumericColumns(:, 10));
PARAMETERS.MU = cell2mat(rawNumericColumns(:, 11));
PARAMETERS.curvature = cell2mat(rawNumericColumns(:, 12));
PARAMETERS.K0 = cell2mat(rawNumericColumns(:, 13));
PARAMETERS.TimeSteps = cell2mat(rawNumericColumns(:, 14));
PARAMETERS.StepsPerPeriod = cell2mat(rawNumericColumns(:, 15));
PARAMETERS.TOL = cell2mat(rawNumericColumns(:, 16));
PARAMETERS.gmres_tol = cell2mat(rawNumericColumns(:, 17));
PARAMETERS.broyden_maxiter = cell2mat(rawNumericColumns(:, 18));
PARAMETERS.gmres_maxiter = cell2mat(rawNumericColumns(:, 19));
PARAMETERS.plot_steps = cell2mat(rawNumericColumns(:, 20));
PARAMETERS.STRAINTWISTX = cell2mat(rawNumericColumns(:, 21));
PARAMETERS.STRAINTWISTY = cell2mat(rawNumericColumns(:, 22));
PARAMETERS.STRAINTWISTZ = cell2mat(rawNumericColumns(:, 23));
clearvars delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns;

%% parse *.dat file and save as arrays
datName = [SimulationLocation, SimulationName, '.dat'];
delimiter = ' ';
fileID = fopen(datName,'r');
formatPerBead = repmat('%s',1,7);
formatAllBeads = repmat(formatPerBead,1,PARAMETERS.Nsw*PARAMETERS.Nworm);
formatSpec = ['%s', formatAllBeads, '%[^\n\r]'];
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);

Nt = size(dataArray{1},1) - 1;
Np = PARAMETERS.Nsw*PARAMETERS.Nworm;

t = zeros(Nt,1);
X = zeros(Nt,Np);
Y = zeros(Nt,Np);
Z = zeros(Nt,Np);
qr = zeros(Nt,Np);
qi = zeros(Nt,Np);
qj = zeros(Nt,Np);
qk = zeros(Nt,Np);

t_in = dataArray{1};
for it = 1:Nt
    itp = it+1;
    t(it)  = str2num(t_in{itp});
end

Ncol = size(dataArray,2) - 2;
assert(Ncol == 7*Np)

BeadCount = 2;
ib = 1;
while (BeadCount < Ncol)
    for it = 1:Nt
        itp = it + 1;
        tmp = dataArray{BeadCount};
        X(it,ib)  = str2num(tmp{itp});

        tmp = dataArray{BeadCount+1};
        Y(it,ib)  = str2num(tmp{itp});

        tmp = dataArray{BeadCount+2};
        Z(it,ib)  = str2num(tmp{itp});

        tmp = dataArray{BeadCount+3};
        qr(it,ib) = str2num(tmp{itp});

        tmp = dataArray{BeadCount+4};
        qi(it,ib) = str2num(tmp{itp});

        tmp = dataArray{BeadCount+5};
        qj(it,ib) = str2num(tmp{itp});

        tmp = dataArray{BeadCount+6};
        qk(it,ib) = str2num(tmp{itp});
    end
    
    BeadCount = BeadCount + 7;
    ib = ib + 1;
end

clearvars delimiter fileID formatPerBead formatAllBeads formatSpec %dataArray

toc

% plotDT = 900;%200;%300;%200;%300;%50;%5;
% plotDT = 600;


% plotDT = 100;
% Nworm = 100;
% Nsw = 2;

%% change units
X = X./(Nworm*2.2);
Y = Y./(Nworm*2.2);
Z = Z./(Nworm*2.2);
%%

%% video
% writerObj = VideoWriter('movieVeryBigSim2.avi');
writerObj = VideoWriter('movieVeryBigSim3.avi');
writerObj.FrameRate=40;
open(writerObj);
%% 

% diagnostic plot
figure
hold on
t = 1;
% % plot3(X(1,:),Y(1,:),Z(1,:),'.','LineWidth',3,'markersize',1,'color',[.7 .7 .7])
% % plot3(X(1,:),Y(1,:),Z(1,:),'k-','LineWidth',1,'markersize',5)
% for nsw=0:(Nsw-1)
%         plot3(X(t,Nworm*nsw+1:Nworm*(nsw+1)),Y(t,Nworm*nsw+1:Nworm*(nsw+1)),Z(t,Nworm*nsw+1:Nworm*(nsw+1)),'-','LineWidth',.5,'markersize',1,'color',[t 0 t]./(Nt+10))
% end

% % plot3(X(t,:),Y(t,:),Z(t,:),'b-','LineWidth',.5)
cc = 1;
for t = 1:plotDT:((Nt-1))
    
    disp('progress: ')
    disp(t/Nt)
    
    
    hold off;
%     plot3(X(t,:),Y(t,:),Z(t,:),'.','LineWidth',.5,'markersize',1,'color',[t/Nt 0 t]./(Nt+10))
% % %     plot3(X(t,:),Y(t,:),Z(t,:),'-','LineWidth',.5,'color',[t 0 t]./(Nt+10))
    for nsw=0:(Nsw-1)
        plot3(X(t,Nworm*nsw+1:Nworm*(nsw+1)),Y(t,Nworm*nsw+1:Nworm*(nsw+1)),Z(t,Nworm*nsw+1:Nworm*(nsw+1)),'k-','LineWidth',.5,'markersize',1)
        hold on;
    end
    
%     for nsw=3
%         plot3(X(t,Nworm*nsw+1:Nworm*(nsw+1)),Y(t,Nworm*nsw+1:Nworm*(nsw+1)),Z(t,Nworm*nsw+1:Nworm*(nsw+1)),'b-','LineWidth',2)%,'markersize',1,'color',[t*cc t*(1-cc) .5*t]./(Nt+10))
%     end
% %     
%     for nsw=5
%         plot3(X(t,Nworm*nsw+1:Nworm*(nsw+1)),Y(t,Nworm*nsw+1:Nworm*(nsw+1)),Z(t,Nworm*nsw+1:Nworm*(nsw+1)),'r-','LineWidth',2)%,'markersize',1,'color',[t*cc t*(1-cc) .5*t]./(Nt+10))
%     end
    
    if (cc == 1)
        cc = 0;
    elseif (cc == 0)
        cc = 1;
    end
    
    xlim([min(min(X))-1 max(max(X))+1])
    ylim([min(min(Y))-1 max(max(Y))+1])
    zlim([min(min(Z))-1 max(max(Z))+1])

    axis normal
    axis equal
    axis on
    
    xlim([min(min(X))-1 max(max(X))+1])
    ylim([min(min(Y))-1 max(max(Y))+1])
    zlim([min(min(Z))-1 max(max(Z))+1])
    
%     % X(t,Nworm*nsw+1:Nworm*(nsw+1))
%     xlim([min(min(X(t,:)))-1 max(max(X(t,:)))+1])
%     ylim([min(min(Y(t,:)))-1 max(max(Y(t,:)))+1])
%     zlim([min(min(Z(t,:)))-1 max(max(Z(t,:)))+1])


    

    view([0 90])
%     view([-85 70]);
    view([-2.8000 7.6000])
    drawnow
    pause(.1)
    
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    
%     F = getframe(ax,[10 10 300 680]);

    F = getframe(gcf);

    writeVideo(writerObj, F.cdata);
    
end


close(writerObj);


% t = Nt-endDt;
% % plot3(X(t,:),Y(t,:),Z(t,:),'.','LineWidth',.5,'markersize',1,'color',[t 0 t]./(Nt+10))
% for nsw=0:(Nsw-1)
%     plot3(X(t,Nworm*nsw+1:Nworm*(nsw+1)),Y(t,Nworm*nsw+1:Nworm*(nsw+1)),Z(t,Nworm*nsw+1:Nworm*(nsw+1)),'-','LineWidth',.5,'markersize',1,'color',[t*cc t*(1-cc) .5*t]./(Nt+10))
% end
% plot3(X(t,:),Y(t,:),Z(t,:),'b-','LineWidth',.5)

% view([90 -90])
% view([10 10])
% view([0 0])
% view([-90 90])
% view([-85 70])
% view([0 0])
% grid on
% axis normal
% axis equal
% 
% view([0 0])
% 
% 
% view([10 80])
% view([-90 90])
% 
% xlim([min(min(X))-1 max(max(X))+1])
% ylim([min(min(Y))-1 max(max(Y))+1])
% zlim([min(min(Z))-1 max(max(Z))+1])
% 
% grid off
% box on
% 
% set(gcf,'renderer','Painters')
% 
% 
% title("$B = 10^2$")
% 
% 
% set(gca,'FontName','Times','FontSize',20);%,24);
% xlim([-5 15])
% zlim([-5 15])
% ylim([-150 14])


% for kk = 1:Np
%     plot3(X(1:plotDT:Nt,kk),Y(1:plotDT:Nt,kk),Z(1:plotDT:Nt,kk),'k-','Marker','.')
%     for t = 1:plotDT:Nt
%        R  = quaternion_rotation_matrix([qr(t,kk) qi(t,kk) qj(t,kk) qk(t,kk)]);
%        e1 = R(:,1);
%        e2 = R(:,2);
%        e3 = R(:,3);
%        quiver3(X(t,kk),Y(t,kk),Z(t,kk), e1(1), e1(2),e1(3),.2,'k-');
%        quiver3(X(t,kk),Y(t,kk),Z(t,kk), e2(1), e2(2),e2(3),.2,'r-');
%        quiver3(X(t,kk),Y(t,kk),Z(t,kk), e3(1), e3(2),e3(3),.2,'b-');
%        e3;
%        norm(e3)
%     end
% end

% 
% % diagnostic plot
% figure
% hold on
% % plot3(X(1,1),Y(1,1),Z(1,1:end-1),'r-')
% % plot3(X(end,:),Y(end,:),Z(end,:),'r-')
% for kk = 1:(Np-1)
%     plot3(X(1:plotDT:Nt,kk),Y(1:plotDT:Nt,kk),Z(1:plotDT:Nt,kk),'k-','Marker','.')
%     for t = 1:plotDT:Nt
%        R  = quaternion_rotation_matrix([qr(t,kk) qi(t,kk) qj(t,kk) qk(t,kk)]);
%        e1 = R(:,1);
%        e2 = R(:,2);
%        e3 = R(:,3);
%        quiver3(X(t,kk),Y(t,kk),Z(t,kk), e1(1), e1(2),e1(3),.2,'k-');
%        quiver3(X(t,kk),Y(t,kk),Z(t,kk), e2(1), e2(2),e2(3),.2,'r-');
%        quiver3(X(t,kk),Y(t,kk),Z(t,kk), e3(1), e3(2),e3(3),.2,'b-');
%     end
% end
% axis equal
% 
% 
% % diagnostic plot
% figure
% hold on
% % plot3(X(1,1),Y(1,1),Z(1,1:end-1),'r-')
% % plot3(X(end,:),Y(end,:),Z(end,:),'r-')
% for kk = 1:(Np-1)
%     plot3(X(1:plotDT:Nt,kk),Y(1:plotDT:Nt,kk),Z(1:plotDT:Nt,kk),'k-','Marker','.')
%     for t = 1:plotDT:Nt
%        R  = quaternion_rotation_matrix([qr(t,kk) qi(t,kk) qj(t,kk) qk(t,kk)]);
%        e1 = R(:,1);
%        e2 = R(:,2);
%        e3 = R(:,3);
%        
%        plot(t,e1(1)-1,'ko');
%        plot(t,e1(2),'k.');
%        plot(t,e1(3),'k+');
%        
%        plot(t,e2(1),'ro');
%        plot(t,e2(2)-1,'r.');
%        plot(t,e2(3),'r+');
%        
%        plot(t,e3(1),'bo');
%        plot(t,e3(2),'b.');
%        plot(t,e3(3)-1,'b+');
%        
%     end
% end
% 
% 
% % diagnostic plot
% figure
% hold on
% plot(qr(1,:).^2 + qi(1,:).^2 + qj(1,:).^2 + qk(1,:).^2 - 1,'r-')
% plot(qr(end,:).^2 + qi(end,:).^2 + qj(end,:).^2 + qk(end,:).^2 - 1,'r-')
% for kk = 1:plotDT:Nt
%     plot(qr(kk,:).^2 + qi(kk,:).^2 + qj(kk,:).^2 + qk(kk,:).^2 - 1,'k.')
% end
% 
% 
% % 3 bead displacement plot
% figure
% hold on
% plot(X(:,1),'ko')
% plot(X(:,2)-2.2,'r.')
% plot(X(:,3)-4.4,'gs')
% title('X')
% 
% figure
% hold on
% plot(Y(:,1),'ko')
% plot(Y(:,2),'r.')
% plot(Y(:,3),'gs')
% title('Y')
% 
% figure
% hold on
% plot(Z(:,1),'ko')
% plot(Z(:,2),'r.')
% plot(Z(:,3),'gs')
% title('Z')
% 
% %% plot angle between subsequent "normal" vectors
% figure
% hold on
% % plot3(X(1,1),Y(1,1),Z(1,1:end-1),'r-')
% % plot3(X(end,:),Y(end,:),Z(end,:),'r-')
% kk = 1;
% t  = 1;
% R  = quaternion_rotation_matrix([qr(t,kk) qi(t,kk) qj(t,kk) qk(t,kk)]);
% e1 = R(:,1);
% e2 = R(:,2);
% e3 = R(:,3);
% for kk = 1
%     for t = 2:Nt
%        R  = quaternion_rotation_matrix([qr(t,kk) qi(t,kk) qj(t,kk) qk(t,kk)]);
%        e1 = R(:,1);
%        e2new = R(:,2);
%        e3 = R(:,3);
%        
%        plot(t, acos(e2new*e2'),'ko');
%        e2 = e2new;
% %        quiver3(X(t,kk),Y(t,kk),Z(t,kk), e1(1), e1(2),e1(3),.2,'k-');
% %        quiver3(X(t,kk),Y(t,kk),Z(t,kk), e2(1), e2(2),e2(3),.2,'r-');
% %        quiver3(X(t,kk),Y(t,kk),Z(t,kk), e3(1), e3(2),e3(3),.2,'b-');
%     end
% end