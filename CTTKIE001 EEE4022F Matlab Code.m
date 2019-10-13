cd '/Users/kierancattell/Documents/MATLAB/EEE4022F'
clear all; clc;

%Averaged data from results
%ave = [0.0153123842344230	0.0166184380830842	0.0164097032171806	0.0162681587550714	0.0135274313624527	0.00883529492523736	0.0142945522758480	0.0155482148773404	0.0119463353521801	0.0409042100219627	0.0426034887053390	0.0301563197090389	0.0553481912177105	0.0638929428945316	0.0545075175306614	0.0754565772961020	0.104701803298076	0.119489593511955	0.151980719608339	0.201754843899574	0.167875149496093	0.121232267965060	0.100285230175563	0.0536470085417544	0.0261018999823997	0.0499505182081873	0.0376987837531750	0.0102456459038541	0.0247309670118932	0.0369686689273559	0.0396415307713794	0.0417193999763622	0.0505161280734018	0.0552097440370425	0.0463861308896106	0.0288097584304199	0.0194083811605454];
%ave = [0.00675540305056303	0.0101960338913966	0.00980067444460197	0.0161637343516815	0.0283256140149799	0.0426543697468764	0.0419782509724377	0.0382850757226573	0.0359479900061052	0.0276860478122698	0.0173477995081473	0.0460008721230769	0.0564367570340506	0.0577133216592659	0.0702670770577873	0.104335799360713	0.131818741984749	0.167114650575728	0.200606260088036	0.159678040553382	0.130326229452697	0.104140699687178	0.0590295668626804	0.0648744890842076	0.0648130344581950	0.0337004684535039	0.0164486004224678	0.0394324587124635	0.0373812235606952	0.0277179062667996	0.0284650590382540	0.0354961535169523	0.0363752916045090	0.0339508814640990	0.0255413782414401	0.0247031365834713	0.0210324750970583];
%ave = [0.0135644555865789,0.0167164713618517,0.0269093336477559,0.0264285397031476,0.0258309881957451,0.0325501511769686,0.0394556781576875,0.0328343710987450,0.0265114755847834,0.0147100109273461,0.0104103930998753,0.0168418067103540,0.0109602702816955,0.0151353441574780,0.0524367896146303,0.0912384446455115,0.125125341875433,0.139834818719428,0.126167869183882,0.105854168180074,0.0922568007773474,0.0621836506075925,0.0335798141156509,0.0284195807670883,0.0472744103485140,0.0388124165420250,0.0273634080196487,0.0213534194922933,0.0181626492149375,0.0204882980891536,0.0279607577505592,0.0324018931734893,0.0288574099764431,0.0218228345797609,0.0115257937709983,0.00870095723445276,0.0125062474942421]];
%ave = [[0.00715819933607817,0.00806708131552422,0.0198533562981865,0.0167217346195905,0.0203082910458305,0.0241148794574872,0.0273035423928384,0.00625924113408726,0.0106102898456025,0.0204076361838597,0.0590773288348936,0.0910593199781340,0.0930185128763223,0.0709570179759627,0.0755123996040829,0.0734585317257516,0.0559461905507376,0.110870267429049,0.128428971798390,0.155852228613387,0.183034067180552,0.214090860801261,0.154994689728515,0.0959152227793916,0.0646633479234907,0.0226800611118996,0.0138456617748774,0.0274402764666611,0.0151155447926166,0.00570478149094250,0.0217893863583328,0.0283541364705240,0.0294301481310372,0.0316149939224665,0.0349043564189400,0.0265589560553121,0.0126019283878619]];

%% Parameters
f = 40e3; % Hz
lambda = 343/f;
T = 1/f;
omega = 2*pi*f;
beta = 2*pi/lambda;

d= 0.0098;
%d= 2*lambda

r = (0):lambda/1:2;
theta = 0:2*pi/3590:2*pi;
%theta= 0:deg2rad(5):2*pi;

ps = 120; %deg
phaseShift = deg2rad(ps); %radians

x=r.'*cos(theta);
y=r.'*sin(theta);

% Antenna locations
tx1y=0; tx1x=-3*d/2;
tx2y=0; tx2x=-d/2;
tx3y=0; tx3x=d/2;
tx4y=0; tx4x=3*d/2;



cd('Single')
for n = 0:36
    fileID = fopen("scope_"+string(n)+".csv",'r');
    dataArray = textscan(fileID, '%f%f%*s%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,2);
    fclose(fileID);
    data = [dataArray{1:end-1}];
    out(n+1) = 2*sqrt(2)*rms(data(:,2));
end
cd '/Users/kierancattell/Documents/MATLAB/EEE4022F'


out2=interp1(0:deg2rad(5):pi,out,theta);
out = out2;


%% Acoustic field
for ix=1:length(r)
    for iy=1:length(theta)
        r_tx1=sqrt( (x(ix,iy)-tx1x)^2 + (y(ix,iy)-tx1y)^2 );
        if theta(iy) <= pi
            E_tx1(ix,iy) = out(iy)*cos(omega - beta*r_tx1 + phaseShift*-3);
        else
            E_tx1(ix,iy) = 0;
        end
        
        r_tx2=sqrt( (x(ix,iy)-tx2x)^2 + (y(ix,iy)-tx2y)^2 );
        if theta(iy) <= pi
            E_tx2(ix,iy) = out(iy)*cos(omega - beta*r_tx2 + phaseShift*-2);
        else
            E_tx2(ix,iy) = 0;
        end
        
        r_tx3=sqrt( (x(ix,iy)-tx3x)^2 + (y(ix,iy)-tx3y)^2 );
        if theta(iy) <= pi
            E_tx3(ix,iy) = out(iy)*cos(omega - beta*r_tx3 + phaseShift*-1);
        else
            E_tx3(ix,iy) = 0;
        end
        
        r_tx4=sqrt( (x(ix,iy)-tx4x)^2 + (y(ix,iy)-tx4y)^2 );
        if theta(iy) <= pi
            E_tx4(ix,iy) = out(iy)*cos(omega - beta*r_tx4 + phaseShift*0);
        else
            E_tx4(ix,iy) = 0;
        end
    end
end

E_total=E_tx1+E_tx2+E_tx3+E_tx4;


figure('Color',[1 1 1]);


expected = E_total(end,:);
expected = abs(expected/max(abs(expected)));

polarplot2=polar(theta,expected); set(polarplot2,'LineWidth',2); hold on;
polarplot3=polar(theta(1:50:1801),expected(1:50:1801)); set(polarplot3,'LineWidth',2);
%polarplot4=polar(theta(1:50:1801),ave/max(ave)); set(polarplot4,'LineWidth',2);

% rad2deg(asin((ps/360)*(lambda/d))) + 90
% (find(expected == max(expected(:)))-1)/10
% (find(expected == max(expected(1:50:1801)))-1)/10
% (find(ave == max(ave))-1)*5
% xcorr(ave/max(ave),expected(1:50:1801),0,'coeff')
% 
% mean((ave/max(ave)-expected(1:50:1801)).^2)

ylim([0,1])
text(0,0,1,'')
text = "Theoretical vs Simulated Results "+string(ps)+"{\circ Phase Shift}";
title(text,'FontSize',20,'Position',[0,1.2,0])


%% Functions

function plotBeamShape(noOfReadings, startAngle, endAngle, amplitudes)
amplitudes = amplitudes/max(amplitudes);
%isNZ =(~amplitudes==0);
dA = (endAngle - startAngle)/(noOfReadings-1);
angles = startAngle:dA:endAngle;
disp = polar (deg2rad(angles), amplitudes);
set(disp,'LineWidth',2);

end

function out = getData(sourcefolder)
cd(sourcefolder)
for n = 0:36
    fileID = fopen("scope_"+string(n)+".csv",'r');
    dataArray = textscan(fileID, '%f%f%*s%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,2);
    fclose(fileID);
    data = [dataArray{1:end-1}];
    
    Fs = 5555555.5555
    L=1000;
    Y = fft(d5(2,1:1000)*11);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    f = Fs*(0:(L/2))/L;
    
    out(n+1) = P1(find(f == 40000));
    
end
cd '/Users/kierancattell/Documents/MATLAB/EEE4022F'
end
