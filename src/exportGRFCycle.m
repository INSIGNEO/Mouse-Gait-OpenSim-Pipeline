function FyNorm = exportGRFCycle( ...
    baseFilename,...
    cycleNumber,...
    f0,...
    f1,...
    Fx,Fy,Fz,...
    Mx,My,Mz,...
    R,T,...
    outputMotionFolder,...
    frequency)

ANALOG_RATIO = 8;

a0 = f0*ANALOG_RATIO ;
a1 = f1*ANALOG_RATIO;

FxCycle = Fx(a0:a1);
FyCycle = Fy(a0:a1);
FzCycle = Fz(a0:a1);

MxCycle = Mx(a0:a1)/1000;
MyCycle = My(a0:a1)/1000;
MzCycle = Mz(a0:a1)/1000;

if max(abs(FzCycle)) < 0.1

    FyNorm = [];

    return

end

%% Rotate forces

F = [FxCycle';FyCycle';FzCycle'];

Frot = R*F;

FxCycle = Frot(1,:)';
FyCycle = Frot(2,:)';
FzCycle = Frot(3,:)';

%% Rotate moments

M = [MxCycle';MyCycle';MzCycle'];

Mrot = R*M;

MxCycle = Mrot(1,:)';
MyCycle = Mrot(2,:)';
MzCycle = Mrot(3,:)';

%% CoP

[CoPx,CoPy,CoPz] = ...
    computeCOP( ...
    FxCycle,FyCycle,FzCycle,...
    MxCycle,MyCycle,MzCycle);

CoP = [CoPx CoPy CoPz];

CoP = CoP*R' + T'/1000;

%% Normalize

FxNorm = normalizeCycle(FxCycle,960);
FyNorm = normalizeCycle(FyCycle,960);
FzNorm = normalizeCycle(FzCycle,960);

MxNorm = normalizeCycle(MxCycle,960);
MyNorm = normalizeCycle(MyCycle,960);
MzNorm = normalizeCycle(MzCycle,960);

CoPxNorm = normalizeCycle(CoP(:,1),960);
CoPyNorm = normalizeCycle(CoP(:,2),960);
CoPzNorm = normalizeCycle(CoP(:,3),960);

%% Export

motFile = sprintf( ...
    '%s_StanceNorm_%02d_grf.mot',...
    baseFilename,...
    cycleNumber);

xmlFile = sprintf( ...
    '%s_StanceNorm_%02d_grf.xml',...
    baseFilename,...
    cycleNumber);

exportGRFMot( ...
    motFile,...
    FxNorm,FyNorm,FzNorm,...
    CoPxNorm,CoPyNorm,CoPzNorm,...
    MxNorm,MyNorm,MzNorm,...
    outputMotionFolder,...
    frequency,...
    a0,...
    a1);

exportGRFXML( ...
    xmlFile,...
    motFile,...
    outputMotionFolder);

FyNorm = FyNorm(:);

TName = sprintf( ...
    '%s_StanceNorm_%02d_TranslationMatrix.mat',...
    baseFilename,...
    cycleNumber);

RName = sprintf( ...
    '%s_StanceNorm_%02d_RotationMatrix.mat',...
    baseFilename,...
    cycleNumber);


TNameFull = fullfile(outputMotionFolder,TName);
RNameFull = fullfile(outputMotionFolder,RName);

save(TNameFull, 'T');
save(RNameFull, 'R');

end
