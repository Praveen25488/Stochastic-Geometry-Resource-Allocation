function [gLOS, gNLOS, gRIS] = sg_gain(model, N)
%SG_GAIN Average direct-path and RIS-assisted link gains.

gLOS = (1+model.r0)^(-model.alphaL);
gNLOS = (1+model.r0)^(-model.alphaN);
if N > 0
    gRIS = model.etaR*N^2/(1+model.r0/2)^(2*model.alphaL);
else
    gRIS = 0;
end
end
