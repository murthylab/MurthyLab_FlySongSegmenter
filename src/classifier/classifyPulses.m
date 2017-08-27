function pulseLabels = classifyPulses(pulses, templates, classifier)
% classifies pulses into P1 and P2
% pulseLabels = classifyPulses(pulses, templates [OPTIONAL], classifier [OPTIONAl])
%
% ARGS
%  pulses     - matrix of normalized pulses (output of normalizePulses), Npulses x Nsamples
%  templates  - for P1 and P2 (2 x Nsamples) - OPTIONAL defaults to NM91 templates
%  classifier - OPTIONAL defaults to QDA classifier trained on NM91 pulses
% RETURNS
%  pulseLabels - 0: P1, 1: P2

if nargin==1 && exist('pulseTypeClassifier.mat', 'file') % load default classifier trained on NM91 data
   load('pulseTypeClassifier.mat')
   classifier = qda;
end
projectedPulses = pulses*templates';                     % project pulses onto templates
pulseLabels = classifier.predict(projectedPulses);       % classify projected data

