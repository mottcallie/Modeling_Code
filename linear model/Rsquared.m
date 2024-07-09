% calculate r2 value
% Y - actual response array
% Yfit - fit based on predictor
%
function R2 = Rsquared(Y,YFit)
%% compute R2
YResid = Y - YFit;
SSresid=sum(YResid.^2,'omitnan');
SStotal=(length(Y)-1)*var(Y,'omitnan');

R2=1-SSresid/SStotal;

end