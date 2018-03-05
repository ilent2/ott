function [p,l]=lglookup(order,index);
%lglookup.m --- gives the lg mode for a given index for the lg->hg, hg->lg
%               series of conversions.
%
% Usage:
%
% [p,l] = lglookup(order,index);
%
% PACKAGE INFO

p=min([order-index(:).';index(:).']);
l=order-2*p;
l(index<round(order/2))=-l(index<round(order/2));

p=reshape(p,size(index));
l=reshape(l,size(index));