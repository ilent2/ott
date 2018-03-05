function [m,n]=hglookup(order,index);
%hglookup.m --- gives the hg mode for a given index for the lg->hg, hg->lg
%               series of conversions.
%
% Usage:
%
% [m,n] = lglookup(order,index);
%
% PACKAGE INFO

n=index;
m=order-index;