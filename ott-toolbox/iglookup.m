function [m,iscosine]=iglookup(order,index);
%lglookup.m --- gives the ig mode for a given index for the lg->ig, ig->lg
%               series of conversions.
%
% Usage:
%
% [m,iscosine] = iglookup(order,index);
%
% PACKAGE INFO

m=(order-2*(floor(index(:)/2)));
iscosine=(rem(index,2)==0);

m=reshape(m,size(index));