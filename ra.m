function [img,dims,scales,bpp,endian] = ra(fname)

% [img, dims,scales,bpp,endian] = RA(fname)
%
% calls read_avw

global homedir;

[img,dims,scales,bpp,endian] = read_avw(fname);