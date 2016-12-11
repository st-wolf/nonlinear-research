function [ u ] = potential( a, b, x )
u = a * (x .^ 4) + b * (x .^ 2);