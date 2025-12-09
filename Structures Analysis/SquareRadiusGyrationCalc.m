function radiusGyration = SquareRadiusGyrationCalc(sideLength, thickness)
shortSide = sideLength - 2*thickness;
moment = (sideLength^4 - shortSide^4) / 12;
crossArea = sideLength^2 - shortSide^2;
radiusGyration = (moment / crossArea) ^ (1/2);
end