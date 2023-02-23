#pragma once

template <typename T> struct Parameters {
  // Domain
  T gridDelta = 0.2;
  T xExtent = 10.;
  T yExtent = 10.;

  // Geometry
  T trenchWidth = 4.;
  T trenchHeight = 8.;
  T taperAngle = -1.;

  // Process
  T processTime = 4.5;
  T sourcePower = 1.;
  T stickingProbability = 0.9;
};
