# Parabolic Origami

This script generates a parabolic origami (reflector) by folding a paraboloid surface into petals. You can adjust its shape by editing the parameters in **[Parabolic_Origami0.m](./Parabolic_Origami0.m)**.

## How to modify the parameters

Open **[Parabolic_Origami0.m](./Parabolic_Origami0.m)** and locate the section marked `Parameter settings`. Change the values as needed:

```matlab
%% Parameter settings

a        = 0.005;            % Parabola coefficient: y = a * x^2
diameter = 100;              % Reflector diameter
ymax     = a * (diameter/2)^2; % Maximum height of the paraboloid surface
nY       = 20;               % Number of height samples
y        = linspace(0, ymax, nY); % Array of heights for horizontal cross‑sections
nPieces  = 20;               % Number of petals
```

- `a`: Controls the “steepness” of the parabola.  
- `diameter`: Sets the total width of the reflector.  
- `ymax`: Computed from `a` and `diameter`, this is the peak height of the paraboloid.  
- `nY`: Number of discrete slices (in height) to sample when building the surface.  
- `y`: A MATLAB vector of heights from 0 up to `ymax`.  
- `nPieces`: The number of fold‑out petals in the origami.

## Further reading

see **[Report.pdf](./Report.pdf)** for a detailed explanation of the underlying code and mathematics alogrithm.
