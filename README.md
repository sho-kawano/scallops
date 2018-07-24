---
title: "Scallops"
author: "Ricardo T. Lemos"
date: "7/21/2018"
output: html_document
---

# Description

Scallops is an R package for geostatistical modeling. 
It is inspired by the paper Lemos and Sanso' (2012).

## Installation

You can install the github development version via

```R
devtools::install_github('rtlemos/Scallops')
```

## Example

Let us load the famous "scallop" data and plot the georeferenced log-catches, 
along with a grid that we will use to interpolate those values.

```R
library(SemiPar)
library(ggmap)
data(scallop)
scallop$logcatch = log(scallop$tot.catch + 1)
dpc_grid = get_grid(c(-73.75, -71.25), c(38.5, 41), 0.5)

map.ny = get_map(location = c(-72.5, 39.75), zoom = 8)
ggmap(map.ny) + 
  geom_point(data = scallop, aes(x = longitude, y = latitude, size = logcatch), shape = 1) +
  geom_point(data = dpc_grid$coord, aes(x = lon, y = lat), shape = 3, color = 'red')
```

![](figs/domain.png)

We will be placing one Gaussian variable over each gridpoint. 
To have an idea of its influence on the interpolation, let us pick one in the
middle of the plot, (39.5N, -72.75E), and depict how its weight changes across space.

```R
get_influence_plot(dpc_grid, lat = 39.5, lon = -72.75) +
    geom_point(data = scallop, aes(x = longitude, y = latitude, size = logcatch), shape = 1)
```

![](figs/isotropic.png)

The plot above is for an "isotropic" kernel (actually, the kernel is not exactly isotropic,
because we are using a latitude-longitude coordinate system). 
By default, the kernel's range corresponds to twice the grid spacing.

Let us now fit the model and depict the spatial variation of the 
posterior weight for the same Gaussian variable.

```R
fit = get_mcmc(nburn = 10, nsample = 1000,
               s = data.frame(lon = scallop$longitude, lat = scallop$latitude), 
               dpc_grid = dpc_grid, y = scallop$logcatch, 
               priors = get_priors(dpc_grid = dpc_grid, precision_diagonal_value = 0.01),
               seed = 1)
get_influence_plot(dpc_grid, lat = 39.5, lon = -72.75, fit = fit) +
    geom_point(data = scallop, aes(x = longitude, y = latitude, size = logcatch), shape = 1)

```

![](figs/anisotropic.png)

We can also look at the shape of the kernels at all gridpoints.

```R
get_ellipses_plot(dpc_grid, fit)
```
![](figs/kernels.png)

And finally, let us look at the interpolation.

```R
get_interpolation_plot(obs_coord = scallop, dpc_grid = dpc_grid, fit = fit, contour_binwidth = 1) +
                       geom_text(data = scallop, aes(x = longitude, y = latitude, label = round(logcatch)))
```

![](figs/interpolation.png)

