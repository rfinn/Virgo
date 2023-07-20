#!/usr/bin/env python


"""
GOAL:
* subtract the continuum using color-correction from legacy images


PROCEDURE:
* create a g-r image
   2.5 log10(flux_r/flux_g)

* reproject onto halpha image

* mask r-band image

* calculate rms in r-band outside the masked regions

* apply color correction to continuum-subtracted images
  
"""

#
