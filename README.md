### US-hydraulic-geometry

The objective of this script is to conduct region-specific hydraulic geometry analyses by HUC2 watersheds.

**Empirical Datasets:**
[EPA NRSA homepage](https://www.epa.gov/national-aquatic-resource-surveys/nrsa)
- Wadeable Streams Assessment (2004-2005)
- National Streams and Rivers Assessment (2008-2009)

**Analysis steps:**
- Download field observations from EPA datasets (WSA and NRSA)
- Join field observations to nearest [NHDPlusV2 flowline segment](https://www.epa.gov/waterdata/nhdplus-national-hydrography-dataset-plus)
- Estimate upstream area for each sample location based on the location of the field sampling point along its corresponding NHDPlus flowline (i.e., proportional location along the reach flowline segment) and the total area at the bottom of the corresponding NHDPlus flowline.
- For each HUC2 region, conduct standardized major axis regression (SMA) to develop watershed area-width relationships.  


