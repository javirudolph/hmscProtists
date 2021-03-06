
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hmscProtists

<!-- badges: start -->

<!-- badges: end -->

## Replicability statement

HMSC-Protists – Rudolph, Resetarits, Leibold

Intro: This project is motivated by two studies: In Resetarits et
al. 2018, we conducted a metacommunity experiment with protists in
microcosms to test the idea of “keystone communities” (Mouquet et
al. 2013). They keystone community concept posits that there may be
certain communities or patches that have a disproportionate effect on
the metacommunity. We did not find evidence for keystone communities but
we did find that different particular arrangements of patches (distinct
landscapes) had large effects on the results of variation partitioning
analysis (Cottenie 2005, Peres-Neto et al. 2006). This was somewhat
surprising since the landscapes were identical in many other respects
(e.g. number of patches, distribution of environmental values, and
actual connectivity network) and involved the same set of species. Part
of this variation was explained by spatial autocorrelation in
environmental conditions between neighboring communities.

Therefore, we conjecture that this result depends on the sensitivity of
patches to the particular identities of neighboring patches, especially
whether these patches were similar or different in environmental values.
We also conjectured that substantial residual variation should be
structured by species interactions (rather than being purely stochastic)
and that the interaction effects should be similar (in pattern even
though they should vary in average magnitude) across landscapes and
replicates.

Ovaskainen et al. (2017) propose that metacommunity data could be better
analyzed by a joint species distribution approach (jSDM) approach than
by variation partitioning (Peres-Neto et al. 2006) and provide a
particular package to do this called HMSC (hierarchical modeling of
species communities). This allows one to understand the contribution of
individual species on the overall pattern of variation in terms of
environmental influence, spatial effects, interactions with other
species and uses the mean value across species to calculate an overall
effect. In this case, the overall environmental and spatial effects are
akin (though technically distinct) from the effects identified in
variation partitioning, but in addition there is the determination of an
‘interaction’ effect that describes residual covariation among species.
In more recent work (Leibold et al. manuscript), we propose using HMSC
to analyze variation within metacommunities and outline an approach that
can compare these patterns across metacommunities. To do this we
modified HMSC to improve interpretation of variation components and we
implemented a variant that allows one to study the separate
contributions of individual patches to the overall metacommunity
structure. Other studies have also studied and justified HMSC in terms
of metacommunity analysis (Norberg et al. 2019, Rybicki et al. 2018).

Here we apply HMSC to the data from Resetarits et al. 2018 to address
the following questions and hypotheses.

Available data: (more details to be filled in) \* 10 individual
landscapes:  
\* Control metacommunity: (1 rep)  
\* Loss of peripheral (low connectivity) patch (2 separate reps) – one
rare habitat type, one common habitat type.  
\* Loss of central (high connectivity) patch (2 separate reps) – one
rare habitat type, one common habitat type.  
\* For the control metacommunities we also have temporal data (9-10
census at weekly intervals). For the treatments we only have the final
census.  
\* For each sample:  
\* Census of protists by abundance (a site-by-species matrix)  
\* For each patch, a set of environmental variables (a
site-by-environment matrix with two columns – habitat type (auto- vs
alloch- production), age of the patch since the last disturbance
(i.e. succession of food source).  
\* A spatial representation of the metacommunity (using Moran
eigenvector maps).  
\* HMSC as we implement in Leibold et al. (in prep) involves the
following steps:  
1\. conduct individual SDMs for each species and calculate variation
attributable to E (environment), S (spatial patterning), I -(covariation
with other species in the latent variables that involve residual
variation) and r^2 (amount of variation explained by the model).  
2\. perform the modification to address site-by-site contributions and
for each patch calculate the same variation components (E, S, I and r^2)
3. recalculate ‘relative E, S, I’ components – these are calculated so
they sum to 1  
4\. derive the species-by-species covariation in the I component – We
want to standardize these by the total I component so we can compare the
pattern of covariation across metacommunities even though they may vary
in total I variation.

Hypotheses and predictions:

Using Low connectivity metacommunities as pseudoreplicates  
1\. Sites that are connected to other sites with the same habitat will
contribute more to E than sites that are not. (spatial autocorrelation
related).  
2\. Species that are more specialized for a given habitat will
contribute more to E than sites that are not. (e.g. paramecium or
colpidium). To investigate this, we will measure niche specialization by
existing data on relative carrying capacities in monocultures (dark v
light) and regress that against the environmental component for each
species. \*we recognize that age of patch could also be an important
niche axis, but we are unsure at the moment of how to quantify this to
regress against environmental component for each species. We can
potentially create an index for each age by habitat where we rank
species in terms of abundance for that patch age by habitat.  
3\. Similarly to e, Species that are less specialized for a given
habitat will contribute more to I than species that are not. (i.e. found
in more patches and thus will interact more).  
4\. Species that have larger typical population size will contribute
more to E than species that have smaller typical population sizes. I.e.
larger populations are less stochastic. To investigate this, we would
regress carrying capacity of each species (mean of both habitat types)
in monoculture (D v C) against E component for that species. (This is
the absolute value of the carrying capacity – not relativized)  
5\. Sites that are more isolated from other sites of the same habitat
will contribute more to I than those that are more connected (highly
dependent on the species interactions already there). To investigate
this, we will find the absolute number of neighbors that are the same
habitat type (0-4) and regress that against the I and E component for
each site.  
6\. How does the age of the patch (i.e. time since last disturbance)
influence its importance for environment structuring? We expect that
older patches will contribute more to E and I and younger patches will
contribute more to S.  
7\. The matrix describing covariation among species in I will be similar
across all landscapes. To investigate this, we will utilize three
pseudoreplicates for each landscape…We will then run a model with aov(I
component \~ species + landscape)

Temporal Variability/time series  
1\. Overall patterns in E,S and I components will stabilize through time
during the experiment. We will test for this by testing for an asymptote
of the values through time using log regression  
2\. How does spatial autocorrelation contribute to time to stability?
(spatial correlation = clumping of habitat ¬patch types). Landscapes
with higher moran’s I (or autocorrelation) would reach stability
quicker.

Effect of metacommunity size  
3\. S and I should increase with decreasing size of the metacommunity
whereas E and r^2 will increase.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("javirudolph/hmscProtists")
```
