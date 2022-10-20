Use of a Global Model to Understand Atmospheric Mercury Observations at Monitoring Sites in Latin America  
=========================================================================================================

Environmental pollution from Hg damages ecosystems through Hg's
transformation into toxic methylmercury and bio-accumulation in food
chains. Further, Hg is highly mobile in the atmosphere, allowing it to
travel to faraway places, resulting in worldwide distribution of its
elemental form, , which can last for as long as six months in the
atmosphere{cite}`horowitz_new_2017; @shah_improved_2021`. Hg in the
atmosphere can be classified as gaseous elemental Hg (GEM), gaseous
oxidized Hg (GOM), and particulate-bound Hg (PBM)
{cite}`lindberg_synthesis_2007; @schroeder_atmospheric_1998; @landis_development_2002`.
In most cases, Hg emissions occur as gaseous elemental , which is
relatively inert and sparingly soluble in water {cite}`horowitz_new_2017`.
Since most Hg entering ecosystems comes from the atmosphere, monitoring
and modeling atmospheric Hg and Hg deposition enables us to understand
its biogeochemical cycle. In addition, a better understanding of Hg's
circulation in the environment would enable effective policies to reduce
its harmful effects.

In this chapter, we compare the outputs of the GEOS-Chem model with
observed Hg at multiple sites in Latin America from two monitoring
networks. We combine simulations of Hg in the atmosphere produced by the
CTM (Sect. [1.2.2](#c2_geos_chem_simulations){reference-type="ref"
reference="c2_geos_chem_simulations"} - Sect.
[1.2.2](#c2_geos_chem_simulations){reference-type="ref"
reference="c2_geos_chem_simulations"}) with ground-based observations of
atmospheric total gaseous mercury (TGM) (Sect.
[1.2.3](#c2_monitoring_site_characteristics){reference-type="ref"
reference="c2_monitoring_site_characteristics"}) from the Global Mercury
Observation System (GMOS){cite}`sprovieri_atmospheric_2016` and gaseous
elemental mercury (GEM) data from a network of passive air samplers
(PAS){cite}`quant_measuring_2021` distributed across Latin America. Then, we
present comparisons between observations and model outputs and discuss
the implications of the current state of atmospheric monitoring and
modeling atmospheric Hg in Latin America(Sect.
[1.3](#c2_results){reference-type="ref" reference="c2_results"}).
Finally, we summarize conclusions from the analysis (Sect.
[1.4](#c2_conclusion){reference-type="ref" reference="c2_conclusion"}).

Background {#chapter2_background}
----------

Hg monitoring networks and atmospheric Hg models are closely
interconnected in the literature. For instance, Gustin et
al's.{cite}`gustin_mercury_2020` describe up-to-date scientific thinking
regarding mercury in the environment, monitoring and modeling
techniques, and how they relate to Minamata Convention on Mercury.
Notably, the mutual dependency between atmospheric modeling and
monitoring is presented in the MC's \"Monitoring Guidance
Document\"{cite}`unep_guidance_2021`, which states that observations are
needed not only to detect and quantify changes but also to improve and
evaluate models of mercury transport, fate, exposure, and
impacts{cite}`unep_guidance_2021`. Similarly, Sprovieri et
al.{cite}`sprovieri_atmospheric_2016` emphasize the importance of consistent
global Hg measurements to validate regional and global-scale models.
According to Brasseur and Jacob{cite}`brasseur_modeling_2017`, it is crucial
to have a large ensemble of observations to evaluate atmospheric
modeling outputs.

Data from passive and active sampling sites are evaluated in this
thesis. The MC's \"Monitoring Guidance Document\"{cite}`unep_guidance_2021`
provides a comprehensive definition and comparison of active and passive
air sampling. Active air samplers, particularly automated ones, can
deliver high-frequency data in a short period from as little as 5
seconds to 5 minutes{cite}`unep_guidance_2021`. However, they may be complex
and costly to set up and operate. In contrast to active samplers, PAS
are low-cost, inexpensive, and do not require electricity, moving parts,
pump operation, or calibration{cite}`macagnano_passive_2018`. Despite their
small size, they can be deployed at background, remote, urban, and
hotspot sites without worrying about media
failures{cite}`unep_guidance_2021`.

Some world regions, such as Latin America and Africa, with high ASGM Hg
emissions, have not been subject to detailed model-observation
comparisons. This may be attributed to these regions' lack of wide
coverage of required high-frequency atmospheric Hg monitoring capacity,
as shown in Figure
[1.1](#fig:global-hg-monitoring-networks){reference-type="ref"
reference="fig:global-hg-monitoring-networks"}, which illustrates the
distribution of different Hg monitoring networks
worldwide{cite}`united_nations_environment_programme_technical_2019]. It is
evident in Figure
[1.1](#fig:global-hg-monitoring-networks){reference-type="ref"
reference="fig:global-hg-monitoring-networks"} that Latin America,
Africa, and South East Asia remain significantly behind Europe and North
America regarding access to large observation ensembles. A majority of
the sampling sites present in these regions are PASs.

![Global map of Hg monitoring networks
{cite}`united_nations_environment_programme_technical_2019`](templates/figures/global-hg-monitoring-networks.pdf){#fig:global-hg-monitoring-networks
width="\\textwidth"}

Methods {#c2_methods}
-------

### GEOS-Chem Description {#c2_geos_chem_description}

The global atmospheric Hg concentration was simulated using version
12.8.1 of GEOS-Chem, whose Hg simulation is described by Horowitz et
al.{cite}`horowitz_new_2017`. All the simulations in this study were run
globally for 47 vertical layers at a resolution of 2.0$\times$2.5, which
is approximately equal to a 222 km$\times$277.5 km grid square at the
equator {cite}`horowitz_new_2017`. Moreover, the MERRA-2 assimilated
meteorological data {cite}`gelaro_modern-era_2017` drive the model's
atmospheric transport, which calculates atmospheric Hg from three
tracers: elemental Hg, Hg^0^, divalent Hg, Hg^2+^, and particulate-bound
divalent Hg, Hg^p^. The Hg chemical scheme in the GEOS-Chem version used
in this study considers bromine (Br) to be the primary
oxidant{cite}`horowitz_new_2017` and employs monthly mean Br oxidant
concentrations from Schmidt et al.{cite}`schmidt_modeling_2016`.

### GEOS-Chem Simulations {#c2_geos_chem_simulations}

The GMA 2018 emissions inventory was used to represent anthropogenic
emissions sources from all sectors{cite}`steenhuisen_development_2019`.
Different inputs to the GEOS-Chem model, such as emissions sources, can
be toggled on or off depending on the research objective; hence a
reference simulation, was created by turning on all Hg emissions sources
globally. Moreover, a was generated by turning off the ASGM source
globally to evaluate the contribution of ASGM to the baseline modeled
Hg$^0$ in the atmosphere by calculating the difference between the and .
Table
[\[tab:geos\_chem\_simulation\_description\]](#tab:geos_chem_simulation_description){reference-type="ref"
reference="tab:geos_chem_simulation_description"} describes the
simulations that were conducted in detail.

[\[tab:geos\_chem\_simulation\_description\]]{#tab:geos_chem_simulation_description
label="tab:geos_chem_simulation_description"}

The simulation frequency was set to output daily averages at the global
scale, while the output for the grid boxes corresponding to the
locations of the GMOS observation sites was set to an hourly frequency.
The GEOS-Chem outputs for all the simulations were in units of parts per
trillion (ppt) and were converted to at standard temperature and
pressure (273 K, 1 atm) to compare them to observations.

### Atmospheric Mercury Monitoring Sites in Latin America {#c2_monitoring_site_characteristics}

The GMOS network is one of a few major projects to develop a global
observing system for Hg pollution. GMOS aims to provide high-quality Hg
data sets in the Northern and Southern hemispheres to enable a more
comprehensive assessment of atmospheric Hg concentrations and their
dependence on meteorology, long-range atmospheric transport, and
atmospheric emissions{cite}`sprovieri_atmospheric_2016`. A vast network of
ground-based monitoring stations, regular oceanographic cruises, and
lower, upper, and stratospheric measurements make up this European
Union-funded project
{cite}`koenig_seasonal_2021; @sprovieri_atmospheric_2016`. More than 40
ground-based monitoring sites constitute the GMOS network, covering many
regions with limited to no observational data before
GMOS{cite}`sprovieri_atmospheric_2016`. The GMOS monitoring network has five
sites in Latin America that actively monitor Hg levels. A detailed
analysis of the Sisal, Calhau, Manaus, Nieuw Nickerie, and Bariloche
sites was conducted by Sprovieri et al.{cite}`sprovieri_atmospheric_2016`,
and the Chalcataya site was analyzed in detail by Koenig et
al.{cite}`koenig_seasonal_2021`. A summary of the sites' characteristics is
shown in Table
[\[tab:gmos\_sites\_info\]](#tab:gmos_sites_info){reference-type="ref"
reference="tab:gmos_sites_info"}. Moreover, the distribution of these
GMOS sites in Latin America is indicated by the red triangles in Figure
[1.2`(#fig:Latam_Passive_SamplerSites){reference-type="ref"
reference="fig:Latam_Passive_SamplerSites"}, which is a map showing the
names and locations of the GMOS Monitoring Network Sites and Passive
Sampler sites in Latin America
{cite}`quant_measuring_2021; @koenig_seasonal_2021`. Moreover, PAS data from
the Latin American Passive Air sampling Network (LAPAN), which was
analyzed in detail by Quant et al.{cite}`quant_measuring_2021`, was used for
the model observation comparison. The respective locations of the PAS
sites are shown by the blue circles in Figure
[1.2](#fig:Latam_Passive_SamplerSites){reference-type="ref"
reference="fig:Latam_Passive_SamplerSites"}.

\[Characteristics of the GMOS sites evaluated\]Characteristics of the
GMOS sites evaluated
{cite}`koenig_seasonal_2021; @sprovieri_atmospheric_2016`.
[\[tab:gmos\_sites\_info\]]{#tab:gmos_sites_info
label="tab:gmos_sites_info"}

The GMOS sites are classified as either secondary or master sites in
Table[\[tab:gmos\_sites\_info\]](#tab:gmos_sites_info){reference-type="ref"
reference="tab:gmos_sites_info"} to indicate the type of data collected
and the type of equipment used at the site. Master stations are those
where Gaseous Elemental Hg (GEM, i.e., the gas phase Hg in its ground
electronic state), Gaseous Oxidized Hg (GOM, i.e., the oxidized gas
phase Hg compounds), Hg associated with suspended particulate matter
(PBM2.5) and Hg in precipitation are continuously measured while
secondary stations are those where only GEM and Hg in precipitation are
continuously measured
{cite}`sprovieri_atmospheric_2016; @gustin_measuring_2015; @koenig_seasonal_2021`.

![Map showing the names and locations of the GMOS Monitoring Network
Sites and Passive Sampler sites in Latin America. GMOS sites are
indicated by the red triangles, and the PAS sites are indicated by the
blue
dots{cite}`quant_measuring_2021; @koenig_seasonal_2021`.](templates/figures/Passive_Samplers/Latam_Passive_SamplerSites.pdf){#fig:Latam_Passive_SamplerSites
width="80%"}

### Pre-processing and Comparison of Observed and Modeled Mercury Concentration in the Atmosphere {#c2_observation_data_manipulation}

Annual average GEM concentration data for 27 PAS sites in Latin America
was obtained from Quant et al.{cite}`quant_measuring_2021`, which included
information about the coordinates of the deployment sites and the period
of measurement. The PAS data was already hence there was no need for
pre-processing before comparison with the modeled Hg concentrations.
Furthermore, the PAS had been deployed for a year; the different
deployment dates ranged from October 20th, 2017, to March 14th, 2020.
Therefore the annual average GEM concentrations from the PASs were
compared to the modeled annual average for 2015 for each site.

Available Hg observation data from the GMOS stations on Figure
[1.2](#fig:Latam_Passive_SamplerSites){reference-type="ref"
reference="fig:Latam_Passive_SamplerSites"} was obtained from the GMOS
online database (http://www.gmos.eu), as well as published studies about
the Hg monitoring data from the different sites
{cite}`sprovieri_atmospheric_2016; @koenig_seasonal_2021`. The data sets were
pre-processed based on the information in Sprovieri et
al.{cite}`sprovieri_atmospheric_2016` and Koenig et
al.{cite}`koenig_seasonal_2021`. Daily and annual averages of the observed
TGM concentration were calculated to compare with the output.

Results and Discussion {#c2_results}
----------------------

![The average annual Hg concentration on Latin America's surface. The
background is the yearly average concentration generated by the Base
(ASGM=ON) simulation for 2015. Circles represent the annual average GEM
concentration at PAS sites, while triangles represent the yearly average
TGM concentration at GMOS
sites{cite}`sprovieri_atmospheric_2016; @quant_measuring_2021; @koenig_seasonal_2021`.](templates/figures/Passive_Samplers/07-27-22_pas_vs_model_Hg0-per-year_001.pdf){#fig:06-12-22_pas_vs_model_Hg0-per-year_001
width="75%"}

Recent publications analyzing global Hg monitoring data highlight an
observed inter-hemispheric gradient of Hg concentration where Hg
concentration in the southern hemisphere is lower than Hg concentration
in the northern
hemisphere{cite}`united_nations_environment_programme_technical_2019; @sprovieri_atmospheric_2016`.
The gradient is evident in the simulated background annual average
concentration as seen in Figure
[1.3](#fig:06-12-22_pas_vs_model_Hg0-per-year_001){reference-type="ref"
reference="fig:06-12-22_pas_vs_model_Hg0-per-year_001"}. Moreover, most
GMOS sites agree with and validate the modeled interhemispheric
gradient. However, a glance at Table
[\[tab:model\_percentage\_overestimation\_of\_mean\]](#tab:model_percentage_overestimation_of_mean){reference-type="ref"
reference="tab:model_percentage_overestimation_of_mean"} shows that the
model overestimates the annual average Hg concentration at all the GMOS
sites. Moreover, Figure
[1.4](#fig:gmos_sites_stats){reference-type="ref"
reference="fig:gmos_sites_stats"} better visualizes the difference
between the modeled and observed concentrations. The bar chart in Figure
[1.4](#fig:gmos_sites_stats){reference-type="ref"
reference="fig:gmos_sites_stats"} compares observed and modeled Hg
concentrations at the GMOS sites. Observed concentrations are indicated
in red, modeled in blue, and ASGM contribution in green. There are
annotations on the bars indicating the average concentration of Hg. In
addition, each bar is annotated above with the data set's standard
deviation and error bars.

\[Comparison of the modeled and the observed TGM concentrations at the
GMOS sites.\]Comparison of the modeled and the observed TGM
concentrations at the GMOS sites. The percentage difference between the
model predictions and the observations depicts the extent to which the
model predicts the observed TGM concentrations
[\[tab:model\_percentage\_overestimation\_of\_mean\]]{#tab:model_percentage_overestimation_of_mean
label="tab:model_percentage_overestimation_of_mean"}

![Bar chart comparing the modeled and observed average Hg concentration
at the respective GMOS Sites. The blue bars indicate the modeled annual
average concentration, the red bars indicate the observed annual average
concentration, and the green bars indicate the ASGM contribution at each
site. The bars are annotated with the average Hg concentration values.
Moreover, the error bars and the annotated value above the bars show the
standard deviation in each data set.
](templates/figures/GMOS_Sites/gmos_sites_stats.pdf){#fig:gmos_sites_stats
width="\\textwidth"}

Figure [1.4](#fig:gmos_sites_stats){reference-type="ref"
reference="fig:gmos_sites_stats"} shows that ASGMs modeled contribution
is low in most sites except for Chacaltaya, Manaus, and Nieuw Nickerie.
The model's behavior regarding the predicted ASGM contribution at these
sites is not surprising since these three sites are in countries
estimated to be among the top 10 Latin American ASGM Hg emitters in the
ASGM emission inventory used for simulation. Even though the model
estimates a notable ASGM contribution at the Manaus (18%) and Nieuw
Nickerie (16%) sites, the sites lack enough data to fully characterize
the ASGM contribution to the modeled concentration over the long term.
However, the predicted ASGM Hg contribution at Chacaltaya is the highest
at 23% as seen in Table
[\[tab:model\_percentage\_overestimation\_of\_mean\]](#tab:model_percentage_overestimation_of_mean){reference-type="ref"
reference="tab:model_percentage_overestimation_of_mean"}.

As far as the PAS observed GEM concentrations are concerned, the modeled
background concentrations in Figure
[1.3](#fig:06-12-22_pas_vs_model_Hg0-per-year_001){reference-type="ref"
reference="fig:06-12-22_pas_vs_model_Hg0-per-year_001"} seem to match
the PAS GEM measurement at Chalcataya, but, according to Quant et
al.{cite}`quant_measuring_2021`, the higher GEM levels observed at Chacaltaya
(1.4) are likely to reflect a known Hg spill near the sampling site and
may not reflect regional Hg concentration values. The Chacaltaya Hg
spill occurred after the GMOS Chacaltaya site had completed its data
collection period. Thus, the GMOS Chacaltaya station's annual average Hg
concentration reflects regional Hg concentration more than the PAS
station. In comparison with PAS data, the GEOS-Chem model also
overestimated atmospheric concentrations. This phenomenon is more
prevalent in inland sites than coastal ones.

In the Amazon region, for example, there is a difference between the
model and inland PAS sites. This may be because the GEOS-Chem model
version used in this study underestimates Hg uptake by plants
{cite}`feinberg_evaluating_2022`, which means that the model predicts a
higher Hg concentration in the atmosphere than it is. Figure
[1.5](#fig:06-12-22_pas_vs_model_Hg0-per-year_by-latitude_001){reference-type="ref"
reference="fig:06-12-22_pas_vs_model_Hg0-per-year_by-latitude_001"}
which shows the modeled (blue circles) and observed (red circles) annual
average plotted as a function of latitude indicates the interhemispheric
gradient observed in GEOS-Chem. The observation error bars represent the
replicate precision of the observations, while the model error bars
represent the bootstrap confidence interval for the mean annual .

![Hg Concentration in the atmosphere as a function of Latitude. The
(blue circles) and observed (red circles) annual average plotted are
plotted as a function of latitude to evaluate spatial trends across the
continent. The observation error bars represent the replicate precision
of the observations while the model error bars represent the 95^th^
bootstrap confidence interval for the mean annual
.](templates/figures/Passive_Samplers/06-12-22_pas_vs_model_Hg0-per-year_by-latitude_001.pdf){#fig:06-12-22_pas_vs_model_Hg0-per-year_by-latitude_001
width="\\textwidth"}

overestimation of Hg concentration in the Amazon region observed above
was also addressed in Feinberg et al.{cite}`feinberg_evaluating_2022` where
simulations were compared with litterfall, throughfall, and flux tower
measurements from 93 forested sites to evaluate vegetation as a Hg sink.
The study concluded that the version, 12.8.1 underestimates dry
deposition, which may explain why measurements of Hg concentration in
Latin America were lower than predicted by .

### Modeled vs. Observed Temporal Trends {#c2_modeled_vs_observed_trends}

![Time series plots of the observed TGM concentrations at different GMOS
sites in red with the corresponding modeled concentration in blue and
the associated ASGM contribution in green. Except for the CHC site,
where the data are from July 2014 and January 2016, the available data
and corresponding model outputs were plotted between January 2013 and
January
2016.](templates/figures/GMOS_Sites/GMOS_Sites.pdf){#fig:GMOSvsGC
width="\\textwidth"}

This study also compared observed and modeled data on a daily resolution
as seen in Figure [1.6](#fig:GMOSvsGC){reference-type="ref"
reference="fig:GMOSvsGC"}, which shows the time series of the modeled in
the atmosphere alongside the observed Hg and the simulated ASGM
contribution to the atmospheric Hg concentration at the GMOS sites. The
model version used in this study overestimated the concentration of Hg
on most days. However, the estimated average over the available
observation period was within one standard deviation of the observed Hg
in most of the sites except for the Manaus and Nieuw Nickerie sites.
overestimates the observed GEM concentrations at the Manaus and
Bariloche master sites by over 25% and the GEM concentration at Nieuw
Nickerie by 20%, which may be indicative of poor parameterization of GEM
in the model. Moreover, the overestimation of GEM concentrations in
Manaus further indicates the model's poor implementation of Hg plant
uptake through dry deposition, as discussed in Feinberg et
al.{cite}`feinberg_evaluating_2022`.

Figure [1.7](#fig:gmos_sites_scatter){reference-type="ref"
reference="fig:gmos_sites_scatter"} displays the scatter plots of
modeled Hg concentrations as a function of the observed Hg
concentrations for each site. Each plot uses the red line to evaluate
the linear relationship between the modeled and observed Hg
concentrations. Equations of the red regression line and the coefficient
of determination for each site are also displayed on the plots. The
general observation is that the model poorly matched the observations,
indicated by the mild and even flat slopes and shallow $R^2$ values. The
low correlations between the model and observations may be attributed to
poor vegetation uptake{cite}`feinberg_evaluating_2022`. Another plausible
hypothesis about the poor model prediction is that the input ASGM
emissions that GEOS-Chem uses are poorly parameterized. Wrong emissions
in the model may reduce the extent to which the model recreates the
observed atmospheric Hg concentrations. Chapter
[\[Chapter3\]](#Chapter3){reference-type="ref" reference="Chapter3"}
investigates this hypothesis in detail and highlights some causes, such
as the underestimation of emissions from high Hg-emitting regions like
Madre de Dios.

![Scatter plots of the modeled Hg concentration as a function of the
observed concentration. The red line is used to investigate the extent
of the linear relationship between the modeled and observed Hg
concentrations. The coefficient of determination($R^2$) and the equation
of the red regression line are shown above each site scatter
plot](templates/figures/GMOS_Sites/gmos_sites_scatter.pdf){#fig:gmos_sites_scatter
width="\\textwidth"}

### Comparison of Model Predictions

GEOS-Chem's skill in reproducing the Hg concentration measured using the
two sampling methods was analyzed using the scatter plots in Figure
[\[fig:pas\_vs\_gmos\]](#fig:pas_vs_gmos){reference-type="ref"
reference="fig:pas_vs_gmos"}. Figure
[\[fig:pas\_vs\_gmos\]](#fig:pas_vs_gmos){reference-type="ref"
reference="fig:pas_vs_gmos"} (a) shows the modeled annual mean for each
GMOS site as a function of the observed at the site, while Figure
[\[fig:pas\_vs\_gmos\]](#fig:pas_vs_gmos){reference-type="ref"
reference="fig:pas_vs_gmos"} (b) shows the modeled annual mean for each
of 6 PAS sites that are the closest to the 6 GMOS sites in (a) as a
function of the observed at the respective site. In each plot, the red
line is the regression line to investigate the strength of the
association between the modeled concentrations and observed
concentrations.

  -- --
     
  -- --

\[Scatter plots comparing model's predicted Hg Concentration means with
the GMOS and PAS average annual Hg Concentration means\]Scatter plots
comparing model's predicted Hg Concentration means with the GMOS and PAS
average annual Hg Concentration means. (a) shows modeled annual mean for
each site as a function of the observed at the site. The red line is the
regression line [\[fig:pas\_vs\_gmos\]]{#fig:pas_vs_gmos
label="fig:pas_vs_gmos"}

Figure [\[fig:pas\_vs\_gmos\]](#fig:pas_vs_gmos){reference-type="ref"
reference="fig:pas_vs_gmos"} shows that GEOS-Chem is better at
reproducing the actively monitored average Hg concentration (steep slope
and large $R^2$) than the data from the passive monitoring. Furthermore,
adding all the remaining Hg concentrations from the PAS did not improve
the slope or $R^2$ but worsened the relationship. Even though Figure
[1.7](#fig:gmos_sites_scatter){reference-type="ref"
reference="fig:gmos_sites_scatter"} informs us that GEOS-Chem poorly
predicted daily averages at the individual GMOS sites, Figure
[\[fig:pas\_vs\_gmos\]](#fig:pas_vs_gmos){reference-type="ref"
reference="fig:pas_vs_gmos"} shows that its predictions of actively
monitored Hg concentrations are better in aggregate than passively
monitored Hg concentrations. This result suggests that modeling studies
would gain better insights from comparisons with actively monitored Hg
concentrations. However, this result does not render PAS monitoring
obsolete. Using PAS technologies, we can obtain high-quality data about
Hg concentrations in the atmosphere and understand regional background
Hg concentrations over long periods. Since PAS networks are relatively
inexpensive and easy to install, countries can use them to understand
their Hg emissions better, providing valuable global and regional
monitoring data. A better understanding of temporal trends can be gained
from active monitoring, as evidenced by the analysis of GMOS TGM
concentrations. Furthermore, a single data set of actively monitored Hg
concentrations can be analyzed to generate metrics such as mean, , and ,
allowing multiple ways to compare modeled and observed concentrations.
Despite not recreating the exact Hg concentrations observed at GMOS
sites, improvements to the model, such as the update in the model's dry
deposition discussed in Feinberg et al.{cite}`feinberg_evaluating_2022` may
improve the predictions of measured Hg concentration in the atmosphere.
Additionally, Shah et al.'s{cite}`shah_improved_2021` improved mechanistic
model of the atmospheric redox chemistry of mercury may also reduce
GEOS-Chems's error in predicting the observed concentrations.

Conclusion {#c2_conclusion}
----------

This chapter explored the relationship between the modeled and observed
Hg concentrations in Latin America. The atmospheric Hg measurements from
six active monitoring sites across Latin America part of the GMOS
network were analyzed and compared with modeled concentrations at the
respective sites{cite}`koenig_seasonal_2021; @sprovieri_atmospheric_2016`.
Furthermore, annual average GEM measurements from 27 PAS sites across
Latin America were compared to 1-year beverages for concentrations at
the respective monitoring sites {cite}`quant_measuring_2021`. A relatively
weak relationship was found between the observed mercury species (GEM
and TGM) and those in the , demonstrating a need for improving the
models. Additionally, GEOS-Chem recreated the average Hg concentration
measurements from active monitoring stations better than PAS. Lastly,
improved vegetation uptake and more accurate emission parameterizations
may improve the model's prediction.
