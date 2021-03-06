### April sprint backlog

We will investigate how networks differ in space using Bs (species composition beta-diversity), Bos (interactions beta-diversity) and Bwn (metaweb beta-diversity) indexes. We will use site-level networks and coarse taxonomic resolution;

- [x] Add co-authors to repository and [Zotero library](https://www.zotero.org/groups/4659890/scalewebs-beta_div)
- [ ] Data wrangling - we need network matrices and geographical data for each network; @aammd to help us organize and navigate.

Site level food webs are in this folder: https://github.com/SrivastavaLab/empirical_food_webs/tree/main/03_food_web_metrics/data/pasted_matrices

Message from Diane: 
> "These site-level webs are then sampled at the bromeliad level using scripts found here: [empirical_food_webs\03_food_web_metrics](https://github.com/SrivastavaLab/empirical_food_webs/tree/main/03_food_web_metrics). Note that we don't save 1500 bromeliad food web matrices, instead we have functions (see [food_web_functions script](https://github.com/SrivastavaLab/empirical_food_webs/blob/main/03_food_web_metrics/food_web_functions.R)) that extract the bromeliad web, calculate a food web metric, and then return that metric. Also note that about 10% of cells in the site-level webs are not 1s or 0s, but fractional probabilities. We therefore assign these cells the value 1 (or 0) with that probability, calculate the food web metrics, and then repeat this many many times to get the mean value of the metric.
> In terms of location, climate and bromeliad environment information for site and bromeliad levels, you can find these files [here](https://github.com/SrivastavaLab/empirical_food_webs/tree/main/04_combine_metrics_climate/data)
> * [bromeliad_level_combined_data.csv](https://github.com/SrivastavaLab/empirical_food_webs/blob/main/04_combine_metrics_climate/data/bromeliad_level_combined_data.csv)
> * [speciesPool_combined_data.csv](https://github.com/SrivastavaLab/empirical_food_webs/blob/main/04_combine_metrics_climate/data/speciesPool_combined_data.csv)
* [Geographical information](https://github.com/SrivastavaLab/empirical_food_webs/04_combine_metrics_climate/data/visitSiteDict.csv)

- [x] Question feasibility of the project - which obstacles we might face in the next days?
  - Result from discussions in the working group: It is not possible to calculate beta-diversity of interactions because interactions are assumed to happen based on co-occurrence.
- [x] ~Beta-diversity Bs, Bos and Bwn calculations - refer to [this code](https://github.com/graciellehigino/the-return-of-the-fleas/blob/master/code.jl). @graciellehigino to lead the analysis.~
- [x] Creative figures time (1 hour) - explore the data however you want, come up with your three best figures. Maps, cluster analysis, gradients, relationships with [environmental data](https://github.com/SrivastavaLab/empirical_food_webs/tree/main/02_climate_data) are some ideas.
  - Some of the best figures are in the [working group report](https://htmlpreview.github.io/?https://github.com/sCaleWebs/wg2022-betadiv/blob/main/report.html).
- [x] Present and discuss your best figures. Group to select 3-5 figures.
- [ ] Write figures legends
- [ ] Paper outline based on figures
- [ ] Working plan moving forward (first draft deadline, submission deadline, journal, expectations and commitment)
- [ ] Populate Zotero with initial references
