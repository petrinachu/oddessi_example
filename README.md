ICTMC 2026 Poster -- example from ODDESSI trial

The non-collapsibility of odds and hazard ratios is well established (i.e. the ratio on population-averaged level from model adjusting for baseline covariates is not the same as the ratio conditional on model baseline covariates). If a population-level estimate of treatment effect is targeted, there are valid methods to produce marginal causal estimands such as standardisation [1-2]. Standardisation can be easily implemented in Stata using the 'margins' command. Presenting both absolute and relative treatment effects for binary outcomes is recommended [3], although this is not always done [4]. Using examples from the ODDESSI trial [5], this walks through how to generate marginal causal log-odd ratios and other marginal causal estimands.

References:
[1] Daniel R, Zhang J, Farewell D. Making apples from oranges: Comparing noncollapsible effect estimators and their standard errors after adjustment for different covariate sets. Biom J. 2021 Mar;63(3):528-557. 
[2] Morris T, Walker A, Williamson E et al. Planning a method for covariate adjustment in individually randomized trials: a practical guide. Trials. 2022; 23, 328
[3] Schulz KF, Altman D, Moher D, et al. CONSORT 2010 statement: updated guidelines for reporting parallel group randomized trials. BMJ. 2010;340:c2332.
[4] Rombach I, Knight R, Peckham N, et al. Current practice in analysing and reporting binary outcome data – a review of randomized controlled trial reports. BMC Med. 2020 Jun;18(1):147.
[5] Pilling S, Clarke K, Parker G, et al. Open Dialogue compared to treatment as usual for adults experiencing a mental health crisis: Protocol for the ODDESSI multi-site cluster randomized controlled trial. Contemp Clin Trials. 2022 Febl113:106664.
