# MS2Decide
![alt text](https://github.com/MejriY/Decide_test/blob/master/image/ms2decide_logo.png?raw=true)

Mass spectrometry-based natural products targeted discovery oftenrelies on a complicated decision-making process involving tediouscomparison of exact masses data and tandemmassspectra-based annotation tools outputs against various spectral reference libraries.
Toaddress this bottleneck we present tandem mass spectrum to decision(MS2DECIDE) which leverages Decision Theory and expert knowledgeto aggregate the annotation outputs of GNPS, Sirius, and ISDB-Lotusand computes a recommendation for targeting natural products inregard to their potential novelty. We demonstrate, through two casestudies, that MS2DECIDE reliably captures the novelty of naturalproducts from their tandem mass spectra.

![alt text](https://github.com/MejriY/Decide_test/blob/master/image/all_article_workflow.png?raw=true)
# How we use tools
## Weighted iterative GNPS analog search (WIGAS)
![alt text](https://github.com/MejriY/Decide_test/blob/master/image/gnps_iterative.png?raw=true)

## Sirius
Sirius annotations were performed in batch mode by using Sirius 6. we utilized the \textit{ConfidenceScoreExact} measure.
Unfortunately, in some cases, Sirius was not able to propose an annotation. To remedy, we associated a value of 0.5 to Sirius matching score. Additionally, in scenarios where there is no match with GNPS or no match with Sirius, the tanimoto between GNPS and SIRIUS cannot be calculated. Hence, a default value of 0.25 was also assigned to $T_{gs}$ in these instances. 

## ISDB-Lotus
For ISDB-Lotus, since a strict library search was applied, we considered a zero answer as an important information regarding our definition of novelty. Hence, no mean value was associated.

# Validation

## Manufactured case
With our function in hand, we first opted to assess the ability of MS2Decide to recommend new NPs among a chemodiverse spectral library of 90 known compounds and 6 unreported ones.
![alt text](https://github.com/MejriY/Decide_test/blob/master/image/case_1.png?raw=true)

## Real case
In our second case study, MS2Decide was applied to assess its ability to streamline the discovery of new NPS from the LC-MS/MS analysis of the alkaloidic extract of a well known monoterpene indole alkaloids-producing plant, \textit{Pleiocarpa mutica} Benth.\ (Apocynaceae).
![alt text](https://github.com/MejriY/Decide_test/blob/master/image/case_2.png?raw=true)

[![PyPI - Version](https://img.shields.io/pypi/v/ms2decide.svg)](https://pypi.org/project/ms2decide)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/ms2decide.svg)](https://pypi.org/project/ms2decide)

-----

**Table of Contents**

- [Installation](#installation)
- [License](#license)

## Installation

```console
pip install ms2decide
```

## License

`ms2decide` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.
