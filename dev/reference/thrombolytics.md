# Thrombolytic treatments data

Data frame containing the results of 50 trials of 8 thrombolytic drugs
(streptokinase, SK; alteplase, t-PA; accelerated alteplase, Acc t-PA;
streptokinase plus alteplase, SK+tPA; reteplase, r-PA; tenocteplase,
TNK; urokinase, UK; anistreptilase, ASPAC) plus per-cutaneous
transluminal coronary angioplasty (PTCA) (Boland et al. 2003; Lu and
Ades 2006; Dias et al. 2011) . The number of deaths in 30 or 35 days
following acute myocardial infarction are recorded.

## Usage

``` r
thrombolytics
```

## Format

A data frame with 102 rows and 5 variables:

- studyn:

  numeric study ID

- trtn:

  numeric treatment code

- trtc:

  treatment name

- r:

  total number of events

- n:

  total number of individuals

## References

Boland A, Dundar Y, Bagust A, Haycox A, Hill R, Mota RM, Walley T,
Dickson R (2003). “Early thrombolysis for the treatment of acute
myocardial infarction: a systematic review and economic evaluation.”
*Health Technology Assessment*, **7**(15).
[doi:10.3310/hta7150](https://doi.org/10.3310/hta7150) .  
  
Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE (2011). “NICE
DSU Technical Support Document 4: Inconsistency in networks of evidence
based on randomised controlled trials.” National Institute for Health
and Care Excellence. <https://sheffield.ac.uk/nice-dsu>.  
  
Lu GB, Ades AE (2006). “Assessing evidence inconsistency in mixed
treatment comparisons.” *Journal of the American Statistical
Association*, **101**(474), 447–459.
[doi:10.1198/016214505000001302](https://doi.org/10.1198/016214505000001302)
.
