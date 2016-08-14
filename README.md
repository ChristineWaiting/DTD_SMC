# Estimating Distance-to-Default with a Sector-Specific Liability Adjustment via Sequential Monte Carlo
Name of QuantLet : DTD_SMC

Published in : Applied Quantitative Finance third Edition

Description : Distance-to-Default (DTD), a widely adopted corporate default predictor, arises
from the classical structural credit risk model of Merton (1974). The modern way of
estimating DTD applies the model on an observed time series of equity values along
with the default point definition made popular by the commercial KMV model. It is
meant to be a default trigger level one year from the evaluation time, and is assumed to
be the short-term debt plus 50% of the long-term debt. This default point assumption,
however, leaves out other corporate liabilities, which can be substantial and particularly
so for financial firms. Duan, et al (2012) rectified it by adding other liabilities after
applying an unknown but estimable haircut. Typical DTD estimation uses a one-year
long daily time series. With at most four quarterly balance sheets, the estimated haircut
is bound to be highly unstable. Post-estimation averaging of the haircuts being applied
to a sector of firms is thus sensible for practical applications. Instead of relying on
post-estimation averaging, we assume a common haircut for all firms in a sector and
devise a novel density-tempered expanding-data sequential Monte Carlo method to
jointly estimate this common and other firm-specific parameters. Joint estimation is
challenging due to a large number of parameters, but the benefits are manifold, for
example, rigorous statistical inference on the common parameter becomes possible
and estimates for asset correlations are a by-product. Four industry groups of US firms
in 2009 and 2014 are used to demonstrate this estimation method. Our results suggest
that this haircut is materially important, and varies over time and across industries;
for example, the estimates are 77.85% in 2009 and 52.99% in 2014 for 40 randomly
selected insurance firms, and 0.94% for all 31 engineering & construction and 79.66%
for 40 randomly selected banks in 2014.

Keywords : 

See also : SFEacfar1, SFEacfar2, SFEacfma1, SFEacfma2, SFEpacfar2, SFEpacfma2

Author : Jin-Chuan Duan and Christine Wei-Ting Wang

Submitted :

Datafiles: Y2009_Sec20018.mat, Y2014_Sec20018.mat, Y2009_Sec20082.mat, Y2014_Sec20082.mat, 
Y2009_Sec20051.mat, Y2014_Sec20051.mat, Y2009_Sec20055.mat, Y2014_Sec20055.mat

Example: 
