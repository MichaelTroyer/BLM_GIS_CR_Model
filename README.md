BLM CRM Model for Colorado

Author:
Michael Troyer

Date:
7/7/18


Summary:
Implments Chi-squared Goodness-of-Fit test for arbitrary input feature classes.
Tests each input feature class for significance and reweights the input classes
according to the standard residuals.

Computes site frequency across input analysis classes within surveyed areas as 'Actuals'.
Computes expected distribution under null hypothesis (no relationship) proportional
to the representation of each class within the surveyed parts of the analysis area.
Only calculates chi2 actuals and expected values for surveyed areas; the standard residuals
are then generalized to the study area in order to identify areas of lower/higher site potential. 
