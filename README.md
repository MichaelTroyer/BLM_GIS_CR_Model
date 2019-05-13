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


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
