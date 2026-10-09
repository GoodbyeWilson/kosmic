# Patient Dotplot

One panel per pathway: each donor's pathway score, as a point, with a box
plot per group (control and disease).

## Scores

The scores are the per-donor pathway scores from the last Pathway DE run,
on the scale of its scoring method. **Standardise per pathway (z-score)**
(default on) puts every pathway on a common scale and shares the y axis.

## The FDR above each panel

The value above each panel is the pathway's FDR from the Pathway DE
results, the same test as the DE workspace reports. A pathway without a
Pathway DE result has no value shown.
