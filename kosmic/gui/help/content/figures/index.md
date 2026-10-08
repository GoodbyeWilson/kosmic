# Figures

Publication-quality matplotlib figure previews with export controls.

## Available figures

Figures are auto-generated based on available results. Gene DE produces volcano plots and gene heatmaps. Pathway DE (scoring mode) produces bubble plots, cascade plots, and dotplots. Enrichment (discovery mode) produces bar plots of enriched terms.

## Exporting

**Export figure…** saves the current figure as PNG, PDF or SVG; the format is set by the file type chosen in the save dialog. Use the size controls to adjust dimensions before export.

Exported figures have a white background in both the light and dark themes. Text in PDF and SVG files stays editable in a vector graphics editor. Scatter plots with 20,000 points or more, such as UMAPs, are embedded as an image inside PDF and SVG files so that the files remain small; axes, labels and other text stay vector. The resolution of PNG files and of these embedded images is the export DPI in Plot Settings (default 300).
