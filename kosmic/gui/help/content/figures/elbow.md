# Elbow Plot

PCA variance explained per principal component. Used to inform the n_neighbors / n_pcs choice in clustering.

## Reading the plot

The grey bars give the percentage of variance each principal component
explains (left axis). The line gives the cumulative percentage (right
axis); dotted lines mark 80% and 90% when the components reach them.

The red vertical line marks the suggested number of components: the
point where the curve bends furthest from the straight line between its
first and last component. Components beyond it each add little variance.

The figure uses the width and height set in the controls.
