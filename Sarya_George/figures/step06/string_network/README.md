# Three-halo network figures

All original selected proteins and both anchors remain in the seed views, including disconnected nodes. Added partners are distinguished by dashed centre outlines. A dagger preserves a supplied pathway-membership review flag.

Ring order: inner = stk/anc; middle = stp/anc; outer = stp/stk. Ring width is constant. Colour is the signed protein interaction effect, with a common -2 to +2 scale. Effects outside that range are colour-clipped only; original values are preserved in nodes.csv. Grey means unavailable/unsuitable quantitative estimate, not zero. Black ring dots indicate protein FDR<0.1.

Seed and full overview figures use the same deterministic expanded-network coordinates, exported in reports/step06/string_network/network_plot_coordinates.csv. The anchor-detail and four module zooms are separately laid out for legibility. Module zooms contain that module's original seeds plus both anchors and retained added partners adjacent to a module seed; they omit other module seeds, so the complete network is always the overview. Membership is exported in module_plot_node_membership.csv. Coordinates, proximity, and components are drawing aids rather than evidence of biological distance or direction.

Node centre colours use the supplied pathway scope; they are not assertions that every membership is canonical. STRING confidence controls edges, whereas protein FDR and effect describe the experiment. An edge is not direct phosphorylation evidence. No NES was used to select or colour proteins.

The source script uses base-R force layout/component packing and ggplot2 annular polygons, avoiding extra graph-package dependencies.
