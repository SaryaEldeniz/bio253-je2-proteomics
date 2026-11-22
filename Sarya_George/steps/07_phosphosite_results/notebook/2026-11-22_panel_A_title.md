# 2026-11-22 - Poster panel A title

## Question
Publish the reviewed separate panel A with the shorter title requested by the students.

## What we did
- Added `panel_A_baseline_sites.svg` and its PNG export under `figures/baseline_vs_interaction/`.
- Replaced the long heading and old panel-A heading with `A | Baseline phosphosite differences`.
- Cropped the unused top margin; retained the vertical genotype legend inside QoxB and the parent-protein-abundance limitation.
- Linked both poster exports from the Step 7 README.

## Results and validation
- SVG viewBox: `0 22 864 316`; PNG: 3456 x 1264 pixels.
- The complete SVG plot group is XML-identical to the reviewed source. Rendered pixels below the heading match the source after the crop.
- Three phosphosites, plotted measurements, group means, detection counts and FDR annotations were preserved. No model, effect estimate or statistical threshold was recalculated.
- Visual inspection confirmed the new title, legend and footer are legible and not clipped.

## Problems
No figure issues found. Publication was prepared from the current GitHub main in a separate worktree because the original local checkout contains unrelated changes and diverged commits.

## AI use
Codex edited the SVG heading and crop, exported the PNG, verified XML and rendered pixels, and prepared this focused publication. A separate Codex reviewer independently checked the title and plot preservation.

## Next steps
Use the SVG for the poster. Keep the separate panel-B interaction statement when interpreting these baseline differences.
