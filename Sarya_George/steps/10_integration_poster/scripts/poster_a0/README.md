# A0 portrait poster foundation: source and assets

`build_poster_foundation.mjs` creates the native editable poster shapes and text in a presentation supplied by the Codex Artifact Session. It uses millimetres for layout, explicit PowerPoint points for fonts, and the user's provided institutional banner. It does not analyse the experimental measurements or invent results.

The current Methods section contains only the **Project analysis** SVG and a separate editable text slot. The experimental-design matrix and the complementary measurement-layer diagram are intentionally excluded from this poster.

- `build_project_analysis_svg.R` creates a poster adaptation of the existing `03_analysis_workflow.svg`. The original diagram is retained. This variant keeps the same analysis branches and selection logic, with larger labels and compact spacing. At its 529.33-mm placed width, the smallest label prints at approximately 25.7 pt.
- `build_poster_foundation.mjs` supplies the complete editable A0 layout, title/authors/affiliations/contact, background, hypothesis, Methods image/text slot, figure frames and reserved later sections.
- `render_foundation.mjs` imports the actual exported PPTX through the bundled presentation runtime and writes a PNG preview plus layout JSON for inspection. Its three command-line arguments are the runtime's `artifact_tool.mjs`, input PPTX and QA output directory. PDF copies can be exported from PowerPoint or LibreOffice at the original A0 page size.
- `assets/ETH_UZH_FGCZ_combined_banner.png` is an unchanged copy of the user-provided `ETHUZHFGCZ.png`. It is placed at 450 mm wide, preserving its aspect ratio, for approximately 170 dpi.

Recreate the SVG from `Sarya_George.Rproj`:

```r
source("steps/10_integration_poster/scripts/poster_a0/build_project_analysis_svg.R", encoding = "UTF-8")
```

In a Codex presentation Artifact Session, import `buildPosterFoundation`, read the PNG as a data URL and the SVG as UTF-8 text, then call:

```js
await artifactSession.run(async ({ presentation }) => {
  buildPosterFoundation(presentation, {
    logoDataUrl,
    analysisSvgText,
  });
  await presentation.export({
    format: "pptx",
    fileName: "JE2_A0_poster_foundation.pptx",
  });
}, { artifactType: "presentation" });
```

The poster exports and editing guide are in `../../presentations/poster_a0_foundation/`. The PowerPoint can be edited normally without the authoring runtime.
