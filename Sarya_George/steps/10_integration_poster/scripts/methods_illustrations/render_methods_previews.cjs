// Render only: original scientific artwork is authored by the adjacent R scripts.
// Run with Node and sharp available. Preview PNGs do not replace the vector SVGs.
const path = require('node:path');
const sharp = require('sharp');
const root = path.resolve(__dirname, '../../figures/methods_svg');
const names = ['00_complete_methods','01_experimental_design','02_proteomics_layers','03_analysis_workflow'];
(async () => {
  for (const name of names) {
    await sharp(path.join(root,name+'.svg'),{density:120})
      .flatten({background:'#ffffff'}).png().toFile(path.join(root,name+'.png'));
    console.log('Rendered '+name);
  }
})().catch(error => { console.error(error.message);process.exitCode=1; });
