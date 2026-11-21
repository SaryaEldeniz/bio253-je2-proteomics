// Render an exported draft for local visual QA. No scientific analyses are run.
// Usage: node render_foundation.mjs <artifact_tool.mjs> <input.pptx> <output_dir>
import fs from 'node:fs/promises';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
const [modulePath,input,outputDir]=process.argv.slice(2);
const {PresentationFile,FileBlob}=await import(pathToFileURL(modulePath).href);
await fs.mkdir(outputDir,{recursive:true});
const presentation=await PresentationFile.importPptx(await FileBlob.load(input));
const slide=presentation.slides.items[0];
const preview=await presentation.export({slide,format:'png',scale:0.5});
await fs.writeFile(path.join(outputDir,'JE2_A0_poster_foundation_preview.png'),new Uint8Array(await preview.arrayBuffer()));
const layout=await slide.export({format:'layout'});
await fs.writeFile(path.join(outputDir,'JE2_A0_poster_foundation.layout.json'),await layout.text());
console.log('Rendered the actual exported PPTX and saved its editable-layout inspection.');
