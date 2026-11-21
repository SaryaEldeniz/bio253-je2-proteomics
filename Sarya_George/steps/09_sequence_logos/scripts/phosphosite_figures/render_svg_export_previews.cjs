// Visual QA only: render converted SVGs and compare them with PDF page previews.
// Usage: node render_svg_export_previews.cjs <svg_pages> <qa_folder> <sharp_module>
const fs = require('node:fs');
const path = require('node:path');
const [pageDir, qaDir, sharpModule] = process.argv.slice(2);
const sharp = require(sharpModule || 'sharp');
(async () => {
  const files = fs.readdirSync(pageDir).filter(x => x.endsWith('.svg')).sort();
  const atlases = new Map();
  for (const file of files) {
    const [, base, number] = file.match(/^(.*)_page(\d+)\.svg$/);
    const out = path.join(qaDir, file.replace('.svg', '.png'));
    await sharp(path.join(pageDir, file), { density: 96 }).flatten({ background: '#fff' }).png().toFile(out);
    if (!atlases.has(base)) atlases.set(base, []);
    atlases.get(base).push({ page: Number(number), svgPreview: out });
  }
  for (const [base, pages] of atlases) {
    const tileWidth = 600, tileHeight = 425;
    const layers = [];
    for (let row = 0; row < pages.length; row++) {
      const p = pages[row];
      for (const [column, file] of [path.join(qaDir, base + '_source-' + p.page + '.png'), p.svgPreview].entries()) {
        const input = await sharp(file).resize(tileWidth, tileHeight, { fit: 'contain', background: '#fff' }).png().toBuffer();
        layers.push({ input, left: column * tileWidth, top: row * tileHeight });
      }
    }
    await sharp({ create: { width: tileWidth * 2, height: tileHeight * pages.length, channels: 3, background: '#fff' } })
      .composite(layers).png().toFile(path.join(qaDir, base + '_PDF_left_SVG_right.png'));
  }
  console.log('Rendered ' + files.length + ' SVG pages and three PDF/SVG comparison sheets.');
})().catch(error => { console.error(error); process.exitCode = 1; });
